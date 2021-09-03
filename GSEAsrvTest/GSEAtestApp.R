library(shiny)
library(shinythemes)
library(shinyjqui)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(dplyr)
library(DT)
library(enrichplot)
library(GSVA)
library(shinycssloaders)
library(ggplot2)
library(ggpubr)



shinytheme("sandstone")

ui <-
    
navbarPage("{GSEA Analysis}",
  tabPanel("GSEA",
           fluidPage(
               mainPanel(
                   h3("This is a test app for GSEA Analysis")
               ))),
  tabPanel("GSEA Analysis",
           fluidPage(
               title = "GSEA Analysis",
               sidebarPanel(
                   width = 2,
                   selectInput("geneSet", "Gene Set", choices = c("Gruber_Lab_MEP", "Gruber_Lab_HSC", "Gruber_Lab_CMP", "Gruber_Lab_MPP", "Gruber_Lab_GMP"))),
               mainPanel(
                   tabsetPanel(
                       id = "dataset",
                       tabPanel("Enriched Signatures Table", div(DT::dataTableOutput("enrich_sig_table"), style = "font-size:12px")),
                       tabPanel("Enrichment Plot", withSpinner(verbatimTextOutput("NESandPval"), type = 6),
                                withSpinner(plotOutput("enrichplot0", width = "100%", height = "800px"), type = 6),
                                div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:10px; height:500px; overflow-y: scroll"),
                                #withSpinner(verbatimTextOutput("LeadingGenes"), type = 6),
                                downloadButton("downloadGeneList", "Download Leading Edge Gene List")),
                       tabPanel("Gene Expression Heatmap", withSpinner(plotOutput("heatmap0", width = "100%", height = "2000px"), type = 6)))
                       ))),
  tabPanel("Single Sample GSEA Analysis",
           fluidPage(
             title = "Single Sample GSEA Analysis",
             sidebarPanel(
               selectInput("ssmethod", "Single Sample GSEA Enrichment Scoring Method", choices = c("ssgsea","gsva","zscore","plage")),
               selectInput("sstype", "Gene Set", choices = c("Gruber_Lab_MEP", "Gruber_Lab_HSC", "Gruber_Lab_CMP", "Gruber_Lab_MPP", "Gruber_Lab_GMP"))
#               fileInput("gmt.u", "Choose GMT File", accept = ".gmt"),
#               fileInput("expr.u", "Choose Expression Data File - tab delimited", accept = c(".tsv", ".txt"))
               ),
             mainPanel(
               tabsetPanel(
                 id = "dataset",
                    tabPanel("Single Sample Enrichment Output", div(DT::dataTableOutput("ssGSEAtab"), style = "font-size:10px")),
                    tabPanel("ssGSEA Heatmap", withSpinner(plotOutput("heatmap1", width = "100%", height = "800px"), type = 6)),
                    tabPanel("ssGSEA Boxplot", withSpinner(plotOutput("boxplot0", width = "500px", height = "500px"), type = 6)))
             ))))



# Gene Set Data
gmt <- read.gmt("~/R/GSEA/Data/Up_Reg_Myeloid_Cells_Gruber_Lab.gmt")
gmt$term <- gsub("MEP", "Gruber_Lab_MEP", gmt$term)
# Expression Data
expr_data <- read.table("~/R/GSEA/LIMMA/BM_Mouse_Model_TPM_20210816_reorder.max.txt", sep="\t", header=T, row.names=1, quote="")



server <- function(input, output, session) {
    output$enrich_sig_table <- DT::renderDataTable({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups
      groupB <- colnames(A)[grep("MIC", colnames(A))]
      ####----Signal-to-Noise Calculation----####
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ####----Reformatting----####
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ####----GSEA----####
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, eps = 0, pvalueCutoff = 1.0)
      gsea.df <- as_tibble(gsea.res@result)
        ## displaying the GSEA results as interactive data table
        DT::datatable(gsea.df,
                      extensions = c("KeyTable", "FixedHeader"),
                      caption = "Enriched Signatures",
                      options = list(keys = T, searchHighlight = T, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
            formatRound(columns = c(2:10), digits = 2)
    })
    
    output$NESandPval <- renderText({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups
      groupB <- colnames(A)[grep("MIC", colnames(A))]
      ####----Signal-to-Noise Calculation----####
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ####----Reformatting----####
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ####----GSEA----####
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, eps = 0, pvalueCutoff = 1.0)
      gsea.df <- as.data.frame(gsea.res@result)
      var <- input$geneSet
      if (var == "Gruber_Lab_MPP"){
        NES = gsea.df$NES[which(gsea.df[,'ID']=='Gruber_Lab_MPP')]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']=='Gruber_Lab_MPP')]
        GS = "Gruber_Lab_MPP"
      }
      if (var == "Gruber_Lab_HSC"){
        NES = gsea.df$NES[which(gsea.df[,'ID']=='Gruber_Lab_HSC')]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']=='Gruber_Lab_HSC')]
        GS = "Gruber_Lab_HSC"
      }
      if (var == "Gruber_Lab_GMP"){
        NES = gsea.df$NES[which(gsea.df[,'ID']=='Gruber_Lab_GMP')]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']=='Gruber_Lab_GMP')]
        GS = "Gruber_Lab_GMP"
      }
      if (var == "Gruber_Lab_CMP"){
        NES = gsea.df$NES[which(gsea.df[,'ID']=='Gruber_Lab_CMP')]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']=='Gruber_Lab_CMP')]
        GS = "Gruber_Lab_CMP"
      }
      if (var == "Gruber_Lab_MEP"){
        NES = gsea.df$NES[which(gsea.df[,'ID']=='Gruber_Lab_MEP')]
        Pval = gsea.df$pvalue[which(gsea.df[,'ID']=='Gruber_Lab_MEP')]
        GS = "Gruber_Lab_MEP"
      }
      NES.o <- paste("NES: ", NES)
      Pval.o <- paste("Pvalue: ", Pval)
      if (NES > 0){
        UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in the Acute Myeloid Leukaemia (AML) group.")
      }
      else {
        UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in the Acute Myeloid Leukaemia (AML) group.")
      }
      paste(NES.o, Pval.o, UpOrDown, sep = '\n')
    })
    
    output$enrichplot0 <- renderPlot({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups
      groupB <- colnames(A)[grep("MIC", colnames(A))]
      ####----Signal-to-Noise Calculation----####
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ####----Reformatting----####
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ####----GSEA----####
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, eps = 0, pvalueCutoff = 1.0)
      gsea.df <- as.data.frame(gsea.res@result)
      var <- input$geneSet
      if (var == "Gruber_Lab_HSC"){
          geneSetID = which(gsea.res[,"ID"] == "Gruber_Lab_HSC")
          title = gsea.res$Description[which(gsea.res[,"ID"] == "Gruber_Lab_HSC")]
      }
      if (var == "Gruber_Lab_CMP"){
        geneSetID = which(gsea.res[,"ID"] == "Gruber_Lab_CMP")
        title = gsea.res$Description[which(gsea.res[,"ID"] == "Gruber_Lab_CMP")]
      }
      if (var == "Gruber_Lab_GMP"){
          geneSetID = which(gsea.res[,"ID"] == "Gruber_Lab_GMP")
          title = gsea.res$Description[which(gsea.res[,"ID"] == "Gruber_Lab_GMP")]
      }
      if (var == "Gruber_Lab_MPP"){
          geneSetID = which(gsea.res[,"ID"] == "Gruber_Lab_MPP")
          title = gsea.res$Description[which(gsea.res[,"ID"] == "Gruber_Lab_MPP")]
      }
      if (var == "Gruber_Lab_MEP"){
          geneSetID = which(gsea.res[,"ID"] == "Gruber_Lab_MEP")
          title = gsea.res$Description[which(gsea.res[,"ID"] == "Gruber_Lab_MEP")]
      }
      gseaplot2(gsea.res,
          geneSetID, title, pvalue_table = F)
    })
    
    output$LeadingEdgeGenes <- DT::renderDataTable({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups
      groupB <- colnames(A)[grep("MIC", colnames(A))]
      ####----Signal-to-Noise Calculation----####
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ####----Reformatting----####
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ####----GSEA----####
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, eps = 0, pvalueCutoff = 1.0)
      gsea.df <- as.data.frame(gsea.res@result)
      GS <- input$geneSet
      ## Subset core enriched genes
      genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
      genes2 <- strsplit(genes1,"/")
      GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
      GeneSymbol$Rank <- rownames(GeneSymbol)
      GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
      DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
      #DT::datatable(GeneSymbol, options = list(lengthMenu = c(5, 10, 20, 100, 1000), pageLength = 10), selection=list(mode = "single", selected = c(1)))
    })
    
    output$downloadGeneList <- downloadHandler({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups
      groupB <- colnames(A)[grep("MIC", colnames(A))]
      ####----Signal-to-Noise Calculation----####
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ####----Reformatting----####
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ####----GSEA----####
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, eps = 0, pvalueCutoff = 1.0)
      gsea.df <- as.data.frame(gsea.res@result)
      GS <- input$geneSet
      ## Subset core enriched genes
      genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
      genes2 <- strsplit(genes1,"/")
      genes3 <- as.data.frame(genes2, col.names = "genes")
      GeneSymbol <- genes3$genes
      GeneSymbol$Rank <- rownames(GeneSymbol)
      GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
      filename = function() {
        paste(gene_symbol,"_",Sys.Date(),".csv",sep = "")
      }
      content = function(file) {
        write.csv(gene_symbol, file)
      }
    })
    
    output$heatmap0 <- renderPlot({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups
      groupB <- colnames(A)[grep("MIC", colnames(A))]
      ####----Signal-to-Noise Calculation----####
      A <- A + 0.00000001
      P = as.matrix(as.numeric(colnames(A) %in% groupA))
      n1 <- sum(P[,1])
      M1 <- A %*% P
      M1 <- M1/n1
      A2 <- A*A
      S1 <- A2 %*% P
      S1 <- S1/n1 - M1*M1 
      S1 <- sqrt(abs((n1/(n1-1)) * S1))
      P = as.matrix(as.numeric(colnames(A) %in% groupB))
      n2 <- sum(P[,1])
      M2 <- A %*% P
      M2 <- M2/n2
      A2 <- A*A
      S2 <- A2 %*% P
      S2 <- S2/n2 - M2*M2
      S2 <- sqrt(abs((n2/(n2-1)) * S2))
      rm(A2)
      # small sigma "fix" as used in GeneCluster
      S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
      S2 <- ifelse(S2 == 0, 0.2, S2)
      S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
      S1 <- ifelse(S1 == 0, 0.2, S1)
      M1 <- M1 - M2
      rm(M2)
      S1 <- S1 + S2
      rm(S2)
      s2n.matrix <- M1/S1
      ####----Reformatting----####
      s2n.df <- as.data.frame(s2n.matrix)
      s2n.df$GeneID <- rownames(s2n.df)
      rownames(s2n.df) <- NULL
      data <- dplyr::select(s2n.df, GeneID, V1)
      data.gsea <- data$V1
      names(data.gsea) <- as.character(data$GeneID)
      s2n.matrix.s <- sort(data.gsea, decreasing = T)
      ####----GSEA----####
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, eps = 0, pvalueCutoff = 1.0)
      gsea.df <- as_tibble(gsea.res@result)
      var <- input$geneSet
      if (var == "Gruber_Lab_HSC"){
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description=="Gruber_Lab_HSC"),"core_enrichment"])
      }
      if (var == "Gruber_Lab_CMP"){
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description=="Gruber_Lab_CMP"),"core_enrichment"])
      }
      if (var == "Gruber_Lab_GMP"){
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description=="Gruber_Lab_GMP"),"core_enrichment"])
      }
      if (var == "Gruber_Lab_MPP"){
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description=="Gruber_Lab_MPP"),"core_enrichment"])
      }
      if (var == "Gruber_Lab_MEP"){
        genes1 <- as.matrix(gsea.df[which(gsea.df$Description=="Gruber_Lab_MEP"),"core_enrichment"])
      }
      genes2 <- strsplit(genes1,"/")
      genes3 <- as.data.frame(genes2, col.names = "genes")
      gene_symbol <- genes3$genes
      ## convert expression matrix to numeric
      class(exp.mat) <- "numeric"
      ## Transforming data
      exp.mat1 = log2(exp.mat + 1) # log
      exp.mat2 = apply(exp.mat1, 1, scale); # z score
      exp.mat3 = apply(exp.mat2, 1, rev); # transpose
      colnames(exp.mat3) = colnames(exp.mat) # set the column name
      exp.mat4 = exp.mat3[-which(is.na(exp.mat3)),] # remove NaN rows
      exp.mat5 = exp.mat4[rownames(exp.mat4) %in% gene_symbol,] # grab gene symbols
      # reassign data
      dataset = exp.mat5
      ## generate color for pheatmap
      if (abs(min(dataset)) > abs(max(dataset))) {
        dataset[dataset < -abs(max(dataset))] = -abs(max(dataset))
      } else {
        dataset[dataset > abs(min(dataset))] = abs(min(dataset))
      }
      zscore_range = 10;
      minimum = -zscore_range;
      maximum = zscore_range;
      bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
      hmcols<- colorRampPalette(c("dark blue","blue","white","red", "dark red"))(length(bk)-1)
      pheatmap(dataset,
               cluster_col = T,
               cluster_row = T,
               fontsize_row = 9,
               fontsize_col = 12,
               show_rownames = T,
               show_colnames = T,
               color=hmcols,
               angle_col = 45,
               cellwidth = 75,
               cluster_cols = F)
     })
    
     output$ssGSEAtab <- DT::renderDataTable({
       exp.mat <- as.matrix(expr_data)
       gs <- list()
       for (i in unique(gmt$term)){
         gs[[i]] <- gmt[gmt$term == i,]$gene
       }
       var <- input$ssmethod
       if (var == "gsva"){
         ssgsea <- gsva(exp.mat, gs, method = "gsva", verbose = F)
       }
       if (var == "ssgsea"){
         ssgsea <- gsva(exp.mat, gs, method = "ssgsea", verbose = F)
       }
       if (var == "zscore"){
         ssgsea <- gsva(exp.mat, gs, method = "zscore", verbose = F)
       }
       if (var == "plage"){
         ssgsea <- gsva(exp.mat, gs, method = "plage", verbose = F)
       }
       ssgsea <- gsva(exp.mat, gs, method = var, verbose = F)
       ssgsea_tib <- as.data.frame(ssgsea)
       DT::datatable(ssgsea,
                     caption = "Single Sample GSEA Enrichment Scores")
     })
     
     output$heatmap1 <- renderPlot({
       exp.mat <- as.matrix(expr_data)
       gs <- list()
       for (i in unique(gmt$term)){
         gs[[i]] <- gmt[gmt$term == i,]$gene
       }
       var <- input$ssmethod
       if (var == "gsva"){
         gsea.u <- gsva(exp.mat, gs, method = "gsva", verbose = F)
       }
       if (var == "ssgsea"){
         gsea.u <- gsva(exp.mat, gs, method = "ssgsea", verbose = F)
       }
       if (var == "zscore"){
         gsea.u <- gsva(exp.mat, gs, method = "zscore", verbose = F)
       }
       if (var == "plage"){
         gsea.u <- gsva(exp.mat, gs, method = "plage", verbose = F)
       }
       gsea.u <- gsva(exp.mat, gs, method = var, verbose = F)
       pheatmap(gsea.u,
                fontsize = 15)
     })
     
     output$boxplot0 <- renderPlot({
       ####----Reformatting GMT File to Gene Set List----####
       gmt <- gmt
       gs <- list()
       for (i in unique(gmt$term)){
         gs[[i]] <- gmt[gmt$term == i,]$gene
       }
       ####----Reformatting Expression Data----####
       expr_data <- expr_data
       exp.mat = as.matrix(expr_data)
       ####----Perform ssGSEA----####
       ssgsea <- gsva(exp.mat, gs, method = "ssgsea", verbose = F)
       ####----Boxplot----####
       # prep data
       ssgsea2 <- as.data.frame(t(ssgsea))
       ssgsea3 <- ssgsea2 %>% 
         mutate(type = case_when(
           grepl("CG2.LUC", rownames(ssgsea2)) ~ "CG2.LUC",
           grepl("CG2.RAS", rownames(ssgsea2)) ~ "CG2.RAS",
           grepl("MIC.LUC", rownames(ssgsea2)) ~ "MIC.LUC"
         ))
       ssgsea3 <- ssgsea3 %>%
         relocate(type)
       var <- input$sstype
       if (var == "Gruber_Lab_HSC"){
         ggplot(ssgsea3, aes(factor(type), Gruber_Lab_HSC, fill = type)) +
           geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
           geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 1) +
           ggtitle("Single Sample GSEA Expression Signature: Gruber_Lab_HSC") +
           theme(plot.title=element_text(size=10),
                 axis.text.x = element_text(size=10, angle=90),
                 axis.text.y = element_text(size=10),
                 axis.title = element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10)) +
           theme_bw() +
           stat_compare_means()
       }
       if (var == "Gruber_Lab_CMP"){
         ggplot(ssgsea3, aes(factor(type), Gruber_Lab_CMP, fill = type)) +
           geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
           geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 1) +
           ggtitle("Single Sample GSEA Expression Signature: Gruber_Lab_CMP") +
           theme(plot.title=element_text(size=10),
                 axis.text.x = element_text(size=10, angle=90),
                 axis.text.y = element_text(size=10),
                 axis.title = element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10)) +
           theme_bw() +
           stat_compare_means()
       }
       if (var == "Gruber_Lab_GMP"){
         ggplot(ssgsea3, aes(factor(type), Gruber_Lab_GMP, fill = type)) +
           geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
           geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 1) +
           ggtitle("Single Sample GSEA Expression Signature: Gruber_Lab_GMP") +
           theme(plot.title=element_text(size=10),
                 axis.text.x = element_text(size=10, angle=90),
                 axis.text.y = element_text(size=10),
                 axis.title = element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10)) +
           theme_bw() +
           stat_compare_means()
       }
       if (var == "Gruber_Lab_MPP"){
         ggplot(ssgsea3, aes(factor(type), Gruber_Lab_MPP, fill = type)) +
           geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
           geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 1) +
           ggtitle("Single Sample GSEA Expression Signature: Gruber_Lab_MPP") +
           theme(plot.title=element_text(size=10),
                 axis.text.x = element_text(size=10, angle=90),
                 axis.text.y = element_text(size=10),
                 axis.title = element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10)) +
           theme_bw() +
           stat_compare_means()
       }
       if (var == "Gruber_Lab_MEP"){
         ggplot(ssgsea3, aes(factor(type), Gruber_Lab_MEP, fill = type)) +
           geom_boxplot(width = 0.5, lwd = 1, fill = "white") +
           geom_dotplot(binaxis = 'y', stackdir = "center", dotsize = 1) +
           ggtitle("Single Sample GSEA Expression Signature: Gruber_Lab_MEP") +
           theme(plot.title=element_text(size=10),
                 axis.text.x = element_text(size=10, angle=90),
                 axis.text.y = element_text(size=10),
                 axis.title = element_text(size=10),
                 legend.text = element_text(size=10),
                 legend.title = element_text(size=10)) +
           theme_bw() +
           stat_compare_means()
       }
     })
}


shinyApp(ui = ui, server = server)


