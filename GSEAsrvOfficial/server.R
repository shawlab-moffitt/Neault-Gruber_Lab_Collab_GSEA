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

# Gene Set Data
gmt <- read.gmt("Up_Reg_Myeloid_Cells_Gruber_Lab.gmt")
gmt$ont <- gsub("MEP", "Gruber_Lab_MEP", gmt$ont)
# Expression Data
expr_data <- read.table("BM_Mouse_Model_TPM_20210816_reorder.max.txt", sep="\t", header=T, row.names=1, quote="")



server <- function(input, output, session) {
    output$enrich_sig_table <- DT::renderDataTable({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- list();
      groupA_name = "CG2.LUC+CG2.RAS"; 
      if (input$Comparisons == "CG2.LUC vs BM_CNTRL") {
          groupA_name = "CG2.LUC"; 
          groupA <- colnames(A)[grep("CG2.LUC", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.RAS vs BM_CNTRL") {
          groupA_name = "CG2.RAS"; 
          groupA <- colnames(A)[grep("CG2.RAS", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.LUC+CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups

      }
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
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = 1.0)
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
      groupA <- list();
      groupA_name = "CG2.LUC+CG2.RAS";
      if (input$Comparisons == "CG2.LUC vs BM_CNTRL") {
          groupA_name = "CG2.LUC";
          groupA <- colnames(A)[grep("CG2.LUC", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.RAS vs BM_CNTRL") {
          groupA_name = "CG2.RAS";
          groupA <- colnames(A)[grep("CG2.RAS", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.LUC+CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups

      }
     
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
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = 1.0)
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
      NES.o <- paste0("NES: ", NES)
      Pval.o <- paste0("Pvalue: ", Pval)
      if (NES > 0){
        UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is upregulated in [", groupA_name, "] group.")
      }
      else {
        UpOrDown <- paste("Based on the normalized enrichment score above, the", GS, "gene set is downregulated in [", groupA_name, "] group.")
      }
      paste(NES.o, Pval.o, UpOrDown, sep = '\n')
    })
    
    output$enrichplot0 <- renderPlot({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- list(); 
      if (input$Comparisons == "CG2.LUC vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2.LUC", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2.RAS", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.LUC+CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups

      }
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
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = 1.0)
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
      groupA <- list(); 
      if (input$Comparisons == "CG2.LUC vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2.LUC", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2.RAS", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.LUC+CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups

      }
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
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = 1.0)
      gsea.df <- as.data.frame(gsea.res@result)
      GS <- input$geneSet
      ## Subset core enriched genes
      genes1 <- as.matrix(gsea.df[which(gsea.df$Description==GS),"core_enrichment"])
      genes2 <- strsplit(genes1,"/")
      GeneSymbol <- as.data.frame(genes2, col.names = "GeneSymbol")
      GeneSymbol$Rank <- rownames(GeneSymbol)
      GeneSymbol <- GeneSymbol[,c("Rank","GeneSymbol")]
      DT::datatable(GeneSymbol, options = list(paging = F), rownames = F)
    })
    
    output$downloadGeneList <- downloadHandler({
      ####----Data----####
      exp.mat <- as.matrix(expr_data)
      A = exp.mat
      groupA <- list(); 
      if (input$Comparisons == "CG2.LUC vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2.LUC", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2.RAS", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.LUC+CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups

      }      
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
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = 1.0)
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
      groupA <- list(); 
      if (input$Comparisons == "CG2.LUC vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2.LUC", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2.RAS", colnames(A))] #assign groups
      } else if (input$Comparisons == "CG2.LUC+CG2.RAS vs BM_CNTRL") {
          groupA <- colnames(A)[grep("CG2", colnames(A))] #assign groups

      }      
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
      gsea.res <- GSEA(s2n.matrix.s, TERM2GENE = gmt, verbose = F, pvalueCutoff = 1.0)
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
      dataset = exp.mat5[,c(groupA,groupB)]
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
               cluster_col = F,
               cluster_row = F,
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
       for (i in unique(gmt$ont)){
         gs[[i]] <- gmt[gmt$ont == i,]$gene
       }
       var <- input$ssmethod
       normflag = TRUE;
       if (input$ssnorm == "no") {
	       normflag = FALSE;
       }
       #if (var == "gsva"){
       #  ssgsea <- gsva(exp.mat, gs, method = "gsva", ssgsea.norm = normflag, verbose = F)
       #}
       #if (var == "ssgsea"){
       #  ssgsea <- gsva(exp.mat, gs, method = "ssgsea", ssgsea.norm = FALSE, verbose = F)
       #}
       #if (var == "zscore"){
       #  ssgsea <- gsva(exp.mat, gs, method = "zscore", ssgsea.norm = normflag, verbose = F)
       #}
       #if (var == "plage"){
       #  ssgsea <- gsva(exp.mat, gs, method = "plage", ssgsea.norm = normflag, verbose = F)
       #}
       ssgsea <- gsva(exp.mat, gs, method = var, ssgsea.norm = normflag, verbose = F)
       ssgsea_tib <- as.data.frame(ssgsea)
       DT::datatable(ssgsea,
                     caption = "Single Sample GSEA Enrichment Scores")
     })
     
     output$heatmap1 <- renderPlot({
       exp.mat <- as.matrix(expr_data)
       gs <- list()
       for (i in unique(gmt$ont)){
         gs[[i]] <- gmt[gmt$ont == i,]$gene
       }
       var <- input$ssmethod
       normflag = TRUE;
       if (input$ssnorm == "no") {
               normflag = FALSE;
       }
       
       #if (var == "gsva"){
       #  gsea.u <- gsva(exp.mat, gs, method = "gsva", verbose = F)
       #}
       #if (var == "ssgsea"){
       #  gsea.u <- gsva(exp.mat, gs, method = "ssgsea", verbose = F)
       #}
       #if (var == "zscore"){
       #  gsea.u <- gsva(exp.mat, gs, method = "zscore", verbose = F)
       #}
       #if (var == "plage"){
       #  gsea.u <- gsva(exp.mat, gs, method = "plage", verbose = F)
       #}
       
       gsea.u <- gsva(exp.mat, gs, method = var, ssgsea.norm = normflag, verbose = F)
       allDat = t(gsea.u)
       scaled = apply(allDat, 2, scale);
       all = apply(scaled, 1, rev)
       colnames(all) = rownames(allDat)
       library(pheatmap)
       minimum = min(all);
       maximum = max(all);
       bk = c(seq(minimum,minimum/2, length=100), seq(minimum/2,maximum/2,length=100),seq(maximum/2,maximum,length=100))
       len = 50
       myBreaks <- c(seq(min(all), 0, length.out=ceiling(len/2) + 1),seq(max(all)/len, max(all), length.out=floor(len/2)))
       hmcols<- colorRampPalette(c("green","black","red"))(length(myBreaks)-1)
       pheatmap(all, cluster_col = F, cluster_row = T, fontsize_row = 12, fonrsize_col = 12, color=hmcols,breaks=myBreaks)
       #pheatmap(gsea.u,
       #         fontsize = 15)
     })
     
     output$boxplot0 <- renderPlot({
       ####----Reformatting GMT File to Gene Set List----####
       gmt <- gmt
       gs <- list()
       for (i in unique(gmt$ont)){
         gs[[i]] <- gmt[gmt$ont == i,]$gene
       }
       ####----Reformatting Expression Data----####
       expr_data <- expr_data
       exp.mat = as.matrix(expr_data)
       ####----Perform ssGSEA----####

       var <- input$ssmethod
       normflag = TRUE;
       if (input$ssnorm == "no") {
               normflag = FALSE;
       }

       pval_flag = TRUE;
       if (input$boxplot_pval_flag == "no") {
               pval_flag = FALSE;
       }
       ssgsea <- gsva(exp.mat, gs, method = var, ssgsea.norm = normflag, verbose = F)
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
       var1 <- input$sstype
       if (var1 == "Gruber_Lab_HSC"){
	 if (pval_flag) {
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
           stat_compare_means(comparisons=list(c("CG2.RAS", "MIC.LUC"), c("CG2.LUC", "MIC.LUC")), method="wilcox.test")

	 } else {
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
           theme_bw()

	 }
       } else if (var1 == "Gruber_Lab_CMP"){
        if (pval_flag) {
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
           stat_compare_means(comparisons=list(c("CG2.RAS", "MIC.LUC"), c("CG2.LUC", "MIC.LUC")), method="wilcox.test")
	} else {
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
           theme_bw()

	}
       } else if (var1 == "Gruber_Lab_GMP"){
	if (pval_flag) {
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
           #stat_compare_means()
           stat_compare_means(comparisons=list(c("CG2.RAS", "MIC.LUC"), c("CG2.LUC", "MIC.LUC")), method="wilcox.test")
	} else {
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
           theme_bw()
	}
       } else if (var1 == "Gruber_Lab_MPP"){
        if (pval_flag) {
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
           #stat_compare_means()
           stat_compare_means(comparisons=list(c("CG2.RAS", "MIC.LUC"), c("CG2.LUC", "MIC.LUC")), method="wilcox.test")
	} else {
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
           theme_bw()
	}
       } else if (var1 == "Gruber_Lab_MEP"){
	if (pval_flag) {
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
           #stat_compare_means()
           stat_compare_means(comparisons=list(c("CG2.RAS", "MIC.LUC"), c("CG2.LUC", "MIC.LUC")), method="wilcox.test")

	} else {
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
           theme_bw()
	}
       }


     })

     output$boxplot1 <- renderPlot({
       ####----Reformatting GMT File to Gene Set List----####
       gmt <- gmt
       gs <- list()
       for (i in unique(gmt$ont)){
         gs[[i]] <- gmt[gmt$ont == i,]$gene
       }
       ####----Reformatting Expression Data----####
       expr_data <- expr_data
       exp.mat = as.matrix(expr_data)
       ####----Perform ssGSEA----####
       var <- input$ssmethod
       normflag = TRUE;
       if (input$ssnorm == "no") {
               normflag = FALSE;
       }
       pval_flag = TRUE;
       if (input$boxplot_pval_flag == "no") {
               pval_flag = FALSE;
       }

       ssgsea <- gsva(exp.mat, gs, method = var, ssgsea.norm = normflag, verbose = F)
       ####----Boxplot----####
       # prep data
       ssgsea2 <- as.data.frame(t(ssgsea))
       ssgsea3 <- ssgsea2 %>%
         mutate(type = case_when(
           grepl("CG2.LUC", rownames(ssgsea2)) ~ "CG2.LUC+CG2.RAS",
           grepl("CG2.RAS", rownames(ssgsea2)) ~ "CG2.LUC+CG2.RAS",
           grepl("MIC.LUC", rownames(ssgsea2)) ~ "MIC.LUC"
         ))
       ssgsea3 <- ssgsea3 %>%
         relocate(type)
       var1 <- input$sstype
       if (var1 == "Gruber_Lab_HSC"){
       min = min(ssgsea3$Gruber_Lab_HSC)
       max = max(ssgsea3$Gruber_Lab_HSC)
       diff = max - min;
       min = min - (diff * 0.05)
       max = max + (diff * 0.1)	
        if (pval_flag) {       
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
           stat_compare_means(method="wilcox.test") + ylim(c(min, max))
	} else {
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
           ylim(c(min, max))

	}
       } else if (var1 == "Gruber_Lab_CMP"){
	 min = min(ssgsea3$Gruber_Lab_CMP)
         max = max(ssgsea3$Gruber_Lab_CMP)
         diff = max - min;
         min = min - (diff * 0.05)
         max = max + (diff * 0.1)
	 if (pval_flag) {
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
           stat_compare_means(method="wilcox.test") + ylim(c(min,max))
	 } else {
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
           ylim(c(min,max))
	 }
       } else if (var1 == "Gruber_Lab_GMP"){

	 min = min(ssgsea3$Gruber_Lab_GMP)
         max = max(ssgsea3$Gruber_Lab_GMP)
         diff = max - min;
         min = min - (diff * 0.05)
         max = max + (diff * 0.1)
         if (pval_flag) {
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
           #stat_compare_means()
           stat_compare_means(method="wilcox.test") + ylim(c(min, max))
	 } else {
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
           #stat_compare_means()
           ylim(c(min, max))
	 }
       } else if (var1 == "Gruber_Lab_MPP"){
         min = min(ssgsea3$Gruber_Lab_MPP)
         max = max(ssgsea3$Gruber_Lab_MPP)
         diff = max - min;
         min = min - (diff * 0.05)
         max = max + (diff * 0.1) 
	 if (pval_flag) {
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
           #stat_compare_means()
           stat_compare_means(method="wilcox.test") + ylim(c(min, max))
	 } else {
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
           #stat_compare_means()
           ylim(c(min, max))
	 }
       } else if (var1 == "Gruber_Lab_MEP"){
         min = min(ssgsea3$Gruber_Lab_MEP)
         max = max(ssgsea3$Gruber_Lab_MEP)
         diff = max - min;
         min = min - (diff * 0.05)
         max = max + (diff * 0.1) 
	 if (pval_flag) {
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
           #stat_compare_means()
           stat_compare_means(method="wilcox.test") + ylim(c(min, max))
	 } else {
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
           #stat_compare_means()
           ylim(c(min, max))
	 }
       }


     })





}


