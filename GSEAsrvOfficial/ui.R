
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
    
navbarPage("{Neault-Gruber Lab Collab}",
  tabPanel("Intro/Methodology",
           fluidPage(
               mainPanel(
                   h2("RNAseq Analysis Method"),
		   p("RNA-seq reads from the mouse models were mapped by (Please insert your own methodology). Mapped reads were counted with (insert quantification strategy here) and normalized to TPM (transcript per kilobase million mapped reads). Differentially up-regulated genes in CG2.LUC and CG2.RAS were computed by LIMMA in R 3.0.1, comparing each murine model to normal BM in murine models. Pathway enrichment analysis was performed by GSEA [PMID: 16199517] and enrichR[PMID: 23586463]. Single-sample gene set enrichment analysis (ssGSEA) [PMID: 19847166, PMID: 30595505] was used to quantify the expression signatures with the up-regulated gene set derived from RNAseq samples of flow-sorted mouse hematopoietic progenitors, including hematopoietic stem cell (HSC), Granulocyte-monocyte progenitors (GMP), multipotent common myeloid progenitor (CMP), megakaryocyte–erythroid progenitor cell (MEP), and multipotent progenitors (MPP) [request Tanja Gruber’s lab to fill out additional details]. LIMMA was used to define differentially expressed genes [see supplemental table for the gene sets]. ")
               ))),
  tabPanel("GSEA Analysis",
           fluidPage(
               title = "GSEA Analysis using Tanja Gruber's Murine Gene Set",
               sidebarPanel(
                   width = 3,
		   selectInput("Comparisons", "Comparisons",
                c("CG2.LUC+CG2.RAS vs BM_CNTRL",
                  "CG2.LUC vs BM_CNTRL",
                  "CG2.RAS vs BM_CNTRL")),
                   selectInput("geneSet", "Gene Set", choices = c("Gruber_Lab_MEP", "Gruber_Lab_HSC", "Gruber_Lab_CMP", "Gruber_Lab_MPP", "Gruber_Lab_GMP"))),
               mainPanel(
                   tabsetPanel(
                       id = "dataset",
                       tabPanel("Enrichment Plot", h3("GSEA Enrichment Plot"), withSpinner(verbatimTextOutput("NESandPval"), type = 6),
                                withSpinner(plotOutput("enrichplot0", width = "500px", height = "450px"), type = 6),
				h3("Leading Edge Genes (~Signal2Noise Ranking)"),
                                div(DT::dataTableOutput("LeadingEdgeGenes"), style = "font-size:10px; height:500px; overflow-y: scroll"),
                                downloadButton("downloadGeneList", "Download Leading Edge Gene List")),
                       tabPanel("Gene Expression Heatmap", withSpinner(plotOutput("heatmap0", width = "100%", height = "2000px"), type = 6)),
                       tabPanel("GSEA Summary Table", div(DT::dataTableOutput("enrich_sig_table"), style = "font-size:12px"))),
                       ))),
  tabPanel("Single Sample GSEA Analysis",
           fluidPage(
             title = "Single Sample GSEA Analysis",
             sidebarPanel(
               #selectInput("ssmethod", "Single Sample GSEA Enrichment Scoring Method", choices = c("ssgsea","gsva","zscore","plage")),
               #selectInput("ssmethod", "Single Sample GSEA Enrichment Scoring Method", choices = c("ssgsea","gsva")),
               selectInput("ssmethod", "Single Sample GSEA Enrichment Scoring Method", choices = c("ssgsea")),
               selectInput("ssnorm", "Single Sample GSEA Normalization", choices = c("yes")),
               selectInput("boxplot_pval_flag", "Show Boxplot Pval - Wilcox.Test", choices = c("yes", "no"))
               #selectInput("ssnorm", "Single Sample GSEA Normalization", choices = c("yes", "no"))
#               fileInput("gmt.u", "Choose GMT File", accept = ".gmt"),
#               fileInput("expr.u", "Choose Expression Data File - tab delimited", accept = c(".tsv", ".txt"))
               ),
             mainPanel(
               tabsetPanel(
                 id = "dataset",
                    tabPanel("ssGSEA Boxplot", 
                             selectInput("sstype", "Gene Set", choices = c("Gruber_Lab_MEP", "Gruber_Lab_HSC", "Gruber_Lab_CMP", "Gruber_Lab_MPP", "Gruber_Lab_GMP")),
			     h3("Boxplot of the ssGSEA Scores"),
			     withSpinner(plotOutput("boxplot1", width = "500px", height = "500px"), type = 6),
			     withSpinner(plotOutput("boxplot0", width = "500px", height = "500px"), type = 6)),
                    tabPanel("ssGSEA Heatmap", 
			     h3("Heatmap of the ssGSEA Scores (Z-score normalized across samples)"),
			     withSpinner(plotOutput("heatmap1", width = "100%", height = "550px"), type = 6)),
                    tabPanel("Single Sample GSEA Raw Output", 
			     h3("ssGSEA Score Table"),
			     div(DT::dataTableOutput("ssGSEAtab"), style = "font-size:10px")))
             ))))




