
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




