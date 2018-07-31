####load packages####
library(shiny)
library(ggplot2)
library(reshape2)
library(grid)
#

#### read data for expression plot####
rpkm_raw <- read.csv("./data/C9_rpkm_ave.csv", stringsAsFactors = FALSE)
rpkm <- melt(rpkm_raw)
rpkm$geno_LS <- c(rep("C9285_4LS", nrow(rpkm) / 3),
                 rep("C9285_6LS", nrow(rpkm) / 3),
                 rep("T65_6LS", nrow(rpkm) / 3))

rpkm$timepoints <- rep(c(rep("00h", nrow(rpkm_raw)),
                        rep("01h", nrow(rpkm_raw)),
                        rep("03h", nrow(rpkm_raw)),
                        rep("06h", nrow(rpkm_raw)),
                        rep("12h", nrow(rpkm_raw)),
                        rep("24h", nrow(rpkm_raw))),3)
rpkm_SD <- read.csv("./data/C9_rpkm_sd.csv", stringsAsFactors = FALSE)
rpkm_SD <- melt(rpkm_SD)
rpkm$SD <- rpkm_SD$value
names(rpkm)[1] <- "ID"

#####

#### read data for chromPlot####

# actual C9 chrom length
chrom_length <- data.frame("chrom" = 1:12,
                           "length" = c(42176892, 35787457, 36858062,
                                        32253451, 28753270, 29951177,
                                        28043875, 27836186, 22786783,
                                        23375750, 28703814, 25723258))

# QTL regions based on NB
QTL_regions <- data.frame("QTL" = c("none","qGTIL3","qGLEI3","qGNEI3","qGLEI8",
                                    "qGTIL9","qGLEI9","qGNEI9","qGNEI10","qGTIL12",
                                    "qTIL1","qRIE1","qLEI3","qTIL12","qNEI12",
                                    "qLEI12","qTIL2","qTIL4"),
                          "chrom"=c(NA,3,3,3,8,9,9,9,10,12,1,1,3,12,12,12,2,4),
                          "QTL_left" = c(NA,6496675,11316497,11316497,9632537,7048026,11759322,
                                         7048026,3356964,10614890,29065902,39466088,10708808,
                                         25203610,25203610,25203610,33636231,23596814),
                          "QTL_right" = c(NA,17514602,17092718,17092718,27994362,19252470,15349185,
                                          19252470,13600514,26133898,44941618,42591824,12954471,
                                          26368115,26368115,26368115,35997665,25863430))

# also based on NB
fixed_GOI <- data.frame("gene" = c("SK1/2","SK3", "SD1"),
                        "chrom" = c(12,3,1),
                        "pos" = c(25380994, 12930367, 38383927))

chrom_pos <- read.csv("data/C9_annotations.csv", stringsAsFactors = FALSE)
names(chrom_pos)[1] <- "ID"
#####

#### theme for plotting####
theme_custom <- theme_grey() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 14, colour = "black", face = "bold"),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks =element_line(colour = "black"),
        panel.grid.major = element_line(colour="grey95"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5,linetype="solid"),
        strip.text = element_text(size = 14, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.position = "bottom"
  )
#####

##### prepare MM table#####
MM_table_option <- list(
  columnDefs = list(list(targets = c(2,3,4,5,6)-1,
                         searchable = FALSE)),
  pageLength = 10)
load("./data/bg_MM_C9V1")
bg_MM <- MM_C9_peptides
#####


shinyServer(function(input, output) {
  
  # return chr with input genes
  inGenes <- reactive({
    
    genesToPlot <- unlist(strsplit(input$genesToPlotDirect,"\n"))
    toupper(unique(genesToPlot))
    
  })
  
  # expression plot
  output$expressionPlot <- renderPlot({
    
    if (input$genesToPlotDirect == "")
      return(NULL)
    if (is.null(input$numFacetCol))
      return(NULL)
    
    rpkmToPlot <- rpkm[rpkm$ID %in% inGenes(),]
    if (nrow(rpkmToPlot) == 0)
      return(NULL)
    if (length(inGenes()) >= 50)
      return(NULL)
    
    rpkmToPlot <- rpkmToPlot[rpkmToPlot$geno_LS %in% input$whichGenoLS,]
    
    if(input$plotSort == "inOrder"){
      rpkmToPlot$ID <- factor(rpkmToPlot$ID, inGenes())
    }
    if (input$facetswitch == "transcript") { 
      ggplot(data=rpkmToPlot, aes(x=timepoints,
                                 y=value,
                                 group = geno_LS,
                                 color = geno_LS)) +
        geom_line(size = 1, stat = "identity")+
        facet_wrap(~ID, scales = "free", ncol = as.integer(input$numFacetCol))+
        geom_errorbar(aes(ymin = value - SD, ymax = value + SD),
                      width = 0.2, colour = "grey25")+
        scale_x_discrete(labels = c("0","1","3","6","12","24"))+
        xlab("time after submergence (h)")+
        ylab("rpkm")+
        expand_limits(y = 0)+
        theme_custom
    } else {
      ggplot(data=rpkmToPlot, aes(x=timepoints,
                                 y=value,
                                 group = ID,
                                 color = ID)) +
        geom_line(size = 0.4, stat = "identity")+
        facet_wrap(~geno_LS, scales = "free", ncol = as.integer(input$numFacetCol))+
        geom_errorbar(aes(ymin = value - SD, ymax = value + SD),
                      width = 0.2, colour = "black")+
        xlab("time after submergence (h)")+
        ylab("rpkm")+
        expand_limits(y = 0)+
        theme_custom
    }
  }, height = reactive({input$exprPlotHeight})
  )
  
  #### table expression plot ####
  output$expr_transcript_info <- renderDataTable(
    
    chrom_pos[chrom_pos$ID %in% inGenes(),],
    options = list(paging = FALSE)
  )
  ##### 
  
  # chromosome/QTL plot  
  output$chromPlot <- renderPlot({
    
    if (input$genesToPlotDirect == "")
      return(NULL)
    if (is.null(input$numFacetCol))
      return(NULL)
    
    chrom_pos_to_plot <- chrom_pos[chrom_pos$ID %in% inGenes(),]
    
    QTL_regions_to_plot <- QTL_regions[QTL_regions$QTL %in% input$whichQTLregion,]
    
    ggplot()+
      geom_segment(data = chrom_length, aes(x = chrom, y = 0, # chrom backbone
                                              xend = chrom, yend = length))+
        geom_rect(data = QTL_regions_to_plot, aes(xmin = chrom - 0.2, xmax = chrom + 0.2, #QTL regions
                                                  ymin = QTL_left, ymax = QTL_right),
                  alpha = 0.2)+
        geom_segment(data = fixed_GOI, aes(x = chrom - 0.1, xend = chrom + 0.1, #fixed GOIs
                                           y = pos, yend = pos), color = "orange", size = 1)+
        geom_text(data = fixed_GOI, aes(x = chrom - 0.3, y = pos, label = gene,
                                        angle = 90))+ #names for fixed GOIs
        
        geom_segment(data = chrom_pos_to_plot, aes(x = chrom - 0.1, xend = chrom + 0.1,
                                                   y = end-(end-start), yend = end-(end-start)))+
        scale_x_discrete(limits = c(1:12))+
        scale_y_reverse(labels = function(x){x/1000000})+
        xlab("chromosome")+
        ylab("Mb")+
        theme_custom
  })
  
  #### transcript info table ####
  output$transcript_info <- renderDataTable(
    
    brushedPoints(chrom_pos[chrom_pos$ID %in% inGenes(),], input$transcript_brush,
                  xvar = "chrom", yvar = "start"),
    options = list(searching = FALSE, paging = FALSE)
  )
  #####
  
  #### expression table ####
  output$expression_table <- renderDataTable({
      if (input$genesToPlotDirect == "")
      return(NULL)
      rpkmToPlot <- rpkm_raw[rpkm_raw[,1] %in% inGenes(),]
      rpkmToPlot
    })
  ######
  
  ##### MM tables #####
  
  output$MM_lev1_out <- renderDataTable({
    
    set_MM <- bg_MM[bg_MM$C9_ID %in% inGenes(),]
    
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name1),
                          "BinBg" = as.vector(table(bg_MM$bin_name1)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name1)),
                          "NotBinBg" = nrow(bg_MM) - as.vector(table(bg_MM$bin_name1)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    if (nrow(testDF) == 0)
      return(NULL)
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    colnames(testDF) <- c("MAPMAN level 1 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichment", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
  }, options = MM_table_option)
  
  #####################
  
  output$MM_lev2_out <- renderDataTable({
    
    set_MM <- bg_MM[bg_MM$C9_ID %in% inGenes(),]
    
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name2),
                          "BinBg" = as.vector(table(bg_MM$bin_name2)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name2)),
                          "NotBinBg" = nrow(bg_MM) - as.vector(table(bg_MM$bin_name2)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    if (nrow(testDF) == 0)
      return(NULL)
    
    FET_P <- NULL
    
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    
    colnames(testDF) <- c("MAPMAN level 2 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichment", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
  }, options = MM_table_option)
  
  #####################
  
  output$MM_lev3_out <- renderDataTable({
    
    set_MM <- bg_MM[bg_MM$C9_ID %in% inGenes(),]
    
    testDF <- data.frame( "BinSet" = table(set_MM$bin_name3),
                          "BinBg" = as.vector(table(bg_MM$bin_name3)),
                          "NotBinSet" =  nrow(set_MM) - as.vector(table(set_MM$bin_name3)),
                          "NotBinBg" = nrow(bg_MM) - as.vector(table(bg_MM$bin_name3)))
    testDF$OR_FC <- round((testDF[,2] / nrow(set_MM)) / (testDF[,3] / nrow(bg_MM)), 2)
    testDF <- testDF[testDF$BinSet.Freq > 0,] #remove bins with 0 genes
    testDF <- testDF[!grepl(pattern = "undefined", testDF$BinSet.Var1),] #remove undefined sub bins
    
    if (nrow(testDF) == 0)
      return(NULL)
    
    FET_P <- NULL
    for(i in 1:nrow(testDF)) {
      
      BinSet <- testDF[i,2]
      BinBg <- testDF[i,3]
      NotBinSet <- testDF[i,4]
      NotBinBg <- testDF[i,5]
      
      FET_P[i] <- fisher.test(matrix(c(BinSet,BinBg,NotBinSet,NotBinBg),nrow = 2),
                              alternative = "greater")$p.value
    }
    testDF$FET_P <- round(FET_P, 3)
    testDF$FDR_P <- round(p.adjust(FET_P, method = "BH"), 3)
    
    colnames(testDF) <- c("MAPMAN level 3 bin",
                          "Genes in set", "Genes in total",
                          "Other genes in set", "Other genes in total",
                          "Fold enrichment", "P-value (FET)", "corrected P-value (BH)")
    testDF[order(testDF[,8]),]
  }, options = MM_table_option)
  
  #####################
  
})

