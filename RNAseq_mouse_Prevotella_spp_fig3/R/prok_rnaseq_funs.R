
## MA_plot function ### keep an eye fc is equal to fc not to 2fold change!!!
MA_plot <- function(x=diff_express, color=color_palette, main_lab=plot_lab, gename="gename"){
  p = ggmaplot(x, main = plot_lab,
               fdr = 0.01, fc = 2, size = 1.5,
               palette = color,
               genenames = as.vector(x[[gename]]),
               legend = "top", top = 25, label.rectangle = FALSE,
               font.label = c("bold", 10),
               font.legend = "bold",
               font.main = "bold",
               ggtheme = ggplot2::theme_minimal(30))
  p1=p+geom_point(alpha = 0.1) +
    theme(axis.text.x = element_text(size = 18), #angle = 1,hjust = 1,
          plot.title = element_text(size = 18, hjust = 0.5),     
          axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),      
          legend.text=element_text(size = 18),
          legend.title=element_blank(), 
          strip.text.x =element_text(size = 18),
          axis.text.y =element_text(size = 18),
          axis.line.x = element_line(colour = "black", size=.5),
          axis.line.y = element_line(colour = "black", size=.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) # + ylim(-6,6)
  return(p1)
}


## returns a subseted and filtered dataset dss for results and contrast
select_dds_pair <- function(x=ddsFull, control = control_con, expcond = exp_con1, 
                            minreads = 10, minsamples = 4) {
ddsexp = x[ , x$condition == control | x$condition == expcond ]
colData(ddsexp) <- colData(ddsexp)[(1:3)]
as.data.frame(colData(ddsexp))
ddsexp$condition <- droplevels(ddsexp$condition)
#ddsexp$sizeFactor <- droplevels(ddsexp$sizeFactor)
#ddsFCseq2$Treatment <- droplevels(ddsFCseq2$Treatment)
design(ddsexp) <- ~ condition
ddsexp$condition<- relevel(ddsexp$condition, ref = control)

dds <- DESeq(ddsexp)

dds_ko10 <- dds[ rowSums(counts(dds, normalized = TRUE) >= minreads) >= minsamples, ]  

return(dds_ko10)
}

###########################

hi_loadings_plot <- function(pcaobj=pcaobj, topN = 10, whichpc = 1) {
  df <- as.data.frame(get_top_pca_genes(pcaobj, topN = topN, whichpc = whichpc))
  df1 <- data.table::setDT(df, keep.rownames = TRUE)[] 
  colnames(df1) <-  c("gene_id", "pca_loadings")
  
  df1$gene_id <- factor(df1$gene_id, levels = c(df[[1]])) 
  
  p <- ggplot(df1, aes(x = df1$gene_id, y = df1$pca_loadings, 
                       fill = ifelse(df1$pca_loadings < 0,"#B31B21", "steelblue"))) +
    geom_col(col= "black") +
    scale_fill_manual(values = c("steelblue", "#B31B21")) +
    theme_minimal(14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          #axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          #axis.title.y = element_blank(),
          legend.position = "",
          legend.title = element_blank(),
          axis.line.y = element_line(color="black", size = 0.2),
          #panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
          panel.grid.minor.x = element_blank()) +
    ylab("PC1 loadings") +
    xlab("GeneID")
  return(p)
}

##########################

## modified function to extract gene features

get_top_pca_genes <- function (pcaobj, whichpc = 1, topN = 10, exprTable = NULL, annotation = NULL, 
                               title = "Top/bottom loadings") 
{
  if (whichpc < 0) 
    stop("Use a positive integer value for the principal component to select")
  if (whichpc > nrow(pcaobj$x)) 
    stop("You can not explore a principal component that is not in the data")
  geneloadings_sorted <- sort(pcaobj$rotation[, whichpc])
  geneloadings_extreme <- c(sort(head(geneloadings_sorted, topN), decreasing = T), 
                            sort(tail(geneloadings_sorted, topN), decreasing = T))
  if (!is.null(exprTable)) {
    tab <- exprTable[names(geneloadings_extreme), ]
    if (!is.null(annotation)) 
      rownames(tab) <- annotation$gene_name[match(rownames(tab), 
                                                  rownames(annotation))]
    return(tab)
  }
  if (!is.null(annotation)) 
    names(geneloadings_extreme) <- annotation$gene_name[match(names(geneloadings_extreme), 
                                                              rownames(annotation))]
  return(geneloadings_extreme)
}


######################

## getting gff features 
require(rtracklayer)
require(GenomicRanges)
get_genproducts <- function(x="gffpath") {
df_gff <- rtracklayer::import.gff3(x)
rangesgff <- as(df_gff, "GRanges")
gene_product_df <- as.data.frame(cbind(rangesgff$ID, rangesgff$product) )
rownames(gene_product_df) <- gene_product_df$gene_id
colnames(gene_product_df) <- c("gene_id", "product")

gene_product_df <- gene_product_df %>%
  unite_(., "GeneID", c("gene_id", "product"), sep = "_", remove = F) 

## finally add gene Names

gene_product_df$Name <- rangesgff$Name[match(gene_product_df$gene_id, unique(rangesgff$ID))]

#quality check, remove NAs in column for rownames
gene_product_df_filtered <- gene_product_df %>%
  filter(., !is.na(gene_id))

rownames(gene_product_df_filtered) <- gene_product_df_filtered$gene_id

return(gene_product_df_filtered)
}

####################################################33
# function to get top differential expressed genes

get_heatmap_IDs <- function(dt1, dt2, dt3, topNb=30){
  x_dt <- dt1 %>%
    drop_na() %>%
    mutate(., "gene_id" = rownames(.)) %>%
    select(., gene_id, 1:ncol(.)) %>%
    arrange_(~ desc(log2FoldChange))
  
  y_dt <- dt2 %>%
    drop_na() %>%
    mutate(., "gene_id" = rownames(.)) %>%
    select(., gene_id, 1:ncol(.)) %>%
    arrange_(~ desc(log2FoldChange))
  
  a <- x_dt %>% top_n(wt =log2FoldChange, n =  topNb)
  b <- x_dt %>% top_n(wt =log2FoldChange, n =  -topNb)
  c <- y_dt %>% top_n(wt =log2FoldChange, n =  topNb)
  
  genes <- unique(c(a$gene_id, b$gene_id, c$gene_id))
  final_dt <- subset(dt3, dt3$gene_id %in% genes)
  rownames(final_dt) <- final_dt$gene_id
  return(final_dt)
}



