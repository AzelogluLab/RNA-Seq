QCftcountsSample <- function(counts){
  ##Quality Check Counts Control
  ##Description:
  ##A quality check function that maskes a heatmap from a counts tabl
  ##Arguments:
  ##counts: A counts table with the rows as genes and the columns as samples
  ##Returns:
  ##heatmap : a heat map comparing the individual group samples to each other
  library(ggpubr)
  library(RColorBrewer)
  library(gplots)
  Pearson <- cor(counts, counts, method = "pearson", use = "complete.obs")
  par(oma = c(4,2,4,4))
  heatmap <- heatmap.2(Pearson, dendrogram = 'col', col = brewer.pal(n = 11, name = "RdBu"), trace = "none", symkey = F)
  return(heatmap)
}
