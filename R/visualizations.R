showTree <- function(object, type = "effect") 0
setGeneric("showTree")
setMethod("showTree", signature("scDiffPop"), function(object, type = "effect") {
  ###ENSURE PACKAGES ARE INSTALLED
  Tree <- object@tree
  TreeIG <- cbind(Tree[,2], Tree[,1])
  TreeIG <- as.matrix(TreeIG)

  G <- igraph::graph_from_edgelist(TreeIG)
  print(V(G))

  counts <- object@counts
  l2total_counts <- log2(counts[1,1] + counts[1,2])
  V(G)$radius <- as.vector(log2(counts[,1] + counts[,2])/l2total_counts)
  c <- 0.2 #Scaling constant
  ncirc <- round(max(V(G)$radius)/min(V(G)$radius))
  V(G)$radius <- V(G)$radius * c

  xy <- layout_as_tree(G)
  V(G)$x <- xy[, 1]
  V(G)$y <- xy[, 2]

  graph_data <- data.frame(name = V(G)$name, x = V(G)$x, y = V(G)$y, radius = V(G)$radius)
  graph_data$color <- c(0, sign(object@results$effect))
  graph_data$color <- ifelse(graph_data$color > 0, "turquoise", "hotpink1")
  graph_data$color[1] <- "white"
  graph_data$intensity <- c(0, object@results$effect/(max(abs(object@results$effect))))
  graph_data$intensity <- graph_data$intensity

  for(i in 2:length(V(G)$name)) {
    s <- out@results$padj[i-1]
    print(s)
    if(s < 0.001) {
      V(G)$name[i] <- paste(graph_data$name[i], "***")
    }
    else if(s < 0.01) {
      V(G)$name[i] <- paste(graph_data$name[i], "**")
    }
    else if(s < 0.05) {
      V(G)$name[i] <- paste(graph_data$name[i], "*")
    }
  }

  if("effect" %in% type) {
    p <- ggraph(G, "manual", x= V(G)$x, y=V(G)$y) + geom_edge_link()
    p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius), colour = "black", show.legend = FALSE, data = graph_data, fill="white")
    p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius, fill = intensity), colour = NA, show.legend = TRUE, data = graph_data)
    #p <- p + scale_fill_manual(values = graph_data$color)
    p <- p + scale_fill_gradient2(low = "turquoise", mid = "white", high = "hotpink1", guide = "colourbar",
                                  name = "", breaks = c(-1, 1), labels = c(object@meta.data$phenotypes[1], object@meta.data$phenotypes[2]), limits = c(-1,1))
    #p <- p + scale_alpha_manual(values = graph_data$intensity)
    p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
    p <- p + theme_graph()
    plot(p)
  }
  if("pies" %in% type) {
    graph_data$pht1 <- counts[,1]; graph_data$pht2 <- counts[,2]
    p <- ggraph(G, "manual", x=  V(G)$x, y=V(G)$y)+ geom_edge_link()
    p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius), colour = NA, show.legend = FALSE, data = graph_data, fill="white")
    p <- p + geom_scatterpie(
      aes(x=x, y=y, r=radius, alpha = forcats::fct_inorder(name)),
      data = graph_data ,
      cols = c("pht1", "pht2"),
      colour = NA,
      legend_name = "Phenotype",
    ) + scale_alpha_manual(values = c(0.3, 1 - object@results$padj), name = NULL, labels = NULL)
    p <- p + scale_fill_manual(values = c("blue", "red"), labels = c("pht1", "pht2"))
    p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
    p <- p + theme_graph()
    plot(p)
  }

})

showMarkerPlots <- function(object, nrow = 3, ncol = 3) 0
setGeneric("showMarkerPlots")
setMethod("showMarkerPlots", signature("scDiffPop"), function(object, nrow = 3, ncol = 3) {
  par(mfrow = c(nrow, ncol))
  nullList <- sapply(out@markers, function(genes) {
    plot(genes$x,genes$y,xlab="Marker Strength", ylab = "Phenotype stat", col = "white", main = genes$main)
    text(genes$x,genes$y,labels=names(genes$x), cex = 0.5)
    abline(lm(genes$y ~ genes$x), col = "blue", lty = 3)
  })
})




