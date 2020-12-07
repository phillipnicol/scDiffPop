showTree <- function(object) 0
setGeneric("showTree")
setMethod("showTree", signature("scDiffPop"), function(object) {
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
  graph_data$intensity <- c(0,(abs(object@results$effect)/(max(abs(object@results$effect)))))
  print(graph_data)

  p <- ggraph(G, "manual", x= V(G)$x, y=V(G)$y) + geom_edge_link()
  p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius), colour = "black", show.legend = FALSE, data = graph_data, fill="white")
  p <- p + geom_node_circle(aes(x0=x,y0=y,r=radius, fill = forcats::fct_inorder(name), alpha = forcats::fct_inorder(name)), colour = NA, show.legend = FALSE, data = graph_data)
  p <- p + scale_fill_manual(values = graph_data$color)
  p <- p + scale_alpha_manual(values = graph_data$intensity)
  p <- p + geom_node_label(aes(label = name, angle = 90), repel = FALSE, nudge_y = 0.25, col = "midnightblue")
  p <- p + theme_graph()
  plot(p)
})
