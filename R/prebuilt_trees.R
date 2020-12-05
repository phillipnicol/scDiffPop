
makeImmuneTree <- function() {
  Immune_edgelist <- c(0, "Lymphocytes", "Lymphocytes", "B",
                       "Lymphocytes", "NK", "Lymphocytes", "T")
  return(matrix(Immune_edgelist, ncol = 2, byrow = TRUE))
}
