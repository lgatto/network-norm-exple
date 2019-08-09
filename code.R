library("gridGraphics")
library("gridExtra")
library("directlabels")
library("ggplot2")

grab_grob <- function(){
  grid.echo()
  grid.grab()
}


## https://www.biostars.org/p/285296/

library("igraph")
library("bioDist")

make_data <- function(n_genes) {
    set.seed(42)
    n_samples <- 5
    e <- matrix(runif(n_samples * n_genes, min = 0, max = 12),
                ncol = n_samples)
    colnames(e) <- paste0("sample_", seq_len(n_samples))
    rownames(e) <- seq_len(n_genes)
    ## add control genes
    a <- c(1, 6, 2, 4, 7)
    b <- a + 4 + rnorm(n_samples, 0, 0.1)
    c <- b / 3 + c(1, 2, -1, 2, 0)
    m <- rbind(e, a, b, c)
    ## rownames(m) <- c(paste0("r", 1:n_genes), letters[1:3])
    rownames(m) <- 1:(n_genes + 3)
    m
}

make_graph <- function(x, w = NULL) {
    x <- as.matrix(x)
    g <- graph.adjacency(x,
                         mode = "undirected",
                         weight = TRUE,
                         diag = FALSE)
    ## if (is.null(w))
    ##     w <- quantile(E(g)$weight, 0.9)
    ## g <- delete_edges(g, E(g)[which(E(g)$weight > w)])
    V(g)$color <- c(rep("lightgrey", nrow(x) - 3),
                    "steelblue",
                    "steelblue",
                    "lightblue")
    g
}

plot_graph <- function(g, main = "") {
    w  <- E(g)$weight / max(E(g)$weight)
    w <- (1 - w) * 5
    plot(g, main = main,
         edge.width = w,
         edge.label = rep("", 10),
         layout = layout.kamada.kawai)
}



## gene networks
e <- make_data(n_genes = 2)
g_euc <- make_graph(e_euc <- euc(e))
g_scaled_euc <- make_graph(e_scaled_euc <- euc(t(scale(t(e)))))
g_cor <- make_graph(e_cor <- cor.dist(e))


## e <- make_data(n_genes = 20)
## e2 <- t(scale(t(e)))
## g_10_cor <- make_graph(cor.dist(e))
## g_10_euc <- make_graph(euc(e))
## g_10_scaled_euc <- make_graph(euc(e2))

## e <- make_data(n_genes = 100)
## e2 <- t(scale(t(e)))
## g_100_cor <- make_graph(cor.dist(e))
## g_100_euc <- make_graph(euc(e))

## g_100_scaled_euc <- make_graph(euc(e2))


edf <- data.frame(expression = as.numeric(e),
                  sample = rep(1:5, each = 5),
                  gene = rep(1:5, 5),
                  type = rep(c("random", "random",
                               "co-expressed", "co-expressed",
                               "expressed"), 5))

## matplot(t(e), type = "b", lwd = 2, lty = 1,
##         ylab = "Expression data",
##         col = c(rep("lightgrey", 2), "steelblue",
##                 "steelblue", "lightblue"))

gr_e <-
    ggplot(edf, aes(x = sample, y = expression, group = gene, colour = type)) +
    geom_line(lwd = 1) +
    scale_colour_manual(values = c("random" = "darkgrey",
                                   "co-expressed" = "steelblue",
                                   "expressed" = "lightblue")) +
    geom_dl(aes(label = gene), method = "last.points") +
    theme_bw()

## plot



plot_graph(g_euc)
gr_euc <- grab_grob()
plot_graph(g_scaled_euc)
gr_scaled_euc <- grab_grob()
plot_graph(g_cor)
gr_cor <- grab_grob()
heatmap(as.matrix(e_cor), scale = "none")
gr_h_cor <- grab_grob()
heatmap(as.matrix(e_euc), scale = "none")
gr_h_euc <- grab_grob()
heatmap(as.matrix(e_scaled_euc), scale = "none")
gr_h_scaled_euc <- grab_grob()


lm <- matrix(c(1, 1:7), ncol = 2, byrow = TRUE)

pdf("scaling.pdf", width = 7, height = 14)
grid.arrange(gr_e,
             gr_h_euc, gr_euc,
             gr_h_scaled_euc, gr_scaled_euc,
             gr_h_cor, gr_cor,
             layout_matrix = lm)
dev.off()