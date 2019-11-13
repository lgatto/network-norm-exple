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
e_euc <- euc(e)
g_euc <- make_graph(e_euc)
e_scaled_euc <- euc(t(scale(t(e))))
e_scaled_euc[8] <- 0.25 ## to avoid overlapping nodes
g_scaled_euc <- make_graph(e_scaled_euc)
e_cor <- cor.dist(e)
e_cor[8] <- 0.08 ## to avoid overlapping nodes
g_cor <- make_graph(e_cor)


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
    geom_dl(aes(label = gene, size = 5), method = list("last.points", cex = 2)) +
    theme_bw() + 
    labs(title = "(a)") +
    theme(plot.title = element_text(size = 16)) 


## save figures separately
ggsave("subfig-a.pdf", gr_e, width = 8, height = 4.5)

pdf("subfig-b.pdf")
heatmap(as.matrix(e_euc), scale = "none")
dev.off()

pdf("subfig-d.pdf")
heatmap(as.matrix(e_scaled_euc), scale = "none")
dev.off()

pdf("subfig-f.pdf")
heatmap(as.matrix(e_cor), scale = "none")
dev.off()

ggsave("subfig-c.pdf", plot_graph(g_euc))
ggsave("subfig-e.pdf", plot_graph(g_scaled_euc))
ggsave("subfig-g.pdf", plot_graph(g_cor))

## compose plots (original file)

## heatmap(as.matrix(e_euc), scale = "none")
## text(-1.2, 1.2, "(b)", cex = 1.5)
## gr_h_euc <- grab_grob()
## plot_graph(g_euc)
## text(-1, 1, "(c)", cex = 1.5)
## gr_euc <- grab_grob()
## heatmap(as.matrix(e_scaled_euc), scale = "none")
## text(-1.2, 1.2, "(d)", cex = 1.5)
## gr_h_scaled_euc <- grab_grob()
## plot_graph(g_scaled_euc)
## text(-1, 1, "(e)", cex = 1.5)
## gr_scaled_euc <- grab_grob()
## heatmap(as.matrix(e_cor), scale = "none")
## text(-1.2, 1.2, "(f)", cex = 1.5)
## gr_h_cor <- grab_grob()
## plot_graph(g_cor)
## text(-1.2, 1.3, "(g)", cex = 1.5)
## gr_cor <- grab_grob()


## lm <- matrix(c(1, 1:7), ncol = 2, byrow = TRUE)

## pdf("scaling.pdf", width = 7, height = 14)
## grid.arrange(gr_e,
##              gr_h_euc, gr_euc,
##              gr_h_scaled_euc, gr_scaled_euc,
##              gr_h_cor, gr_cor,
##              layout_matrix = lm)
## dev.off()
