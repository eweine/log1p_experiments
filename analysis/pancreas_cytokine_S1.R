library(Matrix)
library(fastTopics)
library(log1pNMF)
library(ggplot2)
library(cowplot)
set.seed(1)
load("../data/pancreas_cytokine.RData")
clusters <- read.csv("../data/pancreas_cytokine_S1_tm_k=13_clusters.csv",
                     stringsAsFactors = FALSE)

# Here we will analyze the cells from the untreated mouse only.
i       <- which(samples$mouse == "S1")
samples <- samples[i,]
counts  <- counts[i,]

# Remove three cells that appear to be outliers (one of them appears
# to be an acinar cell based on Eric's analysis).
outliers <- c("TTTGTTGTCGTTAGTG-1","TTTGTTGGTAGAGCTG-1","CCCAACTCACTCATAG-1")
i        <- which(!is.element(samples$barcode,outliers))
samples  <- samples[i,]
counts   <- counts[i,]

# Remove genes that are expressed in fewer than 5 cells.
j      <- which(colSums(counts > 0) > 4)
genes  <- genes[j,]
counts <- counts[,j]

# Get the cluster assignments from the topic model.
samples$cluster <- factor(clusters$cluster)

# Fit a Poisson log1p NMF to the counts, with K = 13.
set.seed(1)
#
# IDEA: Use NNLM to initialize the fit.
#
fit <- fit_poisson_log1p_nmf(counts,K = 13,cc = 1,
                             loglik = "approx",init_method = "rank1",
                             control = list(maxiter = 1000,verbose = TRUE))

# Compare the factors to the clusters.
scale_cols <- function (A, b)
  t(t(A) * b)
set.seed(1)
L <- fit$LL
d <- apply(L,2,max)
L <- scale_cols(L,1/d)
clusters <- samples$cluster
clusters <- factor(clusters,c("beta","alpha","delta+epsilon","gamma","duct",
                              "endothelial-mesenchymal","macrophage"))
i <- c(sample(which(clusters == "beta"),400),
       which(clusters != "beta"))
p1 <- structure_plot(L[i,],grouping = clusters[i],gap = 10,n = Inf,
                     topics = c(1,2,4,5,9,10,11,12,13))
p2 <- structure_plot(L[i,],grouping = clusters[i],gap = 10,n = Inf,
                     topics = c(3,6,7,8))
plot_grid(p1,p2,nrow = 2,ncol = 1)

# TO DO:
#
# * Re-run fastTopics with k = 13.
#
# * Develop this into a workflowr analysis,
#   pancreas_cytokine_S1.Rmd.
#
library(NNLM)
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
shifted_log_counts <- log1p(counts/(a*s))
rownames(shifted_log_counts) <- NULL
shifted_log_counts_dense <- as.matrix(shifted_log_counts)
set.seed(1)
nmf <- nnmf(shifted_log_counts_dense,k = 13,loss = "mse",method = "scd",
            max.iter = 100,verbose = 2,n.threads = 8)
W <- nmf$W
d <- apply(W,2,max)
W <- scale_cols(W,1/d)
p <- structure_plot(W[i,],grouping = clusters[i],gap = 10,n = Inf)

fit2 <- fit_poisson_log1p_nmf(counts,K = 13,cc = 1,
                              init_LL = nmf$W,
                              init_FF = t(nmf$H),
                              loglik = "exact",
                              control = list(maxiter = 20,verbose = TRUE))

ll1 <- logLik(fit,counts)
ll2 <- logLik(fit2,counts)
