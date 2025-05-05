library(Matrix)
library(passPCA)
load("../data/pancreas.RData")
set.seed(1)

# Select the CEL-seq2 data (Muraro et al, 2016).
# This should select 2,285 cells.
i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]
sample_info <- transform(sample_info,celltype = factor(celltype))

# Remove genes that are expressed in fewer than 10 cells.
j      <- which(colSums(counts) > 0 & colSums(counts > 0) > 4)
counts <- counts[,j]

# Load the log1p-NMF model that was previously fit to these data.
load("~/git/passPCA/inst/paper_figures/data/experiment_results.Rdata")
fit_k9 <- res_list$pancreas[["0.001"]]
range(cor(fit_k9$FF))
rm(res_list)

stop()

# Compute the size factors.
s <- rowSums(counts)
s <- s / mean(s)

# Fit a Poisson NMF with the shifted log count link, with c = 0.01.
cc <- 0.001
fit_k9 <- fit_poisson_log1p_nmf(counts,K = 9,cc = cc,s = cc * s,
                                loglik = "exact",init_method = "frob_nmf",
                                control = list(threads = 1,maxiter = 100))

fit_k9 <- fit_poisson_log1p_nmf(counts,K = 9,cc = cc,s = cc * s,
                                loglik = "exact",init_method = "frob_nmf",
                                control = list(threads = 1,maxiter = 100))

j <- which(colMeans(counts) > 10)
cor(fit_k9$FF[j,])
