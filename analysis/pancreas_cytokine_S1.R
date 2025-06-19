library(Matrix)
library(fastTopics)
library(log1pNMF)
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

# Fit a Poisson log1p NMF to the counts, with K = 11.
set.seed(1)
fit <- fit_poisson_log1p_nmf(counts,K = 11,cc = 1,
                             loglik = "approx",init_method = "rank1",
                             control = list(maxiter = 100,verbose = TRUE))

# Compare the factors to the clusters.
scale_cols <- function (A, b)
  t(t(A) * b)
L <- fit$LL
d <- apply(L,2,max)
L <- scale_cols(L,1/d)
p <- structure_plot(L,grouping = samples$cluster,gap = 10,n = Inf)
