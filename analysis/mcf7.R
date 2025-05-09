#
# Downloaded the files from GEO (accession GSE152749).
#
library(data.table)
library(GEOquery)
library(Matrix)
library(fastTopics)
library(flashier)
set.seed(1)

# Load the gene information.
genes <- fread("../data/Human.GRCh38.p13.annot.tsv.gz",
               sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(genes) <- "data.frame"

# Load the RNA-seq data.
counts <- fread("../data/GSE152749_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(counts) <- "data.frame"
rownames(counts) <- counts$GeneID
counts <- counts[,-1]
counts <- as.matrix(counts)
storage.mode(counts) <- "double"
counts <- t(counts)
ids <- rownames(counts)

# Load the sample information.
#
# EtOH = ethanol
# RA   = retinoic acid
# TGFb = TGF-beta
#
geo <- getGEO(filename = "../data/GSE152749_family.soft.gz")
samples <- data.frame(id = names(GSMList(geo)),
                      treatment = sapply(GSMList(geo),
                                         function (x) Meta(x)$title))

samples <- samples[ids,]
rownames(samples) <- NULL
samples <- transform(samples,
                     EtOH = grepl("EtOH",treatment,fixed = TRUE),
                     RA   = grepl("RA",treatment,fixed = TRUE),
                     TGFb = grepl("TGFb",treatment,fixed = TRUE))
samples$label                 <- "EtOH"
samples[samples$RA,"label"]   <- "RA"
samples[samples$TGFb,"label"] <- "TGFb"
samples[with(samples,RA & TGFb),"label"] <- "RA+TGFb"
samples <- transform(samples,
                     label = factor(label,c("EtOH","RA","TGFb","RA+TGFb")))

# Remove genes that are expressed in fewer than 4 samples.
x <- colSums(counts > 0)
i <- which(x > 3)
genes  <- genes[i,]
counts <- counts[,i]

# Fit a topic model to the counts.
fit0 <- fit_poisson_nmf(counts,k = 4,init.method = "random",numiter = 50,
                       verbose = "detailed",
                       control = list(nc = 4,extrapolate = FALSE))
fit <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 50,verbose = "detailed",
                       control = list(nc = 4,extrapolate = TRUE))
fit <- poisson2multinom(fit)
res <- cbind(samples,round(fit$L,digits = 3))
rownames(res) <- NULL

i <- order(res$k1,decreasing = TRUE)
res[i,c("k1","treatment")]

i <- order(res$k2,decreasing = TRUE)
res[i,c("k2","treatment")]

i <- order(res$k3,decreasing = TRUE)
res[i,c("k3","treatment")]

i <- order(res$k4,decreasing = TRUE)
res[i,c("k4","treatment")]

x <- rowSums(fit$L > 0.1)
temp <- cbind(x,samples$treatment)
temp[order(x),]

topic_colors <- c("tomato","darkblue","gold","skyblue")
p <- structure_plot(fit,grouping = samples$label,
                    topics = 1:4       ,colors = topic_colors)
print(p)

# Compute the shifted log counts.
a <- 1
s <- rowSums(counts)
s <- s/mean(s)
Y <- log1p(counts/(a*s))

# Fit an NMF to the shifted log counts.
n <- nrow(counts)
x <- rpois(1e7,1/n)
s1 <- sd(log(x + 1))
fl_nmf <- flash(Y,ebnm_fn = ebnm_point_exponential,var_type = 2, 
                greedy_Kmax = 4,S = s1,backfit = TRUE)
res <- ldf(fl_nmf,type = "i")
p <- structure_plot(res$L,grouping = samples$label,
                    colors = topic_colors)
print(p)
