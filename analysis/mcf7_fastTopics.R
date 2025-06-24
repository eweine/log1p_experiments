# Run this script only after first running the code in mcf7.Rmd.
library(fastTopics)

# First initialize the fit using Eric's code.
set.seed(1)
# scale_cols <- function (A, b)
#   t(t(A) * b)
# fit_init <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e4,
#                                 loglik = "exact",init_method = "rank1",
#                                 control = list(maxiter = 4,verbose = TRUE))
s <- rowSums(counts)
n <- nrow(counts)
L <- matrix(0,n,4)
for (i in 1:n) {
  x <- runif(3)
  x <- x/sum(x)
  L[i,] <- s[i]/2*c(1,x)
}

# I'm saving the data and initial setup (the initial estimate of L) to
# use used in a separate analysis which I will use to illustrate EM
# vs. SCD for Poisson NMF.
library(tools)
save(list = c("samples","genes","counts","L"),
     file = "mcf7.RData")
resaveRdaFiles("mcf7.RData")

# Fit the topic model by performing either (i) 80 EM updates or (ii)
# 80 SCD updates. Both of the fits are first initialized by running 20
# EM updates.
control <- list(extrapolate = FALSE,numiter = 4,nc = 8)
fit_init <- init_poisson_nmf(counts,L = L,init.method = "random")
fit0 <- fit_poisson_nmf(counts,fit0 = fit_init,numiter = 4,method = "em",
                        control = control)
fit1 <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 100,method = "em",
                        control = control)
fit2 <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 100,method = "scd",
                        control = control)
fit3 <- fit_poisson_nmf(counts,fit0 = fit0,numiter = 50,method = "scd",
                        control = control)
control$extrapolate <- TRUE
fit3 <- fit_poisson_nmf(counts,fit0 = fit3,numiter = 50,method = "scd",
                        control = control)
fit4 <- fit_poisson_nmf(counts,fit0 = fit3,numiter = 500,method = "scd",
                        control = control)

# Compare the fits in Structure plots.
topic_colors <- c("olivedrab","dodgerblue","darkblue","tomato")
p0 <- structure_plot(fit0,grouping = samples$label,topics = 1:4,
                     colors = topic_colors) +
  ggtitle("init")
p1 <- structure_plot(fit1,grouping = samples$label,topics = 1:4,
                     colors = topic_colors) +
  ggtitle("EM")
p2 <- structure_plot(fit2,grouping = samples$label,topics = 1:4,
                     colors = topic_colors) +
  ggtitle("SCD")
p3 <- structure_plot(fit3,grouping = samples$label,topics = 1:4,
                     colors = topic_colors) +
  ggtitle("SCD + extrapolate")
p4 <- structure_plot(fit3,grouping = samples$label,topics = 1:4,
                     colors = topic_colors,loadings_order = 1:n) +
  ggtitle("best")
print(plot_grid(p0,p1,p2,p3,p4,nrow = 5,ncol = 1))

# Compare the log-likelihoods.
logliks <- c("initial" = sum(loglik_multinom_topic_model(counts,fit0)),
             "em"      = sum(loglik_multinom_topic_model(counts,fit1)),
             "scd"     = sum(loglik_multinom_topic_model(counts,fit2)),
             "scd+ex"  = sum(loglik_multinom_topic_model(counts,fit3)),
             "best"    = sum(loglik_multinom_topic_model(counts,fit4)))
print(logliks)
n <- nrow(counts)
print(logliks/n)
