# Run this script only after first running the code in mcf7.Rmd.
library(log1pNMF)
library(NNLM)
library(fastTopics)

# Compare the log1pNMF and NNLM optimization algorithms.
set.seed(1)
fit1 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1,
                              loglik = "exact",init_method = "rank1",
                              control = list(maxiter = 100,verbose = TRUE))

nmf <- nnmf(Y,k = 4,init = list(W = fit1$LL,H = t(fit1$F)),
            loss = "mse",method = "scd",max.iter = 20,
            verbose = 2,n.threads = 4)

fit1 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1,loglik = "exact",
                              init_LL = fit1$LL,init_FF = fit1$FF,
                              control = list(maxiter = 30,verbose = TRUE))

fit2 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1,loglik = "exact",
                              init_LL = nmf$W,init_FF = t(nmf$H),
                              control = list(maxiter = 10,verbose = TRUE))

# Compare the two model fits.
print(diff(fit1$objective_trace))
print(diff(fit2$objective_trace))
print(logLik(fit1,counts) - logLik(fit2,counts))

# Compare the memberships in a scatterplot.
scale_cols <- function (A, b)
  t(t(A) * b)
L  <- fit1$LL
d  <- apply(L,2,max)
L1 <- scale_cols(L,1/d)
L  <- fit2$LL
d  <- apply(L,2,max)
L2 <- scale_cols(L,1/d)
plot(L1,L2,pch = 20,xlab = "fit1",ylab = "fit2")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# Compare the log1pNMF and fastTopics optimization algorithms.
set.seed(1)
fit1 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e4,
                              loglik = "exact",init_method = "rank1",
                              control = list(maxiter = 100,verbose = TRUE))

tm <- init_poisson_nmf(counts,F = fit1$FF,L = fit1$LL)
tm <- fit_poisson_nmf(counts,fit0 = tm,numiter = 20,verbose = "detailed",
                      control = list(nc = 4,extrapolate = TRUE))

fit1 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e4,loglik = "exact",
                              init_LL = fit1$LL,init_FF = fit1$FF,
                              control = list(maxiter = 30,verbose = TRUE))

fit2 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e4,loglik = "exact",
                              init_LL = tm$L,init_FF = tm$F,
                              control = list(maxiter = 10,verbose = TRUE))

# Compare the two model fits.
print(diff(fit1$objective_trace))
print(diff(fit2$objective_trace))
print(logLik(fit1,counts) - logLik(fit2,counts))

# Visualize the model fits using Structure plots.
topic_colors <- c("olivedrab","dodgerblue","darkblue","tomato")
L  <- fit1$LL
d  <- apply(L,2,max)
L1 <- scale_cols(L,1/d)
p1 <- structure_plot(L,grouping = samples$label,topics = 1:4,
                     colors = topic_colors)
L  <- fit2$LL
d  <- apply(L,2,max)
L2 <- scale_cols(L,1/d)
p2 <- structure_plot(L,grouping = samples$label,topics = 1:4,
                     colors = topic_colors)
print(plot_grid(p1,p2,nrow = 2,ncol = 1,labels = c("fit1","fit2")))

# Compare the memberships in a scatterplot.
plot(L1,L2,pch = 20,xlab = "fit1",ylab = "fit2")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

