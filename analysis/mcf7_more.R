# Run this script only after first running the code in mcf7.Rmd.
library(log1pNMF)
library(NNLM)
library(fastTopics)

# Compare log1pNMF and NNLM optimization.
set.seed(1)
fit1 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1,
                              loglik = "exact",init_method = "rank1",
                              control = list(maxiter = 100,verbose = TRUE))

nmf <- nnmf(Y,k = 4,init = list(W = fit1$LL,H = t(fit1$F)),
            loss = "mse",method = "scd",max.iter = 20,
            verbose = 2,n.threads = 4)

fit2 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1,loglik = "exact",
                              init_LL = nmf$W,init_FF = t(nmf$H),
                              control = list(maxiter = 10,verbose = TRUE))

fit1 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1,loglik = "exact",
                              init_LL = fit1$LL,init_FF = fit2$FF,
                              control = list(maxiter = 10,verbose = TRUE))

print(diff(fit1$objective_trace))
print(diff(fit2$objective_trace))
print(logLik(fit1,counts))
print(logLik(fit2,counts))

# Compare log1pNMF and fastTopics optimization.
set.seed(1)
fit1 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e-4,
                              loglik = "exact",init_method = "rank1",
                              control = list(maxiter = 10,verbose = TRUE))
print(any(is.na(fit1$objective_trace)))
fit2 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e-4,loglik = "exact",
                              init_LL = fit1$LL,init_FF = fit2$FF,
                              control = list(maxiter = 10,verbose = TRUE))
print(any(is.na(fit2$objective_trace)))

tm <- init_poisson_nmf(counts,F = fit1$FF,L = fit1$LL)
tm <- fit_poisson_nmf(counts,fit0 = tm,numiter = 10,verbose = "detailed",
                      control = list(nc = 4,extrapolate = TRUE))

fit2 <- fit_poisson_log1p_nmf(counts,K = 4,cc = 1e-4,loglik = "exact",
                              init_LL = tm$L,init_FF = tm$F,
                              control = list(maxiter = 10,verbose = TRUE))

print(diff(fit1$objective_trace))
print(diff(fit2$objective_trace))
print(logLik(fit1,counts))
print(logLik(fit2,counts))
