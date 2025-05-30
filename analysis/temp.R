fit_log1p_c10 <- 
  fit_poisson_log1p_nmf(counts,K = 4,cc = 10,loglik = "exact",
                        init_method = "random",
                        control = list(maxiter = 200,verbose = TRUE))
L <- fit_log1p_c10$LL
d <- apply(L,2,max)
L <- scale_cols(L,1/d)
p <- structure_plot(L,grouping = samples$label,topics = 1:4,
                    colors = topic_colors)
p
