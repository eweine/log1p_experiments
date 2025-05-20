# GLM-PCA
# TO DO

# Additive vs. multiplicative.
LF <- tcrossprod(fit$L,fit$F)
err_tm <- apply(abs(counts - LF),2,median)
res <- ldf(fl_nmf,type = "i")
LF <- with(res,L %*% diag(D) %*% t(F))
err_nmf <- apply(abs(Y - LF),2,median)
plot(log10(err_tm),err_nmf,pch = 20)
