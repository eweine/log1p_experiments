# Additive vs. multiplicative.
LF <- tcrossprod(fit$L,fit$F)
mu0 <- colSums(counts)/sum(s)
err_tm <- apply(abs(counts - LF),2,median)
res <- ldf(fl_nmf,type = "i")
LF <- with(res,L %*% diag(D) %*% t(F))
err_nmf <- apply(abs(Y - LF),2,median)
i <- which(mu0 > 100 & mu0 < 500)
plot(log10(err_tm[i]),err_nmf[i],pch = 20)

# NEXT: Look at the results for selected genes.

# Compute the mean expression for each of the four conditions.
n  <- ncol(counts)
mu <- data.frame(EtOH = rep(0,n),
                 RA   = rep(0,n),
                 TGFb = rep(0,n),
                 "RA+TGFb" = rep(0,n),
                 check.names = FALSE)
for (treatment in levels(samples$label)) {
  i <- which(samples$label == treatment)
  mu[[treatment]] <- colSums(counts[i,])/sum(s[i])
}
rownames(mu) <- colnames(counts)

i <- 9172
i <- 24365
genes[i,1:2]
mu[i,]
plot(c(as.matrix(mu[i,])),pch = 20)

