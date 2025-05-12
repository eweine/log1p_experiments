# PCA
library(rsvd)
pca <- rpca(Y,k = 4,center = TRUE,scale = FALSE)
colnames(pca$x) <- paste0("PC",1:4)
pdat <- data.frame(samples,pca$x)
ggplot(pdat,aes(x = PC1,y = PC2,color = label)) +
  geom_point() +
  theme_cowplot(font_size = 10)
  
# GLM-PCA
# TO DO.

# Additive vs. multiplicative.
LF <- tcrossprod(fit$L,fit$F)
err_tm <- apply(abs(counts - LF),2,median)
res <- ldf(fl_nmf,type = "i")
LF <- with(res,L %*% diag(D) %*% t(F))
err_nmf <- apply(abs(Y - LF),2,median)
plot(log10(err_tm),err_nmf,pch = 20)
