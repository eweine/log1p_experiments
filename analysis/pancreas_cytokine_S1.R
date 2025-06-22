# Compare the factors to the clusters.
p1 <- structure_plot(L[i,],grouping = clusters[i],gap = 10,n = Inf,
                     topics = c(1,2,4,5,9,10,11,12,13))
p2 <- structure_plot(L[i,],grouping = clusters[i],gap = 10,n = Inf,
                     topics = c(3,6,7,8))
plot_grid(p1,p2,nrow = 2,ncol = 1)
