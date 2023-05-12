## code to prepare `ord1` dataset goes here
t <- c(0, 1, 2, 3, 4, 5)
clist <- data.frame(matrix(nrow = length(t)))
for (i in 1:3) {
  set.seed(i)
  err <- rnorm(length(t), 0, 0.05)
  c1 <- 15 * exp(-t * 0.7) * (1 + err)
  clist <- cbind.data.frame(clist, c1)
}
concentration <- rowMeans(as.matrix(clist[-1]))
std.error <- apply(as.matrix(clist[-1]), 1, sd) / sqrt(i)

ord1 <- data.frame(t, concentration, std.error)




usethis::use_data(ord1, overwrite = TRUE)
