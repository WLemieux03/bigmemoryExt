# Create big.matrix
A <- matrix(seq(12), 3, 4)
BM <- as.big.matrix(A)
C <- as.big.matrix(A)
head(BM)
# Call in-place cbind
# cbind only 1st and 3rd column
cbindBMIP(BM, C, c(1,3))
head(BM)
