# Create big.matrix
BM <- as.big.matrix(matrix(seq(12), 3, 4))
head(BM)
# Call in-place tranpose
iptBM(BM, direction="r2c")
head(BM)
