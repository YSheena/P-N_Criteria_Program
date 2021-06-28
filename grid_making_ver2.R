### Making objects "M2", "M3", "M4" ##################
# these objects are only used for the moments or cumulants caluation

# "fSupSymIncreasingFirstR " (shortened to "a_v")
# See "Multidimensional Array Indexing and Storage" by Jan de Leeuw, 2017
# DOI:10.13140/RG.2.2.23236.73604
# Input: vector (i1,i2,...,ik), i1 <= i2 <= ,..., <= ik, 
# the index of an element of a symmetric array
# Output: The location in the vector of the input array element, 
# when a symmetric array is vectorized.
a_v <- function (cell) { 
  f <- function (r){choose (r + (cell[r] - 1) - 1, r)} 
  rank <- length (cell)
  cell <- sort (cell)
  return (1 + sum (sapply (1:rank, f)))
}

# the inverse function of the above function
fSupSymIncreasingFirstInverseR <- function (dimension, rank, index) { 
  last.true <- function (x) {
    w <- which (x)
    return (w[length(w)])
  }
  a <- rep (0, rank) 
  v <- index - 1
  for (k in rank:1) {
    s <- choose (k + (0:(dimension - 1)) - 1, k) 
    u <- last.true (v >= s)
    a[k] <- u
    v <- v - s[u]
  }
  return (a)
}


grid_making <- function(p) {

  # Specialized "fSupSymIncreasingFirstInverse" 
  # for a specific "dimension", "rank"
  v_a2 <- function(x) {
    fSupSymIncreasingFirstInverseR(c(p,p),2,x)
  }
  v_a3 <- function(x) {
    fSupSymIncreasingFirstInverseR(c(p,p,p),3,x)
  }
  v_a4 <- function(x) {
    fSupSymIncreasingFirstInverseR(c(p,p,p,p),4,x)
  }
  
  # The size of the vectorizaion of the array
  N2 <- a_v(c(p,p))
  N3 <- a_v(c(p,p,p))
  N4 <- a_v(c(p,p,p,p))
  
  # The matrix whose each row is the 2,3,4-dimensional symmetric array index
  M2 <- lapply(1:N2,v_a2)
  M3 <- lapply(1:N3,v_a3)
  M4 <- lapply(1:N4,v_a4) # takes several minutes

  
  return(list(M2 = M2, M3 = M3, M4 = M4))
}
