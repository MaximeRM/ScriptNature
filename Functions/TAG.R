##%######################################################%##
#                                                          #
####                        TAG                         ####
#                                                          #
##%######################################################%##

TAG <- function(S, R) {
  tag <- function(X, R) {
    as.matrix(R[, !is.na(X)]) %*% 
      as.vector(X[!is.na(X)])/as.matrix(R[, !is.na(X)]) %*% 
      rep(1, sum(!is.na(X)))
  }
  
  abo <- function(X, R) {
    round(apply(R[, !is.na(X)], 1, sum)/apply(R, 1, sum)*100, 1)
  }
  
  T <- apply(S, 2, tag, R)
  A <- apply(S, 2, abo, R)
  rownames(T) <- rownames(A) <- rownames(R)
  
  return(list(tag=T, abo=A))
}
