# Independence, CS and AR(1) corr matrix generating function
# It will be a m x m matrix
corrmat <- function(alpha=0, m, type="Independence"){
  if (type == "Independence") {
    R <- diag(m)
  } else if (type == "CS") {
    one <- rep(1,times=m)
    R <- (1-alpha)*diag(m) + alpha*tcrossprod(one)
  } else if (type == "AR1") {
    times <- c(1:m)
    H <- abs(outer(times, times, "-"))
    R <- alpha^H
    R[cbind(1:m, 1:m)] <- R[cbind(1:m, 1:m)]
    return(R)
  } else {
    stop('Error in corrmat(): Choose type either "Independence", "CS" or "AR1"')
  }
}