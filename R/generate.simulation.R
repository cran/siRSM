#################################################################################
# generate.simulation: creates some simulated data
#################################################################################
generate.simulation <- function(
  pc1 = c(1,0.5), pc2 = c(-0.5,1), weight = c(1,0.5),
  rho=0.4,rho1=0.2,err.var=2,
  beta=c(0.5,-0.5,0.5,-0.3,0.6,-0.3),
  n=1000)
{
# Written by Huan Cheng, 2014
# Modified by Mu Zhu, 2014
  uni <- function(v) { c(v[1]/sqrt(v[1]^2+v[2]^2), v[2]/sqrt(v[1]^2+v[2]^2)) }
  pc1 <- uni(pc1)
  pc2 <- uni(pc2)
  v <- cbind(pc1,pc2)
  x11 <- v[1,1]
  x12 <- v[2,1]
  x21 <- v[1,2]
  x22 <- v[2,2]
  a <- (1-rho^2)*x11^2*x12^2
  b <- 2*x11*x12*x21*x22-rho^2*x11^2*x22^2-rho^2*x21^2*x12^2
  c <- x21^2*x22^2-rho^2*x21^2*x22^2
  k <- (-b+sqrt(b^2-4*a*c))/(2*a)
  eigenval1 <- cbind(c(k,0),c(0,1))
  covmatx <- v%*%eigenval1%*%t(v)

  var1 <- covmatx[1,1]
  var2 <- covmatx[2,2]
  var_vector <- c(var1,var2,var1,var2)
  var_mat <- var_vector%*%t(var_vector)
  rho_mat <- matrix(c(1,rho,rho1,rho1,
                    rho,1,rho1,rho1,
                    rho1,rho1,1,rho,
                    rho1,rho1,rho,1),
                  ncol=4,
                  nrow=4)
  cov_mat <- rho_mat*(var_mat^0.5)

  allx<- mvrnorm(n,c(0,0,0,0),Sigma=cov_mat)
  s <- allx[,1:2]
  e <- allx[,3:4]

  s1 <- s %*% uni(weight)
  e1 <- e %*% uni(weight)

  err <- rnorm(n,0,err.var)

  x <- cbind(rep(1,n),s1,e1,s1^2,s1*e1,e1^2)
  y <- x %*% beta + err 

  simu.dat <- as.data.frame(cbind(y,s,e))
  dimnames(simu.dat)[[2]] <- c("intent","s1","s2","e1","e2")
  return(list(data=simu.dat, rho_mat=rho_mat, cov_mat=cov_mat))
}
