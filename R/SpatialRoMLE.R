#Spatial Robust MLE
RoMLE.error <- function(initial.beta,  initial.s2, initial.lambda, W, y, x, phi.function, converge.v, iter, print.values) {

  initial.par <- data.frame(t(initial.beta), initial.s2, initial.lambda)

  n <- NROW(x)
  m <- NCOL(x)

  X <- cbind(1,x)
  X <- as.matrix(X)
  I <- diag(1,n,n)
  W <- as.matrix(W)

  b <- initial.beta
  s2 <- initial.s2
  l <- initial.lambda

  b.prev <- matrix(c(b), nrow=m+1, ncol=1)
  sl.prev <- matrix(c(s2,l), nrow=2, ncol=1)

  if (phi.function==1) {
    #Cauchy
    c <- 2.385
    fi  <- function(r) (r/(1+(r/c)^2))
    bfi <- function(r) ((c^4-r^2*c^2)/(c^2+r^2)^2)
    #K interal
    integrand  <- function(r) (r/(1+(r/c)^2))^2 * (2*pi)^-0.5*exp(-r^2/2)
    K <- integrate(integrand, lower = -Inf, upper = Inf)$value
    integrand2  <- function(r) ((c^4-r^2*c^2)/(c^2+r^2)^2) * (2*pi)^-0.5*exp(-r^2/2)
    K2 <- integrate(integrand2, lower = -Inf, upper = Inf)$value
  } else if (phi.function==2){
    #Welsch
    c <- 2.985
    fi  <- function(r) (r*exp(-(r/c)^2))
    bfi <- function(r) (exp(-(r/c)^2)*(1-2*(r^2/c^2)))
    #K interal hesab?
    integrand  <- function(r) (r*exp(-(r/c)^2))^2 * (2*pi)^-0.5*exp(-r^2/2)
    K <- integrate(integrand, lower = -Inf, upper = Inf)$value
    integrand2  <- function(r) (exp(-(r/c)^2)*(1-2*(r^2/c^2))) * (2*pi)^-0.5*exp(-r^2/2)
    K2 <- integrate(integrand2, lower = -Inf, upper = Inf)$value
  } else if (phi.function==3) {
    #insha
    c <- 4.685
    fi  <- function(r) (r*(1+(r/c)^4)^-2)
    bfi <- function(r) ((1+(r/c)^4)^-2+r*(-2)*(1+(r/c)^4)^-3*((4*r^3)/c^4))
    #K interal
    integrand  <- function(r) (r*(1+(r/c)^4)^-2)^2 * (2*pi)^-0.5 * exp(-r^2/2)
    K <- integrate(integrand, lower = -Inf, upper = Inf)$value
    integrand2  <- function(r) ((1+(r/c)^4)^-2+r*(-2)*(1+(r/c)^4)^-3*((4*r^3)/c^4)) * (2*pi)^-0.5 * exp(-r^2/2)
    K2 <- integrate(integrand2, lower = -Inf, upper = Inf)$value
  } else if(phi.function==4) {
    #Logistic
    c <- 1.205
    fi  <- function(r) (c*tanh(r/c))
    bfi <- function(r) (sech(r/c)*sech(r/c))
    #K interal
    integrand  <- function(r) (c*tanh(r/c))^2 * (2*pi)^-0.5*exp(-r^2/2)
    K <- integrate(integrand, lower = -Inf, upper = Inf)$value
    integrand2  <- function(r) (sech(r/c)*sech(r/c)) * (2*pi)^-0.5*exp(-r^2/2)
    K2 <- integrate(integrand2, lower = -Inf, upper = Inf)$value
  }

  hata <- 0
  calc_conv <- 1
  i <- 0
  while (calc_conv > converge.v) {
    i <- i+1

    #previous values
    b <- b.prev
    s2 <- as.numeric(sl.prev[1])
    l <- as.numeric(sl.prev[2])

    #estimated values
    bt <- b
    s2t <- s2
    lt <- l

    O <- (I-l*W)
    OO <- O%*%t(O)
    omega <- solve(OO)
    omg <- solve(O)

    Ot <- (I-lt*W)        #by estimated values
    OOt <- Ot%*%t(Ot)     #by estimated values
    omegat <- solve(OOt)  #by estimated values
    omgt <- solve(Ot)     #by estimated values

    #Big Sigma
    SG <- s2*omega
    S <- s2^0.5*omg
    SGi <- s2^-1*OO
    Si <- s2^-0.5*O

    #Phi
    u <- y-X%*%b
    r <- as.vector(Si%*%u)
    rm <- matrix(c(r), nrow = n, ncol = 1)
    fi.r <- apply(rm, 1, fi)
    b.fi.r <- I*apply(rm, 1, bfi)

    #
    A <- 2*l*W%*%t(W)-W-t(W)
    B <- -omega%*%A%*%omega
    C <- -s2*omega%*%A%*%omega
    D <- 2*s2*omega%*%(A%*%omega%*%A-W%*%t(W))%*%omega

    #score functions
    s.b <- (s2t^0.5/s2)*t(X)%*%OO%*%omgt%*%fi.r
    s.s2 <- -(n*K)/(2*s2)+(s2t/(2*s2^2))*t(fi.r)%*%omgt%*%OO%*%omgt%*%fi.r
    s.l <- -K*sum(diag(omg%*%W))+s2t/s2*t(fi.r)%*%omgt%*%O%*%W%*%omgt%*%fi.r

    ssl <- matrix(0,2,1)
    ssl[1,1] <- as.numeric(s.s2)
    ssl[2,1] <- as.numeric(s.l)

    #Fisher Information Matrix
    #Beta
    FIbb <- K2*s2^-1*t(X)%*%OO%*%X

    b.new <- b.prev+solve(FIbb)%*%s.b
    calc_conv.b <- max(abs(b.new-b.prev))
    b.prev <- b.new


    #New Residuals
    bt <- b.new
    u <- y-X%*%bt
    r <- as.vector(Si%*%u)
    rm <- matrix(c(r), nrow = n, ncol = 1)
    fi.r <- apply(rm, 1, fi)
    b.fi.r <- I*apply(rm, 1, bfi)

    s.s2 <- -(n*K)/(2*s2)+(s2t/(2*s2^2))*t(fi.r)%*%omgt%*%OO%*%omgt%*%fi.r
    s.l <- -K*sum(diag(omg%*%W))+s2t/s2*t(fi.r)%*%omgt%*%O%*%W%*%omgt%*%fi.r

    ssl <- matrix(0,2,1)
    ssl[1,1] <- as.numeric(s.s2)
    ssl[2,1] <- as.numeric(s.l)


    #Sigma Lambda
    FI <- matrix(0,2,2)
    FI[1,1] <- (n*K)/(2*s2^2)
    FI[1,2] <- (K/s2)*sum(diag(omg%*%W))
    FI[2,2] <- 2*K*sum(diag(W%*%omg%*%W%*%omg))
    FI[2,1] <- (K/s2)*sum(diag(omg%*%W))

    sl.new <- sl.prev+solve(FI)%*%ssl
    calc_conv.sl <- max(abs(sl.new-sl.prev))
    sl.prev <- sl.new


    #Calculated Converge
    calc_conv <- max(calc_conv.b, calc_conv.sl)


    #print screen
    if (print.values) {
      prdf <- data.frame("Iteration"=i, "Converge"=calc_conv, "S2"=sl.new[1], "Lambda"= sl.new[2], "Beta"=t(as.vector(b.new)))
      names(prdf) <- c("Iteration","Converge","S2","Lambda","intercept", names(x))
      print(prdf)
    }


    #errors

    if (sl.prev[1]<=0) {
      calc_conv <- 0
      hata <- 1
      print('error: negative s2')
    }

    if (abs(sl.new[2])>1) {
      sl.new[2] <- 0
      calc_conv <- 0
      hata <- 1
      print('lambda is out of (-1,1)')
    }

    if (i==iter){
      calc_conv <- 0
      hata <- 1
      print("not converged, increase iteration limit")
    }

  }

  betas <- as.vector(b.new)
  names(betas) <- c("intercept", names(x))
  sigma2 <- sl.new[1]
  lambda <- sl.new[2]

  result <- list(coefficients=betas, lambda=lambda, s2=sigma2)
  return(result)

}
