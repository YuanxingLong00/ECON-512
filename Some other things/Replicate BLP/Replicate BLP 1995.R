# open libraries and define the multiplot function
library(hdm)
library(stargazer)
library(ggplot2)
library(nloptr)
library(AER) 

#define a plot function for problem 1 
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# Problem 1: Load data and explore the data
data(BLP)
tab2 <- c(3.393, 6.711, 8.728, 13.074, 68.597)
stopifnot(all(abs(quantile(BLP$price, c(0,0.25,0.5,0.75,1)) - tab2 + 11.761)<0.005))
BLP$price <- BLP$price + 11.761
BLP$trend <- BLP$trend+ 1971
#Plots
M1 <- aggregate(BLP$model.id,by=list(BLP$trend),FUN=length)
p1 <- ggplot(M1,aes(x=Group.1,y=x)) + geom_point(shape=1) + xlab("year") + ylab("Number of Model")
M2 <- aggregate(BLP$price,by=list(BLP$trend),FUN=mean)
p2 <- ggplot(M2,aes(x=Group.1,y=x)) + geom_point(shape=1) +  geom_smooth(method=lm) + xlab("year") + ylab("Average Price")
M3 <- aggregate(BLP$mpg,by=list(BLP$trend),FUN=mean)
p3 <- ggplot(M3,aes(x=Group.1,y=x)) + geom_point(shape=1) +  geom_smooth(method=lm) + xlab("year") + ylab("MPG")
M4 <- aggregate(BLP$air,by=list(BLP$trend),FUN=mean)
p4 <- ggplot(M4,aes(x=Group.1,y=x)) + geom_point(shape=1) +  geom_smooth(method=lm) + xlab("year") + ylab("Air")
stargazer(BLP,type="latex")
BLP$trend <- BLP$trend - 1970
multiplot(p1,p2,p3,p4,cols=2)

M5 <- aggregate(BLP$space,by=list(BLP$trend),FUN=mean)
M6 <- aggregate(BLP$hpwt,by=list(BLP$trend),FUN=mean)
tab <- cbind(M1,M2$x,M3$x,M4$x,M5$x,M6$x)
colnames(tab) <- c("year","No.model","Price","MPG","Air","Spce","HP/Wt")
stargazer(tab,type="latex",summary=FALSE)


#Problem 2: Logit Model
X <- as.matrix(cbind(1,BLP[,c("hpwt","air","mpg","space","outshr","share")]))
BLP$Z <- hdm:::constructIV(BLP$firm.id,BLP$cdid,BLP$id,X)
ols.reg <- lm(y ~ price + hpwt + air + mpd + space, data = BLP)
iv.reg <- ivreg(y ~ price + hpwt + air + mpd + space |Z + hpwt + air + mpd + space , data = BLP)
ols.reg2 <- lm(log(price) ~  hpwt + air  + mpg + trend + space, data = BLP)
stargazer(ols.reg,iv.reg,ols.reg2,type="latex")
elasticity.iv <- BLP$price  * ols.reg$coefficients[2] 
print(sum(abs(elasticity.iv) < 1))
elasticity.iv <- BLP$price  * iv.reg$coefficients[2] 
print(sum(abs(elasticity.iv) < 1))



#Problem 3: shares and delta

# share.fn(delta,x,log.y,log.yp,v,alpha,sigma )
share.fn <- function(delta,  ## J vector
                     x,      ## J by K
                     log.y,  ## S vector of log(y_i)
                     log.yp, ## S by J of log(y_i - p_j) 
                     v,      ## S by K
                     alpha,  ## scalar
                     sigma)  ## K vector
{
  J <- length(delta)
  S <- length(log.y)
  K <- length(sigma)
  share <- rep(0,J)
  for (i in 1:S){
    share.each <-  exp(x %*% (sigma * v[i,]) + delta + alpha * log.yp[i,])
    share  <-  share + share.each / (sum(share.each) + exp(alpha*log.y[i]))
  }
  share <- share / S
  
  return(share)
}

# delta.fn function 
delta.fn <- function(s,x,log.y,log.yp,v,alpha,sigma,
                     tol=1e-6, tol.s=1e-8,
                     max.iter=100)
{
  delta.new <- log(s) - log(1-sum(s)) ## initial guess
  dd <- 1 ## ||change in delta||
  ds <- 1 ## ||observed shares - share.fn(delta)||
  iter <- 0
  J <- length(delta)
  S <- length(log.y)
  K <- length(sigma)
  exp.delta <- matrix(0,nr=J,nc=S)
  while ((dd>tol) | (ds>tol.s)) {
    delta.old <- delta.new
    
    sm  <- share.fn(delta.old,x,log.y,log.yp,v,alpha,sigma)
    delta.new <- delta.old  + log(s) - log(sm)
    dd <- max(abs(delta.old - delta.new))
    ds <- max(abs(s - sm))
    iter = iter+1
    if (dd < tol | ds < tol.s){
      break
    }
    
    if (iter>max.iter) {
      warning(sprintf("Maximum iterations (%d) reached, returning with norm(delta.new - delta.old) = %.2g", 
                      max.iter, dd))
      break;
    }
  }
  return(delta.new)
}


# test share.fn and delta.fn
#Test code
##construct instruments
##demand instruments
Z <- hdm:::constructIV(BLP$firm.id, BLP$cdid, BLP$id,
                       cbind(1,BLP[,c("hpwt","air","mpd","mpg","space","price")]))
## supply instruments
W <- log(BLP[,c("hpwt","mpg","space","mpd")])
colnames(W) <- paste("log",colnames(W),sep=".")
Wiv <- hdm:::constructIV(BLP$firm.id, BLP$cdid, BLP$id, W)

## Draws for simluting integral
S <- 100 
T <- length(unique(BLP$cdid))
K <- 5
set.seed(41658)
y.s <- matrix(exp(rnorm(S*T,mean=3.082, sd=0.840)),nrow=T,ncol=S)
v.s <- array(rnorm(S*T*K, mean=0,sd=1), dim=c(T,S,K))

## Estimates from Table IV of BLP -- used for comparison and testing
est.blp <- list(alpha=43.501, sigma=c(3.612, 4.628, 1.818, 1.050,
                                      2.056),
                beta=c(-7.061, 2.883, 1.521, -0.122, 3.460),
                gamma=c(0.952, 0.477, 0.619, -.415, -.049, .019))

## Put data into more convenient structure for estimation
est.data <- list(x=as.matrix(cbind(1,BLP[,c("hpwt","air","mpd","space")])),
                 w=as.matrix(cbind(1,log(BLP$hpwt), BLP$air, log(BLP$mpg),
                                   log(BLP$space), BLP$trend )))
est.data$zd <- as.matrix(cbind(est.data$x, Z))
est.data$zs <- as.matrix(cbind(est.data$w, Wiv))
est.data$log.y <- log(y.s)
est.data$log.yp <- list()

## BLP uses log(y - p) in utility function, but some vehicles have
## p>y for some people. BLP do not say what they did in this cases.
## I will take a first order Taylor expansion of log to the left of
## some small number to get an almost log function that is defined
## everywhere. This is very arbitrary though ....
x0 <- 0.1 ## take linear taylor approx to log around x0 for y-p<x0 to
logx0 <- log(x0)
slope <- 1/x0
## avoid log(-)
my.log <- function(x)  ifelse(x>=x0, log(x), logx0 + slope*(x-x0))
dmy.log <- function(x) ifelse(x>=x0, 1/x, slope)
for (t in 1:T) {
  yp <- outer(drop(y.s[t,]), BLP$price[BLP$cdid==t],
              function(x,y) x-y)
  est.data$log.yp[[t]] <- my.log(yp)
  est.data$dlog.yp[[t]] <- dmy.log(yp)
}

## Testing of delta.fn 
t <- 1
inc <- BLP$cdid==t
delta <- rnorm(n=length(BLP$price[inc]))
s <- share.fn(delta, x=drop(est.data$x[inc,]),
              log.y=drop(est.data$log.y[t,]),
              log.yp=est.data$log.yp[[t]],
              v=drop(v.s[t,,]),
              alpha=est.blp$alpha,
              sigma=est.blp$sigma)
d.check <- delta.fn(s, x=drop(est.data$x[inc,]),
                    log.y=drop(est.data$log.y[t,]),
                    log.yp=est.data$log.yp[[t]],
                    v=drop(v.s[t,,]),
                    alpha=est.blp$alpha,
                    sigma=est.blp$sigma, max.iter=10000)

print(summary((abs(delta-d.check))))

#testing deviation function 

dshare.dp <- function(delta,x,log.y, log.yp, dlog.yp, v,alpha,sigma)
{
  J <- length(delta)
  S <- length(log.y)
  K <- length(sigma)
  sj.pj <- rep(0,J) #The diagonal line has a different method 
  sj.pk <- rep(0,J^2)
  for (i in 1:S){
    
    share.each <- exp(x %*% (sigma * v[i,]) + delta + alpha * log.yp[i,])
    share.each <- share.each / (sum(share.each) + exp(alpha*log.y[i]))
    sj.pj <- sj.pj + share.each * (1 - share.each) * ( - alpha * dlog.yp[i,])
    sj.pk <- sj.pk + as.vector(share.each %*% t((share.each * dlog.yp[i,]))) * alpha
  }
  tmp.1 <- diag(as.vector(sj.pj/S))
  tmp.2 <- matrix(sj.pk/S,nc=J)
  conv <- matrix(1,nc=J,nr=J) - diag(J)
  dshare <- t(tmp.2 * conv + tmp.1)
  return(dshare)
}

## Testing dshare.dp
if (!require(numDeriv)) install.packages("numDeriv")
library(numDeriv)
t <- 1
inc <- BLP$cdid==t
dshare.num <- jacobian(function(p) {
  yp <- outer(drop(y.s[t,]), p,
              function(x,y) x-y)
  share.fn(delta, x=drop(est.data$x[inc,]),
           log.y=drop(est.data$log.y[t,]),
           log.yp=my.log(yp),
           v=drop(v.s[t,,]),
           alpha=est.blp$alpha,
           sigma=est.blp$sigma)
}, x=BLP$price[inc])
dshare <- dshare.dp(delta, x=drop(est.data$x[inc,]),
                    log.y=drop(est.data$log.y[t,]),
                    log.yp=est.data$log.yp[[t]], 
                    dlog.yp=est.data$dlog.yp[[t]] ,
                    v=drop(v.s[t,,]),
                    alpha=est.blp$alpha,
                    sigma=est.blp$sigma)
print(summary(diag(dshare-dshare.num)))
print(summary(as.vector(dshare-dshare.num)))


# Problem 4: Code optimization
# use code to compute objective function
moments <- function(alpha,sigma, s,p,x,log.yp, log.y,dlog.yp,
                    v,w,zd,zs, W,
                    market.id,firm.id, delta.tol=1e-10,
                    max.iter=100, compute.variance=FALSE)
{
  ## Find delta and omega for each market
  delta <- rep(NA, length(s))
  omega <- rep(NA, length(s))
  for (t in unique(market.id)) {
    inc <- market.id==t
    delta[inc] <-  delta.fn(s[inc], x[inc,],drop(log.y[t,]),
                            log.yp[[t]] ,
                            v=drop(v[t,,]),
                            alpha,sigma, tol=delta.tol,
                            tol.s=1e-6, max.iter=max.iter)
    
    dShare <- dshare.dp(delta[inc], x[inc,], drop(log.y[t,]),
                        log.yp[[t]], dlog.yp[[t]] , drop(v[t,,]), alpha,sigma)
    dShare <- dShare* outer(firm.id[inc],firm.id[inc],
                            function(x,y) x==y)
    b <- solve(dShare) %*% s[inc]
    omega[inc] <- log(p[inc]-b)
  }
  
  ## Solve for beta and gamma
  X <- as.matrix(rbind(cbind(x,0*w), cbind(0*x,w)))
  Y <- c(delta, omega)
  Z <- as.matrix(rbind(cbind(zd,0*zs), cbind(0*zd, zs)))
  XZ <- t(X) %*% Z 
  B <- solve(XZ %*% W %*% t(XZ)) %*%
    (XZ  %*% W %*% t(Z) %*% Y)
  beta <- B[1:ncol(x)]
  gamma <- B[(1+ncol(x)):(ncol(x)+ncol(w))]
  
  ## Compute GMM objective
  E <- Y - X %*% B
  g <- drop(E)*Z
  G <- colMeans(g)
  obj <- nrow(g)*t(G) %*% W %*% G
  if (compute.variance ){
    ZE <- apply(Z,2,function(x){x*E})
    Omega <- XZ  %*% W %*% t(ZE) %*% ZE %*% W %*% t(XZ)
    H_b <- solve(XZ %*% W %*% t(XZ))
    var <- H_b %*% Omega %*% H_b
  }else{
    var <- NULL
  }
  
  return(list(obj=obj, beta=beta, gamma=gamma,var = var))
}


# use the R profiler to identify which parts of your code are slowest.
Rprof("blp.prof", line.profiling=TRUE) # start the profiler

Z <- as.matrix(rbind(cbind(est.data$zd,0*est.data$zs), cbind(0*est.data$zd, est.data$zs)))
opt.weight <- solve(t(Z)%*% Z)
q <- moments(est.blp$alpha,est.blp$sigma,s=BLP$share,p=BLP$price,
             x=est.data$x, log.y=est.data$log.y,
             log.yp=est.data$log.yp,
             dlog.yp=est.data$dlog.yp,
             v=v.s,w=est.data$w,zd=est.data$zd,zs=est.data$zs,
             W=diag(1,nrow=(ncol(est.data$zd)+ncol(est.data$zs))),
             # W = opt.weight,
             market.id=BLP$cdid, firm.id=BLP$firm.id)
Rprof(NULL) # stop the profiler

summaryRprof("blp.prof", lines="both") # show the results

# modify codes
alpha = est.blp$alpha
sigma = est.blp$sigma
s = BLP$share
p = BLP$price
x=est.data$x 
log.y=est.data$log.y           
log.yp=est.data$log.yp
dlog.yp=est.data$dlog.yp             
v=v.s
w=est.data$w
zd=est.data$zd
zs=est.data$zs
W=diag(1,nrow=(ncol(est.data$zd)+ncol(est.data$zs)))
market.id=BLP$cdid
firm.id=BLP$firm.id
delta.tol=1e-10
max.iter=100


# retest the efficiency of codes 
Rprof("blp.prof", line.profiling=TRUE) # start the profiler

delta <- rep(NA, length(s))
omega <- rep(NA, length(s))
for (t in unique(market.id)) {
  inc <- market.id==t
  delta[inc] <-  delta.fn(s[inc], x[inc,],drop(log.y[t,]),
                          log.yp[[t]] ,
                          v=drop(v[t,,]),
                          alpha,sigma, tol=delta.tol,
                          tol.s=1e-6, max.iter=max.iter)
  
  dShare <- dshare.dp(delta[inc], x[inc,], drop(log.y[t,]),
                      log.yp[[t]], dlog.yp[[t]] , drop(v[t,,]), alpha,sigma)
  dShare <- dShare * outer(firm.id[inc],firm.id[inc],
                           function(x,y) x==y)
  b <- solve(dShare) %*% s[inc]
  omega[inc] <- log(p[inc]-b)
}
Rprof(NULL) # stop the profiler

summaryRprof("blp.prof", lines="both") # show the results

# Problem 5: estimation 
Z <- as.matrix(rbind(cbind(est.data$zd,0*est.data$zs), cbind(0*est.data$zd, est.data$zs)))
opt.weight <- solve(t(Z)%*% Z)
obj.func <- function(alpha.sigma) {  m<-moments(alpha.sigma[1],alpha.sigma[2:6],s=BLP$share,p=BLP$price,
                                                x=est.data$x, log.y=est.data$log.y,
                                                log.yp=est.data$log.yp,
                                                dlog.yp=est.data$dlog.yp,
                                                v=v.s,w=est.data$w,
                                                zd=est.data$zd,
                                                zs=est.data$zs,
                                                W=diag(1,nrow=(ncol(est.data$zd)+ncol(est.data$zs))),
                                                # W = opt.weight,
                                                market.id=BLP$cdid, firm.id=BLP$firm.id)

m$obj
}

x0 <- c(est.blp$alpha,est.blp$sigma)
obj.func(c(est.blp$alpha,est.blp$sigma))
step.a <- nloptr(x0=x0 ,eval_f = obj.func, lb = x0 * 0.1 , opts=list(algorithm="NLOPT_LN_BOBYQA",tol=1e-2)) #print_level=3 ( use for printing each step)
sol <- step.a$solution
step.b <- moments(sol[1],sol[2:6],s=BLP$share,p=BLP$price,
                  x=est.data$x, log.y=est.data$log.y,
                  log.yp=est.data$log.yp,
                  dlog.yp=est.data$dlog.yp,
                  v=v.s,w=est.data$w,
                  zd=est.data$zd,
                  zs=est.data$zs,
                  W=diag(1,nrow=(ncol(est.data$zd)+ncol(est.data$zs))),
                  # W = opt.weight,
                  market.id=BLP$cdid, firm.id=BLP$firm.id, compute.variance = TRUE)


est.result.1 <- matrix(c( step.a$solution[1], step.a$solution[2:6],  step.b$beta,  step.b$gamma),nc=1)
rownames(est.result.1) <- c("alpha",paste("sigma",c("constant","hpwt","air","mpd","space"),sep = "."),paste("beta",c("constant","hpwt","air","mpd","space"),sep = "."),paste("gamma",c("constant","hpwt","air","mpd","space","trend"),sep = "."))

stargazer(est.result.1,type="latex")
#,control = list("algorithm" = "NLOPT_LD_SBPLX" , "xtol_rel"=1.0e-2,
# "#print_level" = 2))

alpha = est.blp$alpha
sigma = est.blp$sigma
s = BLP$share
p = BLP$price
x=est.data$x
w=est.data$w
log.y=est.data$log.y           
log.yp=est.data$log.yp
dlog.yp=est.data$dlog.yp             
v=v.s
# w=est.data$w
zd=est.data$zd
zs=est.data$zs
inc <- BLP$cdid == 20 #Use 1990 data to produce table VI

brand.list <- c("MZ323 ","NISENT","FDESCO","CVCAVA","HDACCO",
                "BKCENT","NIMAXI","ACLEGE","LNTOWN","CDSEVI", "LXLS40","BW735i")

matrix.inc <- c()
els.data <- subset(BLP,cdid==20)

for (each in brand.list){
  matrix.inc <- append(matrix.inc, which(els.data$model.name == each) )
}

t = 20
els.delta <-  delta.fn(s[inc], x[inc,],drop(log.y[t,]),
                       log.yp[[t]] ,
                       v=drop(v.s[t,,]),
                       alpha,sigma, tol=1e-2,
                       tol.s=1e-6, max.iter=1000)

els.dShare <- dshare.dp(els.delta, x[inc,], drop(log.y[t,]),
                        log.yp[[t]], dlog.yp[[t]] , drop(v[t,,]), alpha,sigma)

els.dShare <- els.dShare[matrix.inc,matrix.inc]
els.dShare <- apply(els.dShare,2, function(x){ x  * els.data$price[matrix.inc]})
els.dShare <- apply(els.dShare,1, function(x){ x  / els.data$share[matrix.inc]})

els.table <- els.dShare
colnames(els.table) <- brand.list
rownames(els.table) <- brand.list
stargazer(els.table,type="text")


# Problem 6: Inference

delta.dr <- function(delta,  ## J vector
                     x,      ## J by K
                     log.y,  ## S vector of log(y_i)
                     log.yp, ## S by J of log(y_i - p_j) 
                     v,      ## S by K
                     alpha,  ## scalar
                     sigma)  ## K vector
{
  J <- length(delta)
  S <- length(log.y)
  K <- length(sigma)
  # exp.delta <- matrix(0,nr=J,nc=S)
  delta.dr <- matrix(0,nr=J,nc=K+1)
  for (i in 1:S){
    each.log.yp <- log.yp[i,]
    share.each <-  exp(x %*% (sigma * v[i,]) + delta + alpha * log.yp[i,]) #J by 1 
    share.each <- share.each / (sum(share.each) + exp(alpha*log.y[i])) #J by 1
    delta.dr[,1]  <- delta.dr[,1] + share.each * ( each.log.yp - sum(share.each * each.log.yp))
    svx <-  (share.each ) %*% t(v[i,]) * x #J by K
    delta.dr[,2:(K+1)] <- svx - share.each  %*% t(apply(svx,2,sum) )
  }
  delta.dr <- delta.dr / S
  
  # stop("TODO: compute J vector of shares in this market")
  return(delta.dr)
}


alpha.sigma.variance <- function(alpha,sigma, s,p,x,log.yp, log.y,dlog.yp,
                                 v,w,zd,zs, W,
                                 market.id,firm.id, var)
{
  ## Find delta and omega for each market
  delta.derivative <- matrix(0,nr=length(s),nc=dim(x)[2]+1)
  omega.derivative <- matrix(0,nr=length(s),nc=dim(x)[2]+1)
  delta <- rep(NA, length(s))
  omega <- rep(NA, length(s))
  for (t in unique(market.id)) {
    inc <- market.id==t
    delta.derivative[inc,] <-  delta.dr(s[inc], x[inc,],drop(log.y[t,]),
                                        log.yp[[t]] ,
                                        v=drop(v[t,,]),
                                        alpha,sigma)
    delta.each <-  delta.fn(s[inc], x[inc,],drop(log.y[t,]),
                            log.yp[[t]] ,
                            v=drop(v[t,,]),
                            alpha,sigma, tol=1e-2,
                            tol.s=1e-6, max.iter=1000)
    delta[inc] <- delta.each
    dShare <- dshare.dp(delta.each, x[inc,], drop(log.y[t,]),
                        log.yp[[t]], dlog.yp[[t]] , drop(v[t,,]), alpha,sigma)
    dShare <- dShare* outer(firm.id[inc],firm.id[inc],
                            function(x,y) x==y)
    trans <- solve(dShare)
    
    b <-  trans %*% s[inc]
    omega[inc] <- log(p[inc]-b)
    for (k in 1:(dim(x)[2]+1)){
      each.dr <- diag(delta.derivative[inc,k])
      # omega.derivative[inc,k] <- diag(as.vector(1 / (p[inc]-b)))  %*% trans %*% each.dr %*% dShare %*% t(trans) %*% s[inc]
      omega.derivative[inc,k] <- diag(as.vector(1 / (p[inc]-b)))  %*% each.dr %*% dShare %*% s[inc]
    }
    
  }
  B.dr <- rbind(delta.derivative,omega.derivative)
  X <- as.matrix(rbind(cbind(x,0*w), cbind(0*x,w)))
  Y <- c(delta, omega)
  Z <- as.matrix(rbind(cbind(zd,0*zs), cbind(0*zd, zs)))
  XZ <- t(X) %*% Z 
  M.B <- solve(XZ %*% W %*% t(XZ)) %*%
    (XZ  %*% W %*% t(Z))
  alpha.sigma.var <- t(M.B %*% B.dr) %*% var %*% (M.B %*% B.dr)
  
}


var <- step.b$var # get the variance co-variance matrix.

alpha.sigma.var <- alpha.sigma.variance(sol[1],sol[2:6],s=BLP$share,p=BLP$price,
                                        x=est.data$x, log.y=est.data$log.y,
                                        log.yp=est.data$log.yp,
                                        dlog.yp=est.data$dlog.yp,
                                        v=v.s,w=est.data$w,
                                        zd=est.data$zd,
                                        zs=est.data$zs,
                                        W=diag(1,nrow=(ncol(est.data$zd)+ncol(est.data$zs))),
                                        # W = opt.weight,
                                        market.id=BLP$cdid, firm.id=BLP$firm.id,var=var)

sd <- matrix(c(sqrt(diag(alpha.sigma.var)),sqrt(diag(var))),nc=1) #Calculate standard error
rownames(sd) <- c("alpha",paste("sigma",c("constant","hpwt","air","mpd","space"),sep = "."),paste("beta",c("constant","hpwt","air","mpd","space"),sep = "."),paste("gamma",c("constant","hpwt","air","mpd","space","trend"),sep = "."))


est.result.1 <- matrix(c( step.a$solution[1], step.a$solution[2:6],  step.b$beta,  step.b$gamma),nc=1)
rownames(est.result.1) <- c("alpha",paste("sigma",c("constant","hpwt","air","mpd","space"),sep = "."),paste("beta",c("constant","hpwt","air","mpd","space"),sep = "."),paste("gamma",c("constant","hpwt","air","mpd","space","trend"),sep = "."))

result <- cbind(est.result.1,sd)
colnames(result) <- c("estimate","sd")
stargazer(result,type="latex")





