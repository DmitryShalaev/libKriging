install.packages("../rlibkriging_0.0-0_R_x86_64-pc-linux-gnu.tar.gz",repos=NULL)
library(rlibkriging)
f <- function(X) apply(X, 1, function(x) prod(sin(2*pi*x)))
logn <- 1 #seq(1, 2.5, by=.1)
n <- floor(10^logn)
d <- 2
set.seed(1234)
X <- matrix(runif(n*d),ncol=d)
Y <- f(X)
# DiceView::sectionview3d(f)
# rgl::points3d(cbind(X,Y))
k <- DiceKriging::km(design=X, response=Y, covtype = "gauss")
mll_fun <- function(x) -apply(x,1,
                            function(theta) 
                              DiceKriging::logLikFun(theta,k)
                            )
#DiceView::sectionview3d(mll_fun)
k_C <- rlibkriging::ordinary_kriging(0,NaN,t(Y),X,"gauss")
rlibkriging::ordinary_kriging_model(k_C)
mll_fun_C <- function(x) -apply(x,1,
                                function(theta) 
                                  rlibkriging::ordinary_kriging_logLikelihood(k_C,theta)
                                )
#DiceView::sectionview3d(mll_fun_C,xlim=c(0.0001,1))

x=seq(-1.0001,1,,41)
# contour(x,x,matrix(mll_fun(expand.grid(x,x)),nrow=length(x)),nlevels = 50)
contour(x,x,matrix(mll_fun_C(expand.grid(x,x)),nrow=length(x)),nlevels = 50)

mll_grad <- function(x) -apply(x,1,
                             function(theta) {
                               e=new.env()
                               DiceKriging::logLikFun(theta,k,envir=e) # to correctly setup e matrices for later grad eval
                               return(DiceKriging::logLikGrad(theta,k,envir=e))}
                              )

mll_grad_C <- function(x) -apply(x,1,
                               function(theta) 
                                 rlibkriging::ordinary_kriging_logLikelihoodGrad(k_C,theta)
                                )

delta=.001
for (x in seq(-1.0001,1,,21))
  for (y in seq(-1.0001,1,,21)) {
    g = mll_grad_C(cbind(x,y))
    arrows(x,y,x-delta*g[1,],y-delta*g[2,],length = .1)
  }



set.seed(12345)
X0 = matrix(runif(2),ncol=2)
#write.csv(X0,"X0.csv")

Xn_optimR = matrix(NA,ncol=2,nrow=0)
Xn_optimC = matrix(NA,ncol=2,nrow=0)

# for (ix0 in 1:nrow(X0)) {
  ix0=1
  x0 = X0[ix0,]
  
  if (abs(mll_fun(matrix(x0,ncol=2)) - mll_fun_C(matrix(x0,ncol=2)))>1E-7) stop("Wrong f eval")
  # f_optim = function(x) mll_fun(matrix(x,ncol=2))
  if (any(abs(mll_grad(matrix(x0,ncol=2)) - mll_grad_C(matrix(x0,ncol=2)))>1E-7)) stop("Wrong g eval")
  # grad_optim = function(x) t(ll_grad(matrix(x0,ncol=2)))
  
  hist_x = NULL
  last_x = NULL
  f = function(x) {
    n.f <<- n.f+1
    x=matrix(x,ncol=2); 
    if(exists('last_x'))
      if(!is.null(last_x)) 
        lines(x=c(last_x[,1],x[,1]),y=c(last_x[,2],x[,2]),lty=2,col=rgb(0,0,1,.2)); 
    last_x <<- x;
    hist_x <<- rbind(hist_x,x) 
    #points(x,col=rgb(1,0,0,.2)); 
    mll_fun_C(x) #mll_fun(x)
  }
  g = function(x) {
    n.g <<- n.g+1
    x=matrix(x,ncol=2);
    gr=mll_grad_C(x) #ll_grad(x); 
    arrows(x[,1],x[,2],x[,1]-delta*gr[1],x[,2]-delta*gr[2],length = .1,col=rgb(0,0,1,.2)); 
    return(gr)
  }
  
  n.f <<- n.g <<- 0
  o=optim(x0,f,g,method = "L-BFGS-B", control=list(maxit=10)) #,lower = c(0,0),upper=c(1,1))
  points(o$par[1],o$par[2],col='blue',pch='x')
  text(x=x0[1],y=x0[2],paste0(n.f,",",n.g),col='blue')
  Xn_optimR = rbind(Xn_optimR,o$par)
  
  XG = bench_OptimLL(k_C,t(x0))
  X = abs(XG$fun)
  G = XG$grad
  #points(X[,1],X[,2],col=rgb(0,0,1,.2),pch=20)
  lines(X[,],col=rgb(1,0,0,.5),lty=2)
  points(X[nrow(X),1],X[nrow(X),2],col='red',pch='x')
  n.f <<- nrow(X); n.g <<- nrow(G)
  text(x=x0[1]+.01,y=x0[2]+.01,paste0(n.f,",",n.g),col='red')
  Xn_optimC  = rbind(Xn_optimC,X[nrow(X),])
# }

#write.csv(Xn_optimR,"Xn_optimR.csv")
#write.csv(Xn_optimC,"Xn_optimC.csv")


