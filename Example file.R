## Load the Tensor-mixtures.R file ##

library(abind)

#### Simulate data from EI-EV-EVV TVSEN-Ms ####

set.seed(123)

p <- 2
q <- 3
r <- 4
k <- 2
num <- 200
pis <- rep(1/k,k)

M1 <- array(data = c(1,2,1.5,1.5,1.5,2,
                     3,1,3,0.5,3.5,1,
                     1,0,2,1,3,0,
                     2,1,2.5,1,2.5,2),
            dim = c(p,q,r))
M2 <- M1+1

Sigma1 <- Sigma2 <- list()
Sigma1[[1]] <- matrix(data = c(0.71,0.00,
                               0.00,1.41), p , p)
Sigma1[[2]] <- matrix(data = c(0.81,-0.40,-0.77,
                               -0.40,1.49, -0.53,
                               -0.77,-0.53,2.33), q , q)
Sigma1[[3]] <- matrix(data = c(0.55,-0.11,0.14,0.16,
                               -0.11,0.41,0.02,-0.10,
                               0.14,0.02,0.60,0.16,
                               0.16,-0.10,0.16,0.45), r , r)

Sigma2[[1]] <- matrix(data = c(0.71,0.00,
                               0.00,1.41), p , p)
Sigma2[[2]] <- matrix(data = c(0.82,0.52,-0.50,
                               0.52,1.56,0.77,
                               -0.50,0.77,2.24), q , q)
Sigma2[[3]] <- matrix(data = c(0.63,-0.08,-0.05,0.02,
                               -0.08,0.45,-0.25,-0.04,
                               -0.05,-0.25,0.45,-0.01,
                               0.02,-0.04,-0.01,0.54), r , r)

theta1 <- 0.1
theta2 <- 0.05

temp <- t(rmultinom(num, size = 1, prob = pis))
num1 <- colSums(temp)[1]
num2 <- colSums(temp)[2]

X1 <- TrSenTin(n = num1, M = M1, Sigma = Sigma1, theta = theta1, density = "TVSEN")
X2 <- TrSenTin(n = num2, M = M2, Sigma = Sigma2, theta = theta2, density = "TVSEN")

data <- abind(X1, X2, along = 4)

#### Fit a single model until convergence ####

mod <- list(c("EI", "EV", "EVV"))

init.TVSEN <- TensMixtPars_init(X = data, k = k, nstartR = 50, nThreads = 1, mod.list = mod, density = "TVSEN")
fit.TVSEN <- TensMixtPars_fit(X = data, init.par = init.TVSEN, nThreads = 1)
best.TVSEN <- extract.bestM(results = fit.TVSEN, top = 1)

#### Fit, for varying k, all the models via the fitting strategy ####

k <- 1:3
mod <- "all"
nThreads <- 14 # Modify this argument according to the number of cores
               # you want (and can) use for parallelization
n.top <- as.integer(7^2 * 14 * length(k) * 0.01) # top 1% of models

## fit up to 10 iterations ##

init.TVSEN.pt1 <- TensMixtPars_init(X = data, k = k, nstartR = 50, nThreads = nThreads, mod.list = mod, density = "TVSEN")
fit.TVSEN2.pt1 <- TensMixtPars_fit(X = data, init.par = init.TVSEN.pt1, nThreads = nThreads, maxit = 10)
best.TVSEN2.pt1 <- extract.bestM(results = fit.TVSEN2.pt1, top = n.top)
ext.pt1 <- extract.shortEM(Extract.bestM = best.TVSEN2.pt1)
ext.ini.pt1 <- filt.init(res.init = init.TVSEN.pt1, Extract.shortEM = ext.pt1)

## fit the best models until convergence ## 

res.fin <- win.fin <- list()

for (j in 1:length(ext.ini.pt1)) {
  res.fin[[j]] <- TensMixtPars_fit(X = data, init.par = ext.ini.pt1[[j]], nThreads = nThreads)
  win.fin[[j]] <- extract.bestM(results = res.fin[[j]], top = 1)
}




