library(foreach)
library(doSNOW)
library(progress)

### Other packages that must be installed before using the code ##
# R.utils
# withr
# mclust 
# rTensor
# tidyr
# data.table
# rlist
# expint
# zipfR

TensMixtPars_init <- function(X, k, mod.list = "all", nstartR = 50, nThreads = 1, density = "TVN") {
  TensMixt_init <- function(X, k, mod.list, nstartR = 50, density) {
    dMVnorm <- function(X, M, U, V) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)

      return(pdf)
    }
    tr <- function(x) {
      return(sum(diag(x)))
    }
    split.along.dim <- function(a, n) {
      stats::setNames(
        lapply(base::split(a, arrayInd(seq_along(a), dim(a))[, n]),
          array,
          dim = dim(a)[-n], dimnames(a)[-n]
        ),
        dimnames(a)[[n]]
      )
    }
    dMVsen <- function(X, M, U, V, theta) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      intf <- function(w, gamm) {
        w^((p * r) / 2) * exp(-w * gamm)
      }
      pdfinteg <- sapply(1:num, function(i) stats::integrate(intf, lower = 1, upper = Inf, gamm = delta[i] / 2 + theta)$value)
      pdfconst <- (2 * pi)^(-(p * r) / 2) * theta * exp(theta) * det(U)^(-r / 2) * det(V)^(-p / 2)
      PDF <- pdfconst * pdfinteg

      return(PDF)
    }
    dMVtin <- function(X, M, U, V, theta) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      intf <- function(w, del) {
        w^((p * r) / 2) * exp((-w / 2) * del)
      }
      pdfinteg <- sapply(1:num, function(i) stats::integrate(intf, lower = (1 - theta), upper = 1, del = delta[i])$value)
      pdfconst <- (2 * pi)^(-(p * r) / 2) * (1 / theta) * det(U)^(-r / 2) * det(V)^(-p / 2)
      PDF <- pdfconst * pdfinteg

      return(PDF)
    }

    # Dimensions

    dm <- dim(X)
    d.cov <- length(head(dm, -1))
    num <- tail(dm, 1) # sample size
    X.list <- split.along.dim(X, length(dm))

    # Create some objects

    prior <- matrix(NA, nstartR, k)
    llk <- rep(NA, nstartR)

    M <- Mresh <- Sct <- Sigma <- SigmaW <- vector(mode = "list", length = nstartR)
    for (i in 1:nstartR) {
      M[[i]] <- Mresh[[i]] <- Sct[[i]] <- Sigma[[i]] <- SigmaW[[i]] <- vector(mode = "list", length = k)

      for (j in 1:k) {
        Mresh[[i]][[j]] <- Sct[[i]][[j]] <- Sigma[[i]][[j]] <- SigmaW[[i]][[j]] <- vector(mode = "list", length = d.cov)

        for (d in 1:d.cov) {
          Sct[[i]][[j]][[d]] <- vector(mode = "list", length = num)
        }
      }
    }

    Xresh <- vector(mode = "list", length = d.cov)
    for (i in 1:d.cov) {
      Xresh[[i]] <- R.utils::wrap(X, map = list(i, NA, length(dm)))
    }

    sel <- matrix(NA, k, prod(head(dm, -1)))
    dens <- array(0, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))

    ## Random initialization ##

    eu <- matrix(0, nrow = num, ncol = k)
    classy <- numeric(num)
    rand.start <- matrix(0, nstartR, k)
    nu.start <- matrix(0, nstartR, k)

    withr::with_seed(head(dm, 1) * (prod(head(dm, -1)) - sum(head(dm, -1))) + head(dm, 1), for (i in 1:nstartR) {
      rand.start[i, ] <- sample(c(1:num), k)
      if (density == "TVSEN") {
        nu.start[i, ] <- runif(k, 0.05, 1)
      }
      if (density == "TVTIN") {
        nu.start[i, ] <- runif(k, 0.6, 0.95)
      }
    })

    for (t in 1:nstartR) {
      skip_to_next <- FALSE

      ### part 0 ###

      tryCatch(
        {
          sec <- rand.start[t, ]
          nu <- nu.start[t, ]

          for (j in 1:k) {
            sel[j, ] <- as.vector(X.list[[sec[j]]])
          }

          for (j in 1:k) {
            for (i in 1:(num)) {
              eu[i, j] <- as.vector(dist(rbind(as.vector(X.list[[i]]), sel[j, ])))
            }
          }

          for (i in 1:(num)) {
            classy[i] <- which.min(eu[i, ])
          }

          z <- mclust::unmap(classy)

          ### part 1 ###

          for (j in 1:k) {
            M[[t]][[j]] <- rowSums(X * z[, j][slice.index(X, length(dm))], dims = d.cov) / sum(z[, j])

            for (i in 1:d.cov) {
              Mresh[[t]][[j]][[i]] <- R.utils::wrap(M[[t]][[j]], map = list(i, NA))

              for (n in 1:num) {
                Sct[[t]][[j]][[i]][[n]] <- z[n, j] * (Xresh[[i]][, , n] - Mresh[[t]][[j]][[i]]) %*% t(Xresh[[i]][, , n] - Mresh[[t]][[j]][[i]])
              }

              SigmaW[[t]][[j]][[i]] <- Reduce("+", Sct[[t]][[j]][[i]])
            }
          }

          for (i in 1:d.cov) {
            if (i != d.cov) {
              if (mod.list[i] == "II") {
                for (j in 1:k) {
                  Sigma[[t]][[j]][[i]] <- diag(1, dm[i], dm[i])
                }
              }

              if (mod.list[i] == "EI") {
                if (k == 1) {
                  Sigma[[t]][[1]][[i]] <- diag(diag(SigmaW[[t]][[1]][[i]]), dm[i], dm[i]) / (det(diag(diag(SigmaW[[t]][[1]][[i]]), dm[i], dm[i])))^(1 / dm[i])
                } else {
                  temp <- array(NA, dim = c(dm[i], dm[i], k))

                  for (j in 1:k) {
                    temp[, , j] <- SigmaW[[t]][[j]][[i]]
                  }

                  deltaVX <- diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i]) / (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])

                  for (j in 1:k) {
                    Sigma[[t]][[j]][[i]] <- deltaVX
                  }
                }
              }

              if (mod.list[i] == "EE") {
                if (k == 1) {
                  Sigma[[t]][[1]][[i]] <- SigmaW[[t]][[1]][[i]] / (det(SigmaW[[t]][[1]][[i]]))^(1 / dm[i])
                } else {
                  temp <- array(NA, dim = c(dm[i], dm[i], k))
                  for (j in 1:k) {
                    temp[, , j] <- SigmaW[[t]][[j]][[i]]
                  }

                  for (j in 1:k) {
                    Sigma[[t]][[j]][[i]] <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))
                  }
                }
              }
            } else {
              if (mod.list[i] == "EII") {
                if (k == 1) {
                  phiX <- tr(SigmaW[[t]][[1]][[i]]) / (prod(dm))
                } else {
                  temp <- array(NA, dim = c(dm[i], dm[i], k))
                  for (j in 1:k) {
                    temp[, , j] <- SigmaW[[t]][[j]][[i]]
                  }

                  phiX <- tr(rowSums(temp, dims = 2)) / (prod(dm))
                }

                for (j in 1:k) {
                  Sigma[[t]][[j]][[i]] <- phiX * diag(1, dm[i], dm[i])
                }
              }

              if (mod.list[i] == "EEI") {
                if (k == 1) {
                  deltaUX <- diag(diag(SigmaW[[t]][[1]][[i]]), dm[i], dm[i]) / (det(diag(diag(SigmaW[[t]][[1]][[i]]), dm[i], dm[i])))^(1 / (dm[i]))

                  phiX <- dm[i] * (det(diag(diag(SigmaW[[t]][[1]][[i]]), dm[i], dm[i])))^(1 / dm[i]) / (prod(dm))
                } else {
                  temp <- array(NA, dim = c(dm[i], dm[i], k))
                  for (j in 1:k) {
                    temp[, , j] <- SigmaW[[t]][[j]][[i]]
                  }

                  deltaUX <- diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i]) / (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])

                  phiX <- dm[i] * (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i]) / (prod(dm))
                }

                for (j in 1:k) {
                  Sigma[[t]][[j]][[i]] <- phiX * deltaUX
                }
              }

              if (mod.list[i] == "EEE") {
                if (k == 1) {
                  Sigma[[t]][[1]][[i]] <- dm[i] * SigmaW[[t]][[1]][[i]] / (prod(dm))
                } else {
                  temp <- array(NA, dim = c(dm[i], dm[i], k))
                  for (j in 1:k) {
                    temp[, , j] <- SigmaW[[t]][[j]][[i]]
                  }

                  for (j in 1:k) {
                    Sigma[[t]][[j]][[i]] <- dm[i] * rowSums(temp, dims = 2) / (prod(dm))
                  }
                }
              }
            }
          }

          for (j in 1:k) {
            if (density == "TVN") {
              dens[, j] <- dMVnorm(X = Xresh[[1]], M = Mresh[[t]][[j]][[1]], U = Sigma[[t]][[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[t]][[j]]), -1)))
            } else if (density == "TVSEN") {
              dens[, j] <- dMVsen(X = Xresh[[1]], M = Mresh[[t]][[j]][[1]], U = Sigma[[t]][[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[t]][[j]]), -1)), theta = nu[j])
            } else if (density == "TVTIN") {
              dens[, j] <- dMVtin(X = Xresh[[1]], M = Mresh[[t]][[j]][[1]], U = Sigma[[t]][[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[t]][[j]]), -1)), theta = nu[j])
            }
          }

          if (k == 1) {
            prior[t, ] <- 1
          } else {
            prior[t, ] <- colMeans(z)
          }

          ### part 2 ###

          # mixture density

          numerator <- matrix(rep(prior[t, ], num), num, k, byrow = TRUE) * dens
          mixt.dens <- rowSums(numerator)
          llk[t] <- sum(log(mixt.dens))

          if (any(prior[t, ] <= 0.05)) {
            llk[t] <- NA
          }
        },
        error = function(e) {
          skip_to_next <<- TRUE
        }
      )

      if (skip_to_next) {
        next
      }
    }

    df <- data.frame(llk = llk, pos = c(1:nstartR))
    df <- tidyr::drop_na(df)
    df <- df[!is.infinite(rowSums(df)), ]
    bestR <- head(data.table::setorderv(df, cols = "llk", order = -1), n = 1)$pos

    return(list(
      model = mod.list, prior = prior[bestR, ], M = M[[bestR]], Sigma = Sigma[[bestR]], nu = nu.start[bestR, ]
    ))
  }
  comb <- function(x, ...) {
    lapply(
      seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }

  dm <- dim(X)
  d.cov <- length(head(dm, -1))
  mod7 <- mod14 <- mod.tot <- list()

  if (!is.list(mod.list)) {
    for (i in 1:(d.cov - 1)) {
      mod7[[i]] <- as.character(c("II", "EI", "VI", "EE", "VE", "EV", "VV"))
    }

    req.model7 <- expand.grid(mod7, stringsAsFactors = F)
    mod.last <- as.character(c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV"))

    for (j in 1:14) {
      mod14[[j]] <- rep(mod.last[j], nrow(req.model7))
      mod.tot[[j]] <- cbind(req.model7, mod14[[j]])
    }

    req.model <- rlist::list.rbind(mod.tot)
  } else {
    req.model <- data.frame(matrix(unlist(mod.list), length(mod.list), d.cov, byrow = T), stringsAsFactors = FALSE)
  }

  names(req.model) <- paste(rep("V", each = ncol(req.model)), 1:ncol(req.model), sep = "")

  nest.EII <- c("EII", "VII")
  nest.EEI <- c("EEI", "VEI", "EVI", "VVI")
  nest.EEE <- c("EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")
  nest.II <- ("II")
  nest.EI <- c("EI", "VI")
  nest.EE <- c("EE", "VE", "EV", "VV")

  list.comb <- data.frame(matrix(NA, nrow = nrow(req.model), ncol = ncol(req.model)))
  for (i in 1:nrow(req.model)) {
    for (j in 1:(ncol(req.model) - 1)) {
      if (req.model[i, j] %in% nest.II) {
        list.comb[i, j] <- "II"
      }
      if (req.model[i, j] %in% nest.EI) {
        list.comb[i, j] <- "EI"
      }
      if (req.model[i, j] %in% nest.EE) {
        list.comb[i, j] <- "EE"
      }
    }

    if (req.model[i, ncol(req.model)] %in% nest.EII) {
      list.comb[i, ncol(req.model)] <- "EII"
    }
    if (req.model[i, ncol(req.model)] %in% nest.EEI) {
      list.comb[i, ncol(req.model)] <- "EEI"
    }
    if (req.model[i, ncol(req.model)] %in% nest.EEE) {
      list.comb[i, ncol(req.model)] <- "EEE"
    }
  }
  list.comb2 <- unique(list.comb)

  oper <- vector(mode = "list", length = length(k))

  for (g in 1:length(k)) {
    print(paste(paste("Initializing Parsimonious", density), paste("mixtures with k =", k[g])))

    cluster <- makeCluster(nThreads, type = "SOCK")
    registerDoSNOW(cluster)

    pb <- progress::progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = nrow(list.comb2),
      complete = "=", # Completion bar character
      incomplete = "-", # Incomplete bar character
      current = ">", # Current bar character
      width = 100
    )
    progress <- function(n) {
      pb$tick()
    }
    opts <- list(progress = progress)

    oper[[g]] <- foreach(l = 1:nrow(list.comb2), .combine = "comb", .multicombine = TRUE, .init = list(list()), .options.snow = opts) %dopar% {
      res <- TensMixt_init(X = X, k = k[g], mod.list = as.character(list.comb2[l, ]), nstartR = nstartR, density = density)

      list(res)
    }

    if (g == 1) {
      oper2 <- foreach(i = 1:nrow(req.model), .combine = "comb", .multicombine = TRUE, .init = list(list())) %dopar% {
        for (j in 1:nrow(list.comb2)) {
          if (all(list.comb[i, ] == oper[[1]][[1]][[j]][["model"]])) {
            res <- j
          }
        }

        list(res)
      }
    }

    stopCluster(cluster)
    registerDoSEQ()
  }

  return(list(
    results = oper,
    k = k,
    req.model = req.model,
    init.used = list.comb,
    index = unlist(oper2[[1]]),
    density = density
  ))
} 

TensMixtPars_fit <- function(X, init.par = NULL, tol = 0.001, maxit = 500, nThreads = 1) {
  k <- init.par[[2]]
  density <- init.par[[6]]
  list.comb <- init.par[[3]]
  pt.mod <- nrow(list.comb)

  tol2 <- 0.001
  maxit2 <- 100

  comb <- function(x, ...) {
    lapply(
      seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }
  oper <- vector(mode = "list", length = length(k))
  time <- numeric(length(k))

  TensMixt_fit <- function(X, k, init.par = NULL, mod.list = NULL, tol = 0.001, tol2 = 0.001, maxit = 500, maxit2 = 100, density) {
    ptm <- proc.time()

    # Functions

    dMVnorm <- function(X, M, U, V) {
      tr <- function(x) {
        return(sum(diag(x)))
      }

      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)

      return(pdf)
    }
    tr <- function(x) {
      return(sum(diag(x)))
    }
    split.along.dim <- function(a, n) {
      setNames(
        lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
          array,
          dim = dim(a)[-n], dimnames(a)[-n]
        ),
        dimnames(a)[[n]]
      )
    }
    dMVsen <- function(X, M, U, V, theta) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
        num <- 1
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      intf <- function(w, gamm) {
        w^((p * r) / 2) * exp(-w * gamm)
      }
      pdfinteg <- sapply(1:num, function(i) stats::integrate(intf, lower = 1, upper = Inf, gamm = delta[i] / 2 + theta)$value)
      pdfconst <- (2 * pi)^(-(p * r) / 2) * theta * exp(theta) * det(U)^(-r / 2) * det(V)^(-p / 2)
      PDF <- pdfconst * pdfinteg

      return(PDF)
    }
    dMVtin <- function(X, M, U, V, theta) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X

      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
        num <- 1
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }

      intf <- function(w, del) {
        w^((p * r) / 2) * exp((-w / 2) * del)
      }
      pdfinteg <- sapply(1:num, function(i) stats::integrate(intf, lower = (1 - theta), upper = 1, del = delta[i])$value)
      pdfconst <- (2 * pi)^(-(p * r) / 2) * (1 / theta) * det(U)^(-r / 2) * det(V)^(-p / 2)
      PDF <- pdfconst * pdfinteg

      return(PDF)
    }
    Mstep_AECM <- function(X, k, weights = NULL, M, sigmaU, sigmaV) {
      f1 <- function(par, weights, X, M, sigmaU, sigmaV) {
        theta <- par
        pll <- sum(weights * log(dMVtin(X = X, M = M, U = sigmaU, V = sigmaV, theta = theta)))
        return(pll)
      }
      res <- suppressWarnings(stats::optimize(
        f = f1, interval = c(0, 1), weights = weights,
        X = X, M = M, sigmaU = sigmaU, sigmaV = sigmaV, maximum = TRUE
      ))
      theta <- res$maximum

      return(theta)
    }

    # Dimensions

    dm <- dim(X)
    d.cov <- length(head(dm, -1))
    num <- tail(dm, 1) # sample size
    X.list <- split.along.dim(X, length(dm))

    ## Objects

    Mresh <- Sct <- Sigma.inv <- SigmaW <- delta.cv <- gamma <- ftemp <- tempomega <- temp.numdelta <- V_EVV.UY.K <- vector(mode = "list", length = k)
    for (j in 1:k) {
      Mresh[[j]] <- Sct[[j]] <- Sigma.inv[[j]] <- SigmaW[[j]] <-
        delta.cv[[j]] <- gamma[[j]] <- ftemp[[j]] <- tempomega[[j]] <-
        temp.numdelta[[j]] <- V_EVV.UY.K[[j]] <- vector(mode = "list", length = d.cov)

      for (d in 1:d.cov) {
        Sct[[j]][[d]] <- vector(mode = "list", length = num)
      }
    }

    phiX.k <- temp.phi2 <- numeric(k)
    tempWX_EEV <- tempWX_EV <- vector("list", k) # for EEV, VEV, EV

    Xresh <- vector(mode = "list", length = d.cov)
    for (i in 1:d.cov) {
      Xresh[[i]] <- R.utils::wrap(X, map = list(i, NA, length(dm)))
    }

    ## Other objects

    post <- dens <- array(0, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))
    w <- matrix(0, nrow = num, ncol = k)

    # Preliminary definition of convergence criterions for EM/MM algorithms

    check <- 0
    check2 <- 0
    loglik.old <- -Inf
    loglik.new <- NULL
    ll <- NULL
    MM.old <- -Inf
    m.iter <- 0
    m.iter2 <- 0

    ## Algorithm ###

    M <- init.par$M
    Sigma <- init.par$Sigma
    prior <- init.par$prior
    theta <- init.par$nu

    for (j in 1:k) {
      for (i in 1:d.cov) {
        Mresh[[j]][[i]] <- R.utils::wrap(M[[j]], map = list(i, NA))
      }
    }

    ei <- gammaVX <- vector(mode = "list", length = d.cov)
    for (i in 1:d.cov) {
      ei[[i]] <- vector(mode = "list", length = k)
    }

    for (i in 1:d.cov) {
      for (j in 1:k) {
        ei[[i]][[j]] <- eigen(Sigma[[j]][[i]])

        if (i == d.cov) {
          phiX.k[j] <- prod(ei[[i]][[j]][["values"]])^(1 / dm[i])
          delta.cv[[j]][[i]] <- diag(ei[[i]][[j]][["values"]] / phiX.k[j])
          gamma[[j]][[i]] <- ei[[i]][[j]][["vectors"]]
        } else {
          delta.cv[[j]][[i]] <- diag(ei[[i]][[j]][["values"]])
          gamma[[j]][[i]] <- ei[[i]][[j]][["vectors"]]
        }
      }

      temp <- array(NA, dim = c(dm[i], dm[i], k))

      for (j in 1:k) {
        temp[, , j] <- gamma[[j]][[i]]
      }

      gammaVX[[i]] <- rowSums(temp, dims = 2) / k
    }

    ### Estimation ###

    while (check < 1) {
      m.iter <- m.iter + 1

      ### E - STEP ###

      if (density == "TVN") {
        for (j in 1:k) {
          dens[, j] <- dMVnorm(X = Xresh[[1]], M = Mresh[[j]][[1]], U = Sigma[[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[j]]), -1)))
        }
      }

      if (density == "TVSEN") {
        for (j in 1:k) {
          dens[, j] <- dMVsen(X = Xresh[[1]], M = Mresh[[j]][[1]], U = Sigma[[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[j]]), -1)), theta = theta[j])
        }

        for (j in 1:k) {
          delta <- sapply(1:num, function(i) {
            as.vector(t(as.vector(X.list[[i]]) - as.vector(M[[j]])) %*%
              rTensor::kronecker_list(lapply(rev(Sigma[[j]]), solve)) %*% (as.vector(X.list[[i]]) - as.vector(M[[j]])))
          })

          numer <- expint::gammainc(a = (prod(head(dm, -1)) / 2 + 2), x = (delta / 2 + theta[j]))
          den <- (delta / 2 + theta[j]) * expint::gammainc(a = (prod(head(dm, -1)) / 2 + 1), x = (delta / 2 + theta[j]))

          numer[numer < .Machine$double.xmin] <- .Machine$double.xmin
          den[den < .Machine$double.xmin] <- .Machine$double.xmin

          wtt <- numer / den
          wtt[wtt < 1] <- 1.001
          w[, j] <- wtt
        }
      }

      if (density == "TVTIN") {
        for (j in 1:k) {
          dens[, j] <- dMVtin(X = Xresh[[1]], M = Mresh[[j]][[1]], U = Sigma[[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[j]]), -1)), theta = theta[j])
        }

        for (j in 1:k) {
          delta <- sapply(1:num, function(i) {
            as.vector(t(as.vector(X.list[[i]]) - as.vector(M[[j]])) %*%
              rTensor::kronecker_list(lapply(rev(Sigma[[j]]), solve)) %*% (as.vector(X.list[[i]]) - as.vector(M[[j]])))
          })

          numer <- 2 * (zipfR::Igamma(a = (prod(head(dm, -1)) / 2 + 2), x = (1 - theta[j]) * delta / 2, lower = FALSE) - zipfR::Igamma(a = (prod(head(dm, -1)) / 2 + 2), x = delta / 2, lower = FALSE))
          den <- delta * (zipfR::Igamma(a = (prod(head(dm, -1)) / 2 + 1), x = (1 - theta[j]) * delta / 2, lower = FALSE) - zipfR::Igamma(a = (prod(head(dm, -1)) / 2 + 1), x = delta / 2, lower = FALSE))

          numer[numer < .Machine$double.xmin] <- .Machine$double.xmin
          den[den < .Machine$double.xmin] <- .Machine$double.xmin

          wtt <- numer / den
          wtt[wtt > 1] <- 0.999
          w[, j] <- wtt
        }
      }

      # mixture density

      numerator <- matrix(rep(prior, num), num, k, byrow = TRUE) * dens
      mixt.dens <- rowSums(numerator)

      post <- numerator / mixt.dens

      ### M - STEP ###

      if (density == "TVN") {
        for (j in 1:k) {
          M[[j]] <- rowSums(X * post[, j][slice.index(X, length(dm))], dims = d.cov) / sum(post[, j])

          for (i in 1:d.cov) {
            Mresh[[j]][[i]] <- R.utils::wrap(M[[j]], map = list(i, NA))
          }
        }

        for (i in 1:d.cov) {
          if (i != d.cov) {
            if (mod.list[i] == "II") {
              for (j in 1:k) {
                Sigma[[j]][[i]] <- diag(1, dm[i], dm[i])
              }
            }

            if (mod.list[i] == "EI") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))

              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])

                temp[, , j] <- SigmaW[[j]][[i]]
              }

              deltaVX <- diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i]) / (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])

              for (j in 1:k) {
                Sigma[[j]][[i]] <- deltaVX
              }
            }

            if (mod.list[i] == "VI") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])

                Sigma[[j]][[i]] <- diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]) / (det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i])))^(1 / dm[i])
              }
            }

            if (mod.list[i] == "EE") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])

                temp[, , j] <- SigmaW[[j]][[i]]
              }

              for (j in 1:k) {
                Sigma[[j]][[i]] <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))
              }
            }

            if (mod.list[i] == "VE") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
              }

              while (check2 < 1) {
                m.iter2 <- m.iter2 + 1

                for (j in 1:k) {
                  ftemp[[j]][[i]] <- tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]]) %*% SigmaW[[j]][[i]] - max(eigen(SigmaW[[j]][[i]])$values) * tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]])
                }

                temp2 <- array(NA, dim = c(dm[i], dm[i], k))

                for (j in 1:k) {
                  temp2[, , j] <- ftemp[[j]][[i]]
                }

                f.C <- rowSums(temp2, dims = 2)

                MM.new <- tr(f.C %*% gammaVX[[i]])

                if ((abs(MM.new - MM.old)) < tol2 | m.iter2 == maxit2) {
                  check2 <- 1
                  res.svd.C <- svd(f.C)
                  gammaVX[[i]] <- tcrossprod(res.svd.C$v, res.svd.C$u)
                } else {
                  res.svd.C <- svd(f.C)
                  gammaVX[[i]] <- tcrossprod(res.svd.C$v, res.svd.C$u)
                }

                MM.old <- MM.new
              }

              m.iter2 <- 0
              check2 <- 0
              MM.old <- -Inf

              for (j in 1:k) {
                delta.cv[[j]][[i]] <- diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i]) / (det(diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i])))^(1 / dm[i])
              }

              for (j in 1:k) {
                Sigma[[j]][[i]] <- gammaVX[[i]] %*% tcrossprod(delta.cv[[j]][[i]], gammaVX[[i]])
              }
            }

            if (mod.list[i] == "EV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])

                tempWX_EV[[j]] <- eigen(SigmaW[[j]][[i]])

                gamma[[j]][[i]] <- tempWX_EV[[j]][["vectors"]]

                tempomega[[j]][[i]] <- diag(tempWX_EV[[j]][["values"]], dm[i], dm[i])
              }

              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- tempomega[[j]][[i]]
              }

              deltaVX <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- gamma[[j]][[i]] %*% tcrossprod(deltaVX, gamma[[j]][[i]])
              }
            }

            if (mod.list[i] == "VV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])

                Sigma[[j]][[i]] <- SigmaW[[j]][[i]] / det(SigmaW[[j]][[i]])^(1 / dm[i])
              }
            }
          } else {
            if (mod.list[i] == "EII") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp[, , j] <- SigmaW[[j]][[i]]
              }

              phiX <- tr(rowSums(temp, dims = 2)) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * diag(1, dm[i], dm[i])
              }
            }

            if (mod.list[i] == "VII") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                phiX.k[j] <- tr(SigmaW[[j]][[i]]) / (sum(post[, j]) * prod(head(dm, -1)))
                Sigma[[j]][[i]] <- phiX.k[j] * diag(1, dm[i], dm[i])
              }
            }

            if (mod.list[i] == "EEI") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp[, , j] <- SigmaW[[j]][[i]]
              }

              deltaUX <- diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i]) / (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])

              phiX <- (dm[i] * (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * deltaUX
              }
            }

            if (mod.list[i] == "VEI") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp.numdelta[[j]][[i]] <- (1 / phiX.k[j]) * SigmaW[[j]][[i]]
              }
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- temp.numdelta[[j]][[i]]
              }
              deltaUX <- diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i]) / (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])

              for (j in 1:k) {
                phiX.k[j] <- (tr(solve(deltaUX) %*% SigmaW[[j]][[i]])) / (sum(post[, j]) * prod(head(dm, -1)))

                Sigma[[j]][[i]] <- phiX.k[j] * deltaUX
              }
            }

            if (mod.list[i] == "EVI") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                delta.cv[[j]][[i]] <- diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]) / (det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i])))^(1 / dm[i])

                temp.phi2[j] <- det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]))^(1 / dm[i])
              }

              phiX <- dm[i] * sum(temp.phi2) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * delta.cv[[j]][[i]]
              }
            }

            if (mod.list[i] == "VVI") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                delta.cv[[j]][[i]] <- diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]) / (det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i])))^(1 / dm[i])

                phiX.k[j] <- (dm[i] * det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]))^(1 / dm[i])) / (prod(head(dm, -1)) * sum(post[, j]))

                Sigma[[j]][[i]] <- phiX.k[j] * delta.cv[[j]][[i]]
              }
            }

            if (mod.list[i] == "EEE") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp[, , j] <- SigmaW[[j]][[i]]
              }

              for (j in 1:k) {
                Sigma[[j]][[i]] <- dm[i] * rowSums(temp, dims = 2) / (prod(dm))
              }
            }

            if (mod.list[i] == "VEE") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp.numdelta[[j]][[i]] <- (1 / phiX.k[j]) * SigmaW[[j]][[i]]
              }

              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- temp.numdelta[[j]][[i]]
              }
              deltaUX <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))

              for (j in 1:k) {
                phiX.k[j] <- tr(solve(deltaUX) %*% SigmaW[[j]][[i]]) / (sum(post[, j]) * prod(head(dm, -1)))

                Sigma[[j]][[i]] <- phiX.k[j] * deltaUX
              }
            }

            if (mod.list[i] == "EVE") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
              }
              while (check2 < 1) {
                m.iter2 <- m.iter2 + 1

                for (j in 1:k) {
                  ftemp[[j]][[i]] <- tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]]) %*% SigmaW[[j]][[i]] - max(eigen(SigmaW[[j]][[i]])$values) * tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]])
                }

                temp2 <- array(NA, dim = c(dm[i], dm[i], k))

                for (j in 1:k) {
                  temp2[, , j] <- ftemp[[j]][[i]]
                }

                f <- rowSums(temp2, dims = 2)

                MM.new <- tr(f %*% gammaVX[[i]])

                if ((abs(MM.new - MM.old)) < tol2 | m.iter2 == maxit2) {
                  check2 <- 1
                  res.svd <- svd(f)
                  gammaVX[[i]] <- tcrossprod(res.svd$v, res.svd$u)
                } else {
                  res.svd <- svd(f)
                  gammaVX[[i]] <- tcrossprod(res.svd$v, res.svd$u)
                }

                MM.old <- MM.new
              }

              m.iter2 <- 0
              check2 <- 0
              MM.old <- -Inf

              for (j in 1:k) {
                delta.cv[[j]][[i]] <- diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i]) / (det(diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i])))^(1 / dm[i])

                temp.phi2[j] <- tr(gammaVX[[i]] %*% tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]]) %*% SigmaW[[j]][[i]])
              }

              phiX <- sum(temp.phi2) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * gammaVX[[i]] %*% tcrossprod(delta.cv[[j]][[i]], gammaVX[[i]])
              }
            }

            if (mod.list[i] == "VVE") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
              }

              while (check2 < 1) {
                m.iter2 <- m.iter2 + 1

                for (j in 1:k) {
                  ftemp[[j]][[i]] <- tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]]) %*% SigmaW[[j]][[i]] - max(eigen(SigmaW[[j]][[i]])$values) * tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]])
                }

                temp2 <- array(NA, dim = c(dm[i], dm[i], k))

                for (j in 1:k) {
                  temp2[, , j] <- ftemp[[j]][[i]]
                }

                f <- rowSums(temp2, dims = 2)

                MM.new <- tr(f %*% gammaVX[[i]])

                if ((abs(MM.new - MM.old)) < tol2 | m.iter2 == maxit2) {
                  check2 <- 1
                  res.svd <- svd(f)
                  gammaVX[[i]] <- tcrossprod(res.svd$v, res.svd$u)
                } else {
                  res.svd <- svd(f)
                  gammaVX[[i]] <- tcrossprod(res.svd$v, res.svd$u)
                }

                MM.old <- MM.new
              }

              m.iter2 <- 0
              check2 <- 0
              MM.old <- -Inf

              for (j in 1:k) {
                delta.cv[[j]][[i]] <- diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i]) / (det(diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i])))^(1 / dm[i])
                phiX.k[j] <- (dm[i] * (det(diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i]))^(1 / dm[i]))) / (prod(head(dm, -1)) * sum(post[, j]))
                Sigma[[j]][[i]] <- phiX.k[j] * gammaVX[[i]] %*% tcrossprod(delta.cv[[j]][[i]], gammaVX[[i]])
              }
            }

            if (mod.list[i] == "EEV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                tempWX_EEV[[j]][[i]] <- eigen(SigmaW[[j]][[i]])

                gamma[[j]][[i]] <- tempWX_EEV[[j]][[i]][["vectors"]]

                tempomega[[j]][[i]] <- diag(tempWX_EEV[[j]][[i]][["values"]], dm[i], dm[i])
              }

              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- tempomega[[j]][[i]]
              }

              deltaUX <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))

              phiX <- (dm[i] * ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * gamma[[j]][[i]] %*% tcrossprod(deltaUX, gamma[[j]][[i]])
              }
            }

            if (mod.list[i] == "VEV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                tempWX_EEV[[j]][[i]] <- eigen(SigmaW[[j]][[i]])

                gamma[[j]][[i]] <- tempWX_EEV[[j]][[i]] [["vectors"]]

                tempomega[[j]][[i]] <- diag(tempWX_EEV[[j]][[i]][["values"]], dm[i], dm[i])

                temp.numdelta[[j]][[i]] <- (1 / phiX.k[j]) * tempomega[[j]][[i]]
              }

              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- temp.numdelta[[j]][[i]]
              }

              deltaUX <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))

              for (j in 1:k) {
                phiX.k[j] <- tr(tempomega[[j]][[i]] %*% solve(deltaUX)) / (sum(post[, j]) * prod(head(dm, -1)))

                Sigma[[j]][[i]] <- phiX.k[j] * gamma[[j]][[i]] %*% tcrossprod(deltaUX, gamma[[j]][[i]])
              }
            }

            if (mod.list[i] == "EVV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                V_EVV.UY.K[[j]][[i]] <- SigmaW[[j]][[i]] / ((det(SigmaW[[j]][[i]]))^(1 / dm[i]))

                temp.phi2[j] <- det(SigmaW[[j]][[i]])^(1 / dm[i])
              }

              phiX <- (dm[i] * sum(temp.phi2)) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * V_EVV.UY.K[[j]][[i]]
              }
            }

            if (mod.list[i] == "VVV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }
                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                Sigma[[j]][[i]] <- dm[i] * SigmaW[[j]][[i]] / (sum(post[, j]) * prod(head(dm, -1)))
              }
            }
          }
        }

        for (j in 1:k) {
          dens[, j] <- dMVnorm(X = Xresh[[1]], M = Mresh[[j]][[1]], U = Sigma[[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[j]]), -1)))
        }
      } else {
        for (j in 1:k) {
          M[[j]] <- rowSums(X * (w[, j] * post[, j])[slice.index(X, length(dm))], dims = d.cov) / sum((w[, j] * post[, j]))

          for (i in 1:d.cov) {
            Mresh[[j]][[i]] <- R.utils::wrap(M[[j]], map = list(i, NA))
          }
        }

        for (i in 1:d.cov) {
          if (i != d.cov) {
            if (mod.list[i] == "II") {
              for (j in 1:k) {
                Sigma[[j]][[i]] <- diag(1, dm[i], dm[i])
              }
            }

            if (mod.list[i] == "EI") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))

              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp[, , j] <- SigmaW[[j]][[i]]
              }

              deltaVX <- diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i]) / (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])

              for (j in 1:k) {
                Sigma[[j]][[i]] <- deltaVX
              }
            }

            if (mod.list[i] == "VI") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                Sigma[[j]][[i]] <- diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]) / (det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i])))^(1 / dm[i])
              }
            }

            if (mod.list[i] == "EE") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp[, , j] <- SigmaW[[j]][[i]]
              }

              for (j in 1:k) {
                Sigma[[j]][[i]] <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))
              }
            }

            if (mod.list[i] == "VE") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
              }
              while (check2 < 1) {
                m.iter2 <- m.iter2 + 1

                for (j in 1:k) {
                  ftemp[[j]][[i]] <- tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]]) %*% SigmaW[[j]][[i]] - max(eigen(SigmaW[[j]][[i]])$values) * tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]])
                }

                temp2 <- array(NA, dim = c(dm[i], dm[i], k))

                for (j in 1:k) {
                  temp2[, , j] <- ftemp[[j]][[i]]
                }

                f.C <- rowSums(temp2, dims = 2)

                MM.new <- tr(f.C %*% gammaVX[[i]])

                if ((abs(MM.new - MM.old)) < tol2 | m.iter2 == maxit2) {
                  check2 <- 1
                  res.svd.C <- svd(f.C)
                  gammaVX[[i]] <- tcrossprod(res.svd.C$v, res.svd.C$u)
                } else {
                  res.svd.C <- svd(f.C)
                  gammaVX[[i]] <- tcrossprod(res.svd.C$v, res.svd.C$u)
                }

                MM.old <- MM.new
              }

              m.iter2 <- 0
              check2 <- 0
              MM.old <- -Inf

              for (j in 1:k) {
                delta.cv[[j]][[i]] <- diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i]) / (det(diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i])))^(1 / dm[i])
              }

              for (j in 1:k) {
                Sigma[[j]][[i]] <- gammaVX[[i]] %*% tcrossprod(delta.cv[[j]][[i]], gammaVX[[i]])
              }
            }

            if (mod.list[i] == "EV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                tempWX_EV[[j]] <- eigen(SigmaW[[j]][[i]])

                gamma[[j]][[i]] <- tempWX_EV[[j]][["vectors"]]

                tempomega[[j]][[i]] <- diag(tempWX_EV[[j]][["values"]], dm[i], dm[i])
              }

              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- tempomega[[j]][[i]]
              }

              deltaVX <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- gamma[[j]][[i]] %*% tcrossprod(deltaVX, gamma[[j]][[i]])
              }
            }

            if (mod.list[i] == "VV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                Sigma[[j]][[i]] <- SigmaW[[j]][[i]] / det(SigmaW[[j]][[i]])^(1 / dm[i])
              }
            }
          } else {
            if (mod.list[i] == "EII") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp[, , j] <- SigmaW[[j]][[i]]
              }

              phiX <- tr(rowSums(temp, dims = 2)) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * diag(1, dm[i], dm[i])
              }
            }

            if (mod.list[i] == "VII") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                phiX.k[j] <- tr(SigmaW[[j]][[i]]) / (sum(post[, j]) * prod(head(dm, -1)))
                Sigma[[j]][[i]] <- phiX.k[j] * diag(1, dm[i], dm[i])
              }
            }

            if (mod.list[i] == "EEI") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp[, , j] <- SigmaW[[j]][[i]]
              }

              deltaUX <- diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i]) / (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])

              phiX <- (dm[i] * (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * deltaUX
              }
            }

            if (mod.list[i] == "VEI") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp.numdelta[[j]][[i]] <- (1 / phiX.k[j]) * SigmaW[[j]][[i]]
              }
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- temp.numdelta[[j]][[i]]
              }
              deltaUX <- diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i]) / (det(diag(diag(rowSums(temp, dims = 2)), dm[i], dm[i])))^(1 / dm[i])

              for (j in 1:k) {
                phiX.k[j] <- (tr(solve(deltaUX) %*% SigmaW[[j]][[i]])) / (sum(post[, j]) * prod(head(dm, -1)))

                Sigma[[j]][[i]] <- phiX.k[j] * deltaUX
              }
            }

            if (mod.list[i] == "EVI") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                delta.cv[[j]][[i]] <- diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]) / (det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i])))^(1 / dm[i])

                temp.phi2[j] <- det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]))^(1 / dm[i])
              }

              phiX <- dm[i] * sum(temp.phi2) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * delta.cv[[j]][[i]]
              }
            }

            if (mod.list[i] == "VVI") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                delta.cv[[j]][[i]] <- diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]) / (det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i])))^(1 / dm[i])

                phiX.k[j] <- (dm[i] * det(diag(diag(SigmaW[[j]][[i]]), dm[i], dm[i]))^(1 / dm[i])) / (prod(head(dm, -1)) * sum(post[, j]))

                Sigma[[j]][[i]] <- phiX.k[j] * delta.cv[[j]][[i]]
              }
            }

            if (mod.list[i] == "EEE") {
              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp[, , j] <- SigmaW[[j]][[i]]
              }

              for (j in 1:k) {
                Sigma[[j]][[i]] <- dm[i] * rowSums(temp, dims = 2) / (prod(dm))
              }
            }

            if (mod.list[i] == "VEE") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                temp.numdelta[[j]][[i]] <- (1 / phiX.k[j]) * SigmaW[[j]][[i]]
              }

              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- temp.numdelta[[j]][[i]]
              }
              deltaUX <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))

              for (j in 1:k) {
                phiX.k[j] <- tr(solve(deltaUX) %*% SigmaW[[j]][[i]]) / (sum(post[, j]) * prod(head(dm, -1)))

                Sigma[[j]][[i]] <- phiX.k[j] * deltaUX
              }
            }

            if (mod.list[i] == "EVE") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
              }
              while (check2 < 1) {
                m.iter2 <- m.iter2 + 1

                for (j in 1:k) {
                  ftemp[[j]][[i]] <- tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]]) %*% SigmaW[[j]][[i]] - max(eigen(SigmaW[[j]][[i]])$values) * tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]])
                }

                temp2 <- array(NA, dim = c(dm[i], dm[i], k))

                for (j in 1:k) {
                  temp2[, , j] <- ftemp[[j]][[i]]
                }

                f <- rowSums(temp2, dims = 2)

                MM.new <- tr(f %*% gammaVX[[i]])

                if ((abs(MM.new - MM.old)) < tol2 | m.iter2 == maxit2) {
                  check2 <- 1
                  res.svd <- svd(f)
                  gammaVX[[i]] <- tcrossprod(res.svd$v, res.svd$u)
                } else {
                  res.svd <- svd(f)
                  gammaVX[[i]] <- tcrossprod(res.svd$v, res.svd$u)
                }

                MM.old <- MM.new
              }

              m.iter2 <- 0
              check2 <- 0
              MM.old <- -Inf

              for (j in 1:k) {
                delta.cv[[j]][[i]] <- diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i]) / (det(diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i])))^(1 / dm[i])

                temp.phi2[j] <- tr(gammaVX[[i]] %*% tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]]) %*% SigmaW[[j]][[i]])
              }

              phiX <- sum(temp.phi2) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * gammaVX[[i]] %*% tcrossprod(delta.cv[[j]][[i]], gammaVX[[i]])
              }
            }

            if (mod.list[i] == "VVE") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
              }
              while (check2 < 1) {
                m.iter2 <- m.iter2 + 1

                for (j in 1:k) {
                  ftemp[[j]][[i]] <- tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]]) %*% SigmaW[[j]][[i]] - max(eigen(SigmaW[[j]][[i]])$values) * tcrossprod(solve(delta.cv[[j]][[i]]), gammaVX[[i]])
                }

                temp2 <- array(NA, dim = c(dm[i], dm[i], k))

                for (j in 1:k) {
                  temp2[, , j] <- ftemp[[j]][[i]]
                }

                f <- rowSums(temp2, dims = 2)

                MM.new <- tr(f %*% gammaVX[[i]])

                if ((abs(MM.new - MM.old)) < tol2 | m.iter2 == maxit2) {
                  check2 <- 1
                  res.svd <- svd(f)
                  gammaVX[[i]] <- tcrossprod(res.svd$v, res.svd$u)
                } else {
                  res.svd <- svd(f)
                  gammaVX[[i]] <- tcrossprod(res.svd$v, res.svd$u)
                }

                MM.old <- MM.new
              }

              m.iter2 <- 0
              check2 <- 0
              MM.old <- -Inf

              for (j in 1:k) {
                delta.cv[[j]][[i]] <- diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i]) / (det(diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i])))^(1 / dm[i])
                phiX.k[j] <- (dm[i] * (det(diag(diag(crossprod(gammaVX[[i]], SigmaW[[j]][[i]]) %*% gammaVX[[i]]), dm[i], dm[i]))^(1 / dm[i]))) / (prod(head(dm, -1)) * sum(post[, j]))
                Sigma[[j]][[i]] <- phiX.k[j] * gammaVX[[i]] %*% tcrossprod(delta.cv[[j]][[i]], gammaVX[[i]])
              }
            }

            if (mod.list[i] == "EEV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                tempWX_EEV[[j]][[i]] <- eigen(SigmaW[[j]][[i]])

                gamma[[j]][[i]] <- tempWX_EEV[[j]][[i]][["vectors"]]

                tempomega[[j]][[i]] <- diag(tempWX_EEV[[j]][[i]][["values"]], dm[i], dm[i])
              }

              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- tempomega[[j]][[i]]
              }

              deltaUX <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))

              phiX <- (dm[i] * ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * gamma[[j]][[i]] %*% tcrossprod(deltaUX, gamma[[j]][[i]])
              }
            }

            if (mod.list[i] == "VEV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                tempWX_EEV[[j]][[i]] <- eigen(SigmaW[[j]][[i]])

                gamma[[j]][[i]] <- tempWX_EEV[[j]][[i]] [["vectors"]]

                tempomega[[j]][[i]] <- diag(tempWX_EEV[[j]][[i]][["values"]], dm[i], dm[i])

                temp.numdelta[[j]][[i]] <- (1 / phiX.k[j]) * tempomega[[j]][[i]]
              }

              temp <- array(NA, dim = c(dm[i], dm[i], k))
              for (j in 1:k) {
                temp[, , j] <- temp.numdelta[[j]][[i]]
              }

              deltaUX <- rowSums(temp, dims = 2) / ((det(rowSums(temp, dims = 2)))^(1 / dm[i]))

              for (j in 1:k) {
                phiX.k[j] <- tr(tempomega[[j]][[i]] %*% solve(deltaUX)) / (sum(post[, j]) * prod(head(dm, -1)))

                Sigma[[j]][[i]] <- phiX.k[j] * gamma[[j]][[i]] %*% tcrossprod(deltaUX, gamma[[j]][[i]])
              }
            }

            if (mod.list[i] == "EVV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                V_EVV.UY.K[[j]][[i]] <- SigmaW[[j]][[i]] / ((det(SigmaW[[j]][[i]]))^(1 / dm[i]))

                temp.phi2[j] <- det(SigmaW[[j]][[i]])^(1 / dm[i])
              }

              phiX <- (dm[i] * sum(temp.phi2)) / (prod(dm))

              for (j in 1:k) {
                Sigma[[j]][[i]] <- phiX * V_EVV.UY.K[[j]][[i]]
              }
            }

            if (mod.list[i] == "VVV") {
              for (j in 1:k) {
                for (n in 1:num) {
                  Sct[[j]][[i]][[n]] <- w[n, j] * post[n, j] * (Xresh[[i]][, , n] - Mresh[[j]][[i]]) %*% rTensor::kronecker_list(lapply(rev(Sigma[[j]][-i]), solve)) %*% t(Xresh[[i]][, , n] - Mresh[[j]][[i]])
                }

                SigmaW[[j]][[i]] <- Reduce("+", Sct[[j]][[i]])
                Sigma[[j]][[i]] <- dm[i] * SigmaW[[j]][[i]] / (sum(post[, j]) * prod(head(dm, -1)))
              }
            }
          }
        }
      }

      if (density == "TVSEN") {
        for (j in 1:k) {
          theta[j] <- sum(post[, j]) / (sum(post[, j] * (w[, j] - 1)))
        }

        for (j in 1:k) {
          dens[, j] <- dMVsen(X = Xresh[[1]], M = Mresh[[j]][[1]], U = Sigma[[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[j]]), -1)), theta = theta[j])
        }
      }

      if (density == "TVTIN") {
        for (j in 1:k) {
          theta[j] <- Mstep_AECM(X = Xresh[[1]], weights = post[, j], M = Mresh[[j]][[1]], sigmaU = Sigma[[j]][[1]], sigmaV = rTensor::kronecker_list(head(rev(Sigma[[j]]), -1)))
        }

        for (j in 1:k) {
          dens[, j] <- dMVtin(X = Xresh[[1]], M = Mresh[[j]][[1]], U = Sigma[[j]][[1]], V = rTensor::kronecker_list(head(rev(Sigma[[j]]), -1)), theta = theta[j])
        }
      }

      # Mixture weights #

      if (k == 1) {
        prior <- 1
      } else {
        prior <- colMeans(post)
      }

      # mixture density

      numerator <- matrix(rep(prior, num), num, k, byrow = TRUE) * dens
      mixt.dens <- rowSums(numerator)
      loglik.new <- sum(log(mixt.dens))
      ll <- c(ll, loglik.new)

      # stopping rule

      if ((loglik.new - loglik.old) <= tol) {
        check <- 1
      }

      if (m.iter == maxit) {
        check <- 1
      }

      loglik.old <- loglik.new
    }

    #### Output ####

    # Classification #

    if (k == 1) {
      classification <- rep(1, num)
    } else {
      colnames(post) <- c(1:k)
      classification <- as.numeric(colnames(post)[max.col(post, ties.method = "first")])
    }

    # -------------------- #
    # Information criteria #
    # -------------------- #

    # Number of parameters

    mean.par <- prod(head(dm, -1)) * k
    cov.par <- list()
    for (i in 1:d.cov) {
      if (mod.list[i] == "EII") {
        cov.par[[i]] <- 1
      }
      if (mod.list[i] == "VII") {
        cov.par[[i]] <- k
      }
      if (mod.list[i] == "EEI") {
        cov.par[[i]] <- dm[i]
      }
      if (mod.list[i] == "VEI") {
        cov.par[[i]] <- k + (dm[i] - 1)
      }
      if (mod.list[i] == "EVI") {
        cov.par[[i]] <- 1 + k * (dm[i] - 1)
      }
      if (mod.list[i] == "VVI") {
        cov.par[[i]] <- k * dm[i]
      }
      if (mod.list[i] == "EEE") {
        cov.par[[i]] <- dm[i] * (dm[i] + 1) / 2
      }
      if (mod.list[i] == "VEE") {
        cov.par[[i]] <- k - 1 + dm[i] * (dm[i] + 1) / 2
      }
      if (mod.list[i] == "EVE") {
        cov.par[[i]] <- 1 + k * (dm[i] - 1) + dm[i] * (dm[i] - 1) / 2
      }
      if (mod.list[i] == "VVE") {
        cov.par[[i]] <- k * dm[i] + dm[i] * (dm[i] - 1) / 2
      }
      if (mod.list[i] == "EEV") {
        cov.par[[i]] <- dm[i] + k * dm[i] * (dm[i] - 1) / 2
      }
      if (mod.list[i] == "VEV") {
        cov.par[[i]] <- k + (dm[i] - 1) + (k * dm[i] * (dm[i] - 1) / 2)
      }
      if (mod.list[i] == "EVV") {
        cov.par[[i]] <- 1 + k * (dm[i] * ((dm[i] + 1) / 2) - 1)
      }
      if (mod.list[i] == "VVV") {
        cov.par[[i]] <- k * dm[i] * (dm[i] + 1) / 2
      }

      if (mod.list[i] == "II") {
        cov.par[[i]] <- 0
      }
      if (mod.list[i] == "EI") {
        cov.par[[i]] <- dm[i] - 1
      }
      if (mod.list[i] == "VI") {
        cov.par[[i]] <- k * (dm[i] - 1)
      }
      if (mod.list[i] == "EE") {
        cov.par[[i]] <- dm[i] * ((dm[i] + 1) / 2) - 1
      }
      if (mod.list[i] == "VE") {
        cov.par[[i]] <- k * (dm[i] - 1) + dm[i] * (dm[i] - 1) / 2
      }
      if (mod.list[i] == "EV") {
        cov.par[[i]] <- (dm[i] - 1) + k * dm[i] * (dm[i] - 1) / 2
      }
      if (mod.list[i] == "VV") {
        cov.par[[i]] <- k * (dm[i] * ((dm[i] + 1) / 2) - 1)
      }
    }
    cv.par <- sum(unlist(cov.par))

    weights <- k - 1

    npar <- mean.par + cv.par + weights
    if (density == "TVSEN" | density == "TVTIN") {
      npar <- npar + k
    }

    name <- as.character(mod.list)

    # to be minimized

    BIC <- -2 * loglik.new + npar * log(num)

    ptm2 <- proc.time() - ptm
    time <- ptm2[3]
    if (any(round(prior, digits = 2) <= 0.05)) {
      sp <- 1
    } else {
      sp <- 0
    }

    return(list(
      name = name, prior = prior, M = M, Sigma = Sigma, theta = theta,
      loglik = loglik.new, npar = npar, iter = m.iter, time = time, BIC = BIC,
      class = classification, sp = sp
    ))
  }

  for (g in 1:length(k)) {
    ptm <- proc.time()

    print(paste(paste("Fitting Parsimonious", density), paste("mixtures with k =", k[g])))

    cluster <- makeCluster(nThreads, type = "SOCK")
    registerDoSNOW(cluster)

    pb <- progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = pt.mod,
      complete = "=", # Completion bar character
      incomplete = "-", # Incomplete bar character
      current = ">", # Current bar character
      width = 100
    )
    progress <- function(n) {
      pb$tick()
    }
    opts <- list(progress = progress)

    oper[[g]] <- foreach(l = 1:pt.mod, .combine = "comb", .multicombine = TRUE, .init = list(list()), .options.snow = opts) %dopar% {
      res <- tryCatch(TensMixt_fit(
        X = X, k = k[g], init.par = init.par[[1]][[g]][[1]][[init.par[[5]][l]]],
        mod.list = list.comb[l, ], tol = tol, tol2 = tol2, maxit = maxit, maxit2 = maxit2, density = density
      ), error = function(e) {
        NA
      })

      list(res)
    }

    stopCluster(cluster)
    registerDoSEQ()

    ptm2 <- proc.time() - ptm
    time[g] <- ptm2[3]
  }

  return(list(results = oper, c.time = time, models = list.comb))
} 

Trnorm <- function(n, m.vec, mu = array(0, m.vec), Sigma.list = NULL) {
  K <- length(m.vec)
  kcov <- 1
  for (k in K:1) {
    kcov <- kronecker(kcov, Sigma.list[[k]])
  }
  Data <- array(0, c(m.vec, n))
  ncols <- ncol(kcov)
  vecdata <- matrix(rnorm(n * ncols), ncol = ncols) %*% chol(kcov)
  Data <- array(rep(mu, n), c(m.vec, n)) + array(
    t(vecdata),
    c(m.vec, n)
  )
  return(Data)
}

TrSenTin <- function(n, M, Sigma, theta, density) {
  X <- vector(mode = "list", length = n)
  d.cov <- dim(M)
  M0 <- array(0, dim = d.cov)

  if (density == "TVSEN") {
    for (i in 1:n) {
      w <- 1 + stats::rexp(n = 1, theta)

      X[[i]] <- array(M, dim = c(d.cov, 1)) + Trnorm(n = 1, m.vec = d.cov, mu = M0, Sigma.list = Sigma) / sqrt(w)
    }
  }
  if (density == "TVTIN") {
    for (i in 1:n) {
      w <- stats::runif(n = 1, min = 1 - theta, 1)

      X[[i]] <- array(M, dim = c(d.cov, 1)) + Trnorm(n = 1, m.vec = d.cov, mu = M0, Sigma.list = Sigma) / sqrt(w)
    }
  }

  return(array(unlist(X), dim = c(d.cov, n)))
}

extract.bestM <- function(results, top = 1, short = FALSE) {
  k <- length(results[["results"]])
  num.mod <- length(results[["results"]][[1]][[1]])
  list.mod <- results[["models"]]
  list.mod2 <- do.call("rbind", replicate(k, list.mod, simplify = FALSE))
  count.k <- sort(rep(1:k, num.mod))
  count.mod <- rep(1:num.mod, k)
  list.mod3 <- data.frame(list.mod2, count.k, count.mod)

  allBIC <- numeric(k * num.mod)
  cont <- 0

  if (short == FALSE) {
    for (j in 1:k) {
      for (i in 1:num.mod) {
        if (!all(is.na(results[["results"]][[j]][[1]][[i]]))) {
          if (results[["results"]][[j]][[1]][[i]]$sp == 0) {
            cont <- cont + 1
            allBIC[cont] <- -results[["results"]][[j]][[1]][[i]][["BIC"]]
          } else {
            cont <- cont + 1
            allBIC[cont] <- NA
          }
        } else {
          cont <- cont + 1
          allBIC[cont] <- NA
        }
      }
    }
  } else {
    for (j in 1:k) {
      for (i in 1:num.mod) {
        if (!all(is.na(results[["results"]][[j]][[1]][[i]]))) {
          cont <- cont + 1
          allBIC[cont] <- -results[["results"]][[j]][[1]][[i]][["BIC"]]
        } else {
          cont <- cont + 1
          allBIC[cont] <- NA
        }
      }
    }
  }

  topBIC <- which(allBIC >= sort(allBIC, decreasing = T)[top])
  topBIC.order <- order(allBIC[topBIC], decreasing = T)
  tempBIC <- list.mod3[topBIC[topBIC.order], ]
  bestBIC <- vector(mode = "list", length = top)
  for (i in 1:top) {
    bestBIC[[i]] <- results[["results"]][[as.numeric(tempBIC[i, ncol(tempBIC) - 1])]][[1]][[as.numeric(tempBIC[i, ncol(tempBIC)])]]
  }

  return(bestBIC = bestBIC)
} 

extract.shortEM <- function(Extract.bestM) {
  d.cov <- length(Extract.bestM[[1]][["name"]])
  container <- as.data.frame(matrix(NA, nrow = length(Extract.bestM), ncol = d.cov + 1))

  for (j in 1:length(Extract.bestM)) {
    container[j, 1:d.cov] <- Extract.bestM[[j]][["name"]]
    container[j, d.cov + 1] <- length(Extract.bestM[[j]][["prior"]])
  }

  return(res = container)
} 

filt.init <- function(res.init, Extract.shortEM) {
  dmv <- ncol(Extract.shortEM)
  ext <- list()
  colsToUse <- intersect(colnames(res.init[["req.model"]]), colnames(Extract.shortEM[, 1:(dmv - 1)]))

  pt1 <- split(Extract.shortEM, Extract.shortEM[, dmv])
  ext <- list()
  for (j in 1:length(pt1)) {
    ext[[j]] <- which(!is.na(match(do.call("paste", res.init[["req.model"]][, colsToUse]), do.call("paste", pt1[[j]][, 1:(dmv - 1)][, colsToUse]))))
  }

  grp <- sort(unique(Extract.shortEM[, dmv]))

  new.init <- vector("list", length(ext))
  for (i in 1:length(ext)) {
    new.init[[i]] <- vector("list", 6)

    new.init[[i]][[1]] <- list()
    new.init[[i]][[1]][[1]] <- list()
    new.init[[i]][[1]][[1]][[1]] <- list()
    if (length(res.init[[1]]) == 1) {
      new.init[[i]][[1]][[1]][[1]][unique(res.init[["index"]][ext[[i]]])] <- res.init[[1]][[1]][[1]][unique(res.init[["index"]][ext[[i]]])]
    } else {
      new.init[[i]][[1]][[1]][[1]][unique(res.init[["index"]][ext[[i]]])] <- res.init[[1]][[grp[i]]][[1]][unique(res.init[["index"]][ext[[i]]])]
    }

    new.init[[i]][[2]] <- grp[i]
    new.init[[i]][[3]] <- res.init[[3]][ext[[i]], ]
    new.init[[i]][[4]] <- res.init[[4]][ext[[i]], ]
    new.init[[i]][[5]] <- res.init[[5]][ext[[i]]]
    new.init[[i]][[6]] <- res.init[[6]]
  }

  return(new.init)
} 
