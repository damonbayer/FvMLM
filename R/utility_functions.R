# Insert Zeros -----------------------------------------------------------------
.InsertZeros <- function(x, redundant.reduced, index = FALSE) {
  # index example
  # redundant.reduced = 1, 4, 5
  # x = 1, 2, 3
  # return 2, 3, 6
  p <- dim(x)[2]
  if (is.null(p)) {
    return(x)
  }

  if (index) {
    new <-  (1:p)[-redundant.reduced][x]

  } else {
    new <- rep(0, length = length(x) + length(redundant.reduced))
    new[-redundant.reduced] <- x
  }
  new
}

# NormalizeG -------------------------------------------------------------------
.NormalizeG <- function(mu, G, mu.G, k, p){

  ifelse(G  == 1,
         mu.new <- rep(mu.G, k),
         mu.new <- cbind(mu[,1:(G - 1)], rep(mu.G, k)))

  # Insert a column for mu.G in the Gth column of mu
  # mu <- cbind(mu[,1:(G - 1)], rep(mu.G, k))
  if (G < p) {
    mu.new <- cbind(unname(mu.new), mu[,G:(p - 1)])
  }

  # Scale everything except the Gth column
  mu.new[,-G] <- mu.new[,-G] * sqrt(1 - mu.new[,G]^2) / sqrt(rowSums(mu.new[,-G]^2))

  return(mu.new)
}

# Likelihood -------------------------------------------------------------------
.CalculateLikelihood <- function(kappa, mu, alpha, x){
  densities <- dmovMF(x = x, theta = kappa * mu, alpha = alpha, log = TRUE)
  sum(densities)
}

# BIC --------------------------------------------------------------------------
.BICVonMises <- function(kappa, mu, alpha, x, method = "original") {
  n <- dim(x)[1]
  p <- dim(x)[2]
  k <- length(alpha)

  param.count <- switch(method,
                        original  = k - 1 + k * p + k,
                        redundant = k - 1 + k * (p - 1) + k,
                        noise     = k - 1 + k * (p - 1) + k + 1)
  likelihood <- .CalculateLikelihood(kappa, mu, alpha, x)
  BIC <- log(n) * param.count - 2 * likelihood

  list(param.count = param.count,
       likelihood = likelihood,
       BIC = BIC)
}

# Q Function -------------------------------------------------------------------
.QFunction <- function(kappa, mu, x, pi){
  -sum(pi * dmovMF(x = x, theta = kappa * mu, log = TRUE))
}


# QfunctionNoiseMuLambda -------------------------------------------------------
.QfunctionNoiseMuLambda <- function(pars, kappa, G, pi, x, n, p, k){
  lambda <- pars[1:k]
  mu <- matrix(pars[(k + 1):(k*p)], nrow = k, ncol = p - 1, byrow = TRUE)
  mu.G <- pars[k * p + 1]

  mu <- .NormalizeG(mu, G, mu.G, k, p)

  res <- 0
  for (i in 1:n) {
    for (h in 1:k) {
      res <- res + pi[i,h] * kappa[h] * (mu.G * x[i,G] + sum(mu[h,-G] * x[i,-G])) +
        lambda[h] * (1 - mu.G^2 - sum(mu[h,-G]^2))
    }
  }
  - res
}

# Expectation Step  ------------------------------------------------------------
.EStep <- function(kappa, mu, alpha, x){
  k <- length(alpha)
  n <- dim(x)[1]


  p <- matrix(nrow = n, ncol = k)

  for (h in 1:k) {
    p[,h] <- alpha[h] * dmovMF(x = x, theta = kappa[h] * mu[h,], alpha = 1)
  }
  p <- p / rowSums(p)
  # work around dmovMfreturning Inf by setting to 1 and 0
  p[is.na(p)] <- 1

  return(p)
}

# Maximization Step ------------------------------------------------------------
.MStep <- function(x, pi, k){
  n <- dim(x)[1]
  p <- dim(x)[2]

  kappa <- NULL
  # mu <- list()
  mu <- matrix(nrow = k, ncol = p)
  alpha <- NULL
  lambda <-  NULL

  r <- NULL
  r.norm <- NULL
  r.bar <- NULL

  for (h in 1:k) {
    alpha[h] <- sum(pi[,h]) / n

    r <- colSums(x * pi[,h])
    r.norm <- sqrt(sum(r^2))
    r.bar <- r.norm / n
    mu[h,] <- r / r.norm

    kappa[h] <- optimize(f = .QFunction,
                         interval = c(0,1500),
                         # interval = c(0,100),
                         mu = mu[h,],
                         x = x,
                         pi = pi[,h])$minimum
    lambda[h] <- kappa[h] / 2 * r.norm
  }
  list(alpha = alpha,
       mu = mu,
       kappa = kappa,
       lambda = lambda,
       pi = pi)
}

# .MStepRedundant --------------------------------------------------------------
.MStepRedundant <- function(x, pi, k, G){
  n <- dim(x)[1]
  p <- dim(x)[2]

  kappa <- vector(mode = "numeric", length = k)
  mu <- matrix(nrow = k, ncol = p)
  alpha <- vector(mode = "numeric", length = k)

  for (h in 1:k) {
    alpha[h] <- sum(pi[,h]) / n

    r <- colSums(x * pi[,h])
    mu.G <- sum(r[G]) / (sqrt(sum(r[G])^2/length(G) + sum(r[-G]^2)) * length(G))
    mu[h,] <- r / (sqrt(sum(r[G])^2/length(G) + sum(r[-G]^2)))
    mu[h,G] <- mu.G

    # kappa[h] <- optimize(f = .QFunction, interval = c(0,100), mu = mu[h,], x = x, pi = pi[,h])$minimum
    kappa[h] <- optimize(f = .QFunction, interval = c(0,1500), mu = mu[h,], x = x, pi = pi[,h])$minimum
  }
  list(alpha = alpha,
       mu = mu,
       kappa = kappa,
       pi = pi)
}

# .MStepNoise ------------------------------------------------------------------
.MStepNoise <- function(x, pi, G, mu, kappa, lambda, k){
  n <- dim(x)[1]
  p <- dim(x)[2]
  k <- dim(pi)[2]

  # alpha
  alpha <- colSums(pi) / n

  # Initial estimates for optimization
  mu.G <- mean(mu[,G])
  mu <- mu[,-G]

  # Optimize to find mu and lambda
  mu.lambda.res <- optim(par = c(lambda, mu, mu.G),
                         f = .QfunctionNoiseMuLambda,
                         kappa = kappa, G = G, pi = pi, x = x,
                         n = n, p = p, k = k,
                         lower = c(rep(0,k), rep(-1,(k * p - 1))),
                         upper = c(rep(Inf,k), rep(1,(k * p - 1))),
                         method = "L-BFGS-B")

  # Unpack optimization results
  param <- mu.lambda.res$par
  lambda <- param[1:k]
  mu <- matrix(param[(k + 1):(k*p)], nrow = k, ncol = p - 1, byrow = TRUE)
  mu.G <- param[k * p + 1]
  mu <- .NormalizeG(mu = mu, G = G, mu.G = mu.G, p = p, k = k)

  # Optimize to find kappa
  for (h in 1:k) {
    kappa.res <- optimize(f = .QFunction,
                          mu = mu[h,],
                          pi = pi[,h], x = x ,interval = c(0, 1500))
    # pi = pi[,h], x = x ,interval = c(0, 100))

    kappa[h] <- kappa.res$minimum
  }
  list(alpha = alpha,
       mu = mu,
       kappa = kappa,
       lambda = lambda,
       pi = pi)
}

# Initialization ---------------------------------------------------------------
.Initialize <- function(x, k, runs = 10, max.iter = 3, tol = 1e-2) {
  n <- dim(x)[1]
  p <- dim(x)[2]

  for (i in 1:runs) {
    ids <- 0
    while (min(table(ids)) < 3) {
      # Randomly choose means
      mu <- x[sample(nrow(x), k),]

      # For each point, compute the distance from the mean
      distances <- skmeans::skmeans_xdist(x, mu)

      # Assign each point to the mean it is closest to
      ids <- apply(distances, 1, which.min)
    }

    # .Initialize new parameters
    best.likelihood <- 0
    alpha <- NULL
    kappa <- NULL
    lambda <- NULL

    # for every cluster
    for (h in 1:k) {
      mu[h,] <- colSums(x[ids == h,]) / dim(x[ids == h,])[1]
      r.norm <- sqrt(sum((colSums(x[ids == h,]))^2))
      r.bar <-  r.norm / sum(ids == h)
      alpha[h] <- sum(ids == h) / length(ids)
      kappa[h] <- (r.bar * p - r.bar^3) / (1 - r.bar^2)
      lambda[h] <- kappa[h] / 2 * r.norm
    }
    parameters <- list(alpha = alpha, mu = mu, kappa = kappa, lambda = lambda)

    # Run EM a few times
    new.likelihood <- .CalculateLikelihood(kappa = parameters$kappa, mu = parameters$mu, alpha = parameters$alpha, x = x)
    old.likelihood <- 9e-10
    i <- 0

    while (abs((new.likelihood - old.likelihood) / old.likelihood) > tol & i < max.iter) {
      pi <- .EStep(kappa = parameters$kappa, mu = parameters$mu, alpha = parameters$alpha, x = x)
      parameters <- .MStep(x = x, pi = pi, k = k)
      i <- i + 1
    }
    current.likelihood <- .CalculateLikelihood(kappa = parameters$kappa, mu = parameters$mu, alpha = parameters$alpha, x = x)

    if (current.likelihood > best.likelihood) {
      best.parameters <- parameters
      best.likelihood <- current.likelihood
    }
  }
  best.parameters
}

# EM ---------------------------------------------------------------------------
.EM <- function(x, k, max.iter = 200, tol = 1e-6, G = NULL, method = "original", initial.params = NULL, use.movMF = F) {
  if (use.movMF & method == 'original') {
    parameters <- list()
    movMF.results <- movMF::movMF(x = x, k = k)
    parameters$alpha <- movMF.results$alpha
    parameters$kappa <- unname(sqrt(rowSums(movMF.results$theta^2)))
    parameters$mu <- movMF.results$theta / parameters$kappa
    parameters$pi <- movMF.results$P
    ids <- unname(apply(parameters$pi, 1, which.max))
    parameters$lambda <- numeric(k)
    for (h in 1:k) {
      r.norm <- sqrt(sum((colSums(x[ids == h,]))^2))
      parameters$lambda[h] <- parameters$kappa[h] / 2 * r.norm
    }
  }
  else {
    if (is.null(initial.params)) {
      parameters <- .Initialize(x = x, k = k)
    } else {
      parameters <- initial.params
    }

    new.likelihood <- .CalculateLikelihood(kappa = parameters$kappa, mu = parameters$mu, alpha = parameters$alpha, x = x)
    old.likelihood <- 9e-10
    i <- 0

    while (abs((new.likelihood - old.likelihood) / old.likelihood) > tol & i < max.iter) {
      pi <- .EStep(kappa = parameters$kappa, mu = parameters$mu, alpha = parameters$alpha, x = x)

      flag <- FALSE

      tryCatch(parameters <- switch(method,
                                    original = .MStep(x = x, pi = pi, k = k),
                                    redundant = .MStepRedundant(x = x, pi = pi, k = k, G = G),
                                    noise = .MStepNoise(x = x, pi = pi, k = k, G = G,
                                                        mu = parameters$mu,
                                                        kappa = parameters$kappa,
                                                        lambda = parameters$lambda)),
               error = function(e) flag <<- TRUE)
      if (flag) {
        # print('I Broke!')
        # print(head(x))
        break}
      old.likelihood <- new.likelihood
      new.likelihood <- .CalculateLikelihood(kappa = parameters$kappa, mu = parameters$mu, alpha = parameters$alpha, x = x)

      i <- i + 1
    }
    # print(paste('Finished', G))
  } # end else
  BIC.results <- .BICVonMises(kappa = parameters$kappa, mu = parameters$mu, alpha = parameters$alpha, x = x, method = method)

  list(parameters = parameters,
       # alpha = parameters$alpha,
       # mu = parameters$mu,
       # kappa = parameters$kappa,
       BIC = BIC.results$BIC,
       likelihood = BIC.results$likelihood,
       param.count = BIC.results$param.count)
}

# Remove Redundant -------------------------------------------------------------
.RemoveRedundant <- function(redundant){
  to.remove <- vector(mode = "numeric")
  freq <- tabulate(redundant)

  for (i in 1:dim(redundant)[1]) {
    ifelse(test = freq[redundant[i,1]] < freq[redundant[i,2]],
           yes = to.remove <- c(to.remove, redundant[i,1]),
           no = to.remove <- c(to.remove, redundant[i,2]))
  }
  unique(to.remove)
}

# Redundant Cutoff Consideration -----------------------------------------------
.RedundantConsideration <- function(x, k, max.iter = 100, tol = 1e-6, cutoff = 0.237, max.consid = 0, use.movMF = F) {
  # print('Entered .RedundantConsideration')
  p <- dim(x)[2]
  G <- t(combn(p, 2))

  # If no max considered, we return all
  if (max.consid == 0) {
    max.consid <- choose(p,2)
  }

  initial.results <- .EM(x = x, k = k, max.iter = max.iter, tol = tol, use.movMF = use.movMF)
  # print('Finished Initial Results in .RedundantConsideration')
  # print('Finished Initial Results in .RedundantConsideration')

  # Vectorized distance
  # dist <- apply(X = G,
  #               MARGIN = 1,
  #               FUN = function(z) {
  #                 print(z)
  #                 as.numeric(dist(rbind(initial.results$parameters$mu[,z[1]], initial.results$parameters$mu[,z[2]])))
  #               }
  # )

  # Loop distance
  dist <- numeric(nrow(G))
  for (i in 1:nrow(G)) {
    dist[i] <- as.numeric(
      dist(rbind(
        initial.results$parameters$mu[,G[i,1]],
        initial.results$parameters$mu[,G[i,2]])))
  }

  # print('Finished dist in .RedundantConsideration')
  # Get the index of the closest distances
  dist.ord <- order(dist)[1:max.consid]
  # print('Finished ord in .RedundantConsideration')
  # Cut any distances beyond cutoff
  final.G.index <- dist.ord[which(dist[dist.ord] < cutoff)]
  # print('Finished cut in .RedundantConsideration')
  matrix(G[final.G.index,], ncol = 2)
}
# Noise Cutoff Consideration ---------------------------------------------------
.NoiseConsideration <- function(x, noise.max.consid) {
  head(order(apply(x, 2, var), decreasing = TRUE), noise.max.consid)
}
