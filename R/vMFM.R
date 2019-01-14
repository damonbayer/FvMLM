#' Variable Selection for Von Mises–Fisher Mixtures
#'
#' @param x data.frame or matrix to perform clustering and variable selection
#' @param k Integer. The number of clusters.
#' @param max.iter Integer. The maximum number of iterations for each EM run.
#' @param tol Numeric. The tolerance value used to determine if convergence has been reached.
#' @param red.cutoff (optional) Numeric. the cutoff value to use for selecting redundant columns by distance.
#' @param red.max.consid (optional) integer. The number of redundant variable pairs to consider.
#' @param red.method "dist" (default), to select red.max.consid by minimum distance, "all" to select all variables, or "rand" to randomly select red.max.consid pairs
#' @param noise.max.consid (optional) integer. The number of noise variables to consider.
#' @param do.redundant logical. True (default) if redundant variables should be removed.
#' @param do.noise logical. True (default) if noise variables should be removed.
#' @param use.movMF logical. True if the movMF package should be used when possible. False (Default) to use built-in EM algorithm.
#' @return A list object with the redundant variable pairs (if any), noise variables (if any), and final clustering.
vMFM <- function(x, k, max.iter = 100, tol = 1e-6,
                 red.cutoff = Inf, red.max.consid = 0, red.method = 'dist',
                 noise.max.consid = 0,
                 do.redundant = T, do.noise = T, use.movMF = F) {
  # .Initialize everything NULL incase some operations are not performed
  G <- initial.results <- BIC.redundant <- param.count.redundant <-
    likelihood.redundant <- redundant <- redundant.reduced <-
    after.redundant.results <- noise.to.check <- BIC.noise <-
    param.count.noise <- likelihood.noise <- noise <- final.results <-
    timer <- NULL

  timer <- Sys.time()
  initial.results <- .EM(x = x, k = k, max.iter = max.iter, tol = tol, use.movMF = use.movMF)
  # print('Intial Results Finished')

  p <- dim(x)[2]
  G <- NULL

  if (do.redundant){
    # TmpSave
    if (!is.null(tmp.save.path)) {
      # print('SAVING!')
      save(list = ls(envir = environment(), all.names = TRUE),
           file = tmp.save.path,
           envir = environment())
    }
    # print('red.method')
    # print(red.method)
    G <- switch(red.method,
                all  = t(combn(p, 2)),
                dist = .RedundantConsideration(x = x, k = k,
                                               max.iter = max.iter,
                                               tol = tol,
                                               cutoff = red.cutoff,
                                               max.consid = red.max.consid, use.movMF = use.movMF),
                rand = t(combn(p, 2))[sample(1:choose(p, 2), red.max.consid, replace = F),])
    # print(G[1:20,])
    # print('Finished G')

    # TmpSave
    if (!is.null(tmp.save.path)) {
      # print('SAVING!')
      save(list = ls(envir = environment(), all.names = TRUE),
           file = tmp.save.path,
           envir = environment())
    }

    # Beginning Redundant
    BIC.redundant <- vector(mode = "numeric", dim(G)[2])
    param.count.redundant <- vector(mode = "numeric", dim(G)[2])
    likelihood.redundant <- vector(mode = "numeric", dim(G)[2])

    redundant <- NULL
    redundant.reduced <- NULL

    # print('Starting Redundant Loop')

    for (i in 1:dim(G)[1]) {
      cols <- G[i,]
      res <- .EM(x = x, k = k, method = "redundant", G = cols, max.iter = max.iter, tol = tol, initial.params = initial.results$parameters)
      BIC.redundant[i] <- res$BIC
      param.count.redundant[i] <- res$param.count
      likelihood.redundant[i] <- res$likelihood

      if (BIC.redundant[i] < initial.results$BIC) {
        redundant <- rbind(redundant, G[i,])
      }

      # TmpSave
      if (!is.null(tmp.save.path)) {
        # print('SAVING!')
        save(list = ls(envir = environment(), all.names = TRUE),
             file = tmp.save.path,
             envir = environment())
      }
    }

    if (!is.null(redundant)) {
      redundant.reduced <- .RemoveRedundant(redundant)
      # print(c('redundamt.reduced:', redundant.reduced))
      x <- x[, -redundant.reduced, drop = FALSE]
    }

    # Rescale x to be on unit hypershphere
    x <- x / sqrt(rowSums(x^2))
  }
  p <- dim(x)[2]
  if (p <= 1 | !do.noise) {
    return(list(
      last.complete = 'redundant',
      G = G,
      initial.results = initial.results,
      BIC.redundant = BIC.redundant,
      param.count.redundant = param.count.redundant,
      likelihood.redundant = likelihood.redundant,
      redundant = redundant,
      redundant.reduced = redundant.reduced,
      x = x,
      timer = Sys.time() - timer))
  }
  after.redundant.results <- .EM(x = x, k = k, max.iter = max.iter, tol = tol, use.movMF = use.movMF)
  # Beginning noise
  param.count.noise <- vector(mode = "numeric", p)
  likelihood.noise <- vector(mode = "numeric", p)
  BIC.noise <- vector(mode = "numeric", p)

  noise <- NULL
  # print('Starting Noise!')
  if (noise.max.consid == 0) {
    noise.to.check <- 1:p
  } else {
    noise.to.check <- .NoiseConsideration(x = x, noise.max.consid = noise.max.consid)
  }
  # print('Starting Noise!')
  for (i in 1:length(noise.to.check)) {
    # print(paste('Noise column', i, 'of', length(noise.to.check)))
    ntc <- noise.to.check[i]
    res <- .EM(x = x, k = k, method = "noise", G = ntc, tol = tol, max.iter = max.iter, initial.params = after.redundant.results$parameters)
    BIC.noise[ntc] <- res$BIC
    param.count.noise[ntc] <- res$param.count
    likelihood.noise[ntc] <- res$likelihood

    if (BIC.noise[ntc] < after.redundant.results$BIC) {
      noise <- c(noise, ntc)
    }

    # TmpSave
    if (!is.null(tmp.save.path)) {
      # print('SAVING!')
      save(list = ls(envir = environment(), all.names = TRUE),
           file = tmp.save.path,
           envir = environment())
    }
  }

  BIC.noise <- .InsertZeros(x = BIC.noise, redundant.reduced = redundant.reduced)
  param.count.noise <- .InsertZeros(x = param.count.noise,
                                    redundant.reduced = redundant.reduced)
  likelihood.noise <- .InsertZeros(x = likelihood.noise,
                                   redundant.reduced = redundant.reduced)
  noise <-  .InsertZeros(x = noise,
                         redundant.reduced = redundant.reduced,
                         index = TRUE)

  if (!is.null(noise)) {
    x <- x[, -noise, drop = FALSE]
  }

  p <- dim(x)[2]
  # print(paste('p', p))

  if (p <= 1) {
    return(list(
      last.complete = 'noise',
      G = G,
      initial.results = initial.results,
      BIC.redundant = BIC.redundant,
      param.count.redundant = param.count.redundant,
      likelihood.redundant = likelihood.redundant,
      redundant = redundant,
      redundant.reduced = redundant.reduced,
      after.redundant.results = after.redundant.results,
      noise.to.check = noise.to.check,
      BIC.noise = BIC.noise,
      param.count.noise = param.count.noise,
      likelihood.noise = likelihood.noise,
      noise = noise,
      x = x,
      timer = Sys.time() - timer))
  }

  # Rescale x to be on unit hypershphere
  x <- x / sqrt(rowSums(x^2))

  # print('beginning final EM Run')
  final.results <- .EM(x = x, k = k, max.iter = max.iter, tol = tol, use.movMF = use.movMF)

  list(
    last.complete = 'final',
    G = G,
    initial.results = initial.results,
    BIC.redundant = BIC.redundant,
    param.count.redundant = param.count.redundant,
    likelihood.redundant = likelihood.redundant,
    redundant = redundant,
    redundant.reduced = redundant.reduced,
    after.redundant.results = after.redundant.results,
    noise.to.check = noise.to.check,
    BIC.noise = BIC.noise,
    param.count.noise = param.count.noise,
    likelihood.noise = likelihood.noise,
    noise = noise,
    final.results = final.results,
    x = x,
    timer = Sys.time() - timer)
}
