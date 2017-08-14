mledist = function (data, distr, start = NULL, fix.arg = NULL, optim.method = "default", 
          lower = -Inf, upper = Inf, custom.optim = NULL, weights = NULL, 
          silent = TRUE, gradient = NULL, ...) 
{
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else distname <- distr
  ddistname <- paste("d", distname, sep = "")
  if (!exists(ddistname, mode = "function")) 
    stop(paste("The ", ddistname, " function must be defined"))
  if (is.null(custom.optim)) 
    optim.method <- match.arg(optim.method, c("default", 
                                              "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", 
                                              "Brent"))
  start.arg <- start
  if (is.vector(start.arg)) 
    start.arg <- as.list(start.arg)
  txt1 <- "data must be a numeric vector of length greater than 1 for non censored data"
  txt2 <- "or a dataframe with two columns named left and right and more than one line for censored data"
  if (!is.null(weights)) {
    if (any(weights < 0)) 
      stop("weights should be a vector of integers greater than 0")
    # if (!is.allint.w(weights)) 
    #   stop("weights should be a vector of (strictly) positive integers")
    if (length(weights) != NROW(data)) 
      stop("weights should be a vector with a length equal to the observation number")
    warning("weights are not taken into account in the default initial values")
  }
  if (is.vector(data)) {
    cens <- FALSE
    if (!(is.numeric(data) & length(data) > 1)) 
      stop(paste(txt1, txt2))
  }
  else {
    cens <- TRUE
    censdata <- data
    if (!(is.vector(censdata$left) & is.vector(censdata$right) & 
          length(censdata[, 1]) > 1)) 
      stop(paste(txt1, txt2))
    pdistname <- paste("p", distname, sep = "")
    if (!exists(pdistname, mode = "function")) 
      stop(paste("The ", pdistname, " function must be defined to apply maximum likelihood to censored data"))
  }
  if (cens) {
    irow.lcens <- is.na(censdata$left)
    lcens <- censdata[irow.lcens, ]$right
    if (any(is.na(lcens))) 
      stop("An observation cannot be both right and left censored, coded with two NA values")
    irow.rcens <- is.na(censdata$right)
    rcens <- censdata[irow.rcens, ]$left
    irow.ncens <- censdata$left == censdata$right & !is.na(censdata$left) & 
      !is.na(censdata$right)
    ncens <- censdata[irow.ncens, ]$left
    irow.icens <- censdata$left != censdata$right & !is.na(censdata$left) & 
      !is.na(censdata$right)
    icens <- censdata[irow.icens, ]
    data <- c(rcens, lcens, ncens, (icens$left + icens$right)/2)
  }
  argddistname <- names(formals(ddistname))
  chfixstt <- checkparam(start.arg = start.arg, fix.arg = fix.arg, 
                         argdistname = argddistname, errtxt = NULL, data10 = head(data, 
                                                                                  10), distname = distname)
  if (!chfixstt$ok) 
    stop(chfixstt$txt)
  if (is.function(chfixstt$start.arg)) 
    vstart <- unlist(chfixstt$start.arg(data))
  else vstart <- unlist(chfixstt$start.arg)
  if (is.function(fix.arg)) {
    fix.arg.fun <- fix.arg
    fix.arg <- fix.arg(data)
  }
  else fix.arg.fun <- NULL
  if (distname == "unif") {
    par <- c(min = min(data), max = max(data))
    res <- list(estimate = par[!names(par) %in% names(fix.arg)], 
                convergence = 0, loglik = NA, hessian = NA, optim.function = NA, 
                fix.arg = fix.arg)
    return(res)
  }
  if (!cens && is.null(weights)) {
    if ("log" %in% argddistname) {
      fnobj <- function(par, fix.arg, obs, ddistnam) {
        -sum(do.call(ddistnam, c(list(obs), as.list(par), 
                                 as.list(fix.arg), log = TRUE)))
      }
    }
    else {
      fnobj <- function(par, fix.arg, obs, ddistnam) {
        -sum(log(do.call(ddistnam, c(list(obs), as.list(par), 
                                     as.list(fix.arg)))))
      }
    }
  }
  else if (cens && is.null(weights)) {
    argpdistname <- names(formals(pdistname))
    if (("log" %in% argddistname) & ("log.p" %in% argpdistname)) 
      fnobjcens <- function(par, fix.arg, rcens, lcens, 
                            icens, ncens, ddistnam, pdistnam) -sum(do.call(ddistnam, 
                                                                           c(list(x = ncens), as.list(par), as.list(fix.arg), 
                                                                             list(log = TRUE)))) - sum(do.call(pdistnam, 
                                                                                                               c(list(q = lcens), as.list(par), as.list(fix.arg), 
                                                                                                                 list(log = TRUE)))) - sum(do.call(pdistnam, 
                                                                                                                                                   c(list(q = rcens), as.list(par), as.list(fix.arg), 
                                                                                                                                                     list(lower.tail = FALSE), list(log = TRUE)))) - 
        sum(log(do.call(pdistnam, c(list(q = icens$right), 
                                    as.list(par), as.list(fix.arg))) - do.call(pdistnam, 
                                                                               c(list(q = icens$left), as.list(par), as.list(fix.arg)))))
    else fnobjcens <- function(par, fix.arg, rcens, lcens, 
                               icens, ncens, ddistnam, pdistnam) -sum(log(do.call(ddistnam, 
                                                                                  c(list(x = ncens), as.list(par), as.list(fix.arg))))) - 
        sum(log(do.call(pdistnam, c(list(q = lcens), as.list(par), 
                                    as.list(fix.arg))))) - sum(log(1 - do.call(pdistnam, 
                                                                               c(list(q = rcens), as.list(par), as.list(fix.arg))))) - 
        sum(log(do.call(pdistnam, c(list(q = icens$right), 
                                    as.list(par), as.list(fix.arg))) - do.call(pdistnam, 
                                                                               c(list(q = icens$left), as.list(par), as.list(fix.arg)))))
  }
  else if (!cens && !is.null(weights)) {
    fnobj <- function(par, fix.arg, obs, ddistnam) {
      -sum(weights * log(do.call(ddistnam, c(list(obs), 
                                             as.list(par), as.list(fix.arg)))))
    }
  }
  else if (cens && !is.null(weights)) {
    fnobjcens <- function(par, fix.arg, rcens, lcens, icens, 
                          ncens, ddistnam, pdistnam) {
      p1 <- log(do.call(ddistnam, c(list(x = ncens), as.list(par), 
                                    as.list(fix.arg))))
      p2 <- log(do.call(pdistnam, c(list(q = lcens), as.list(par), 
                                    as.list(fix.arg))))
      p3 <- log(1 - do.call(pdistnam, c(list(q = rcens), 
                                        as.list(par), as.list(fix.arg))))
      p4 <- log(do.call(pdistnam, c(list(q = icens$right), 
                                    as.list(par), as.list(fix.arg))) - do.call(pdistnam, 
                                                                               c(list(q = icens$left), as.list(par), as.list(fix.arg))))
      -sum(weights[irow.ncens] * p1) - sum(weights[irow.lcens] * 
                                             p2) - sum(weights[irow.rcens] * p3) - sum(weights[irow.icens] * 
                                                                                         p4)
    }
  }
  owarn <- getOption("warn")
  if (is.null(custom.optim)) {
    hasbound <- any(is.finite(lower) | is.finite(upper))
    if (optim.method == "default") {
      meth <- ifelse(length(vstart) > 1, "Nelder-Mead", 
                     "BFGS")
    }
    else meth <- optim.method
    if (meth == "BFGS" && hasbound && is.null(gradient)) {
      meth <- "L-BFGS-B"
      txt1 <- "The BFGS method cannot be used with bounds without provided the gradient."
      txt2 <- "The method is changed to L-BFGS-B."
      warning(paste(txt1, txt2))
    }
    options(warn = ifelse(silent, -1, 0))
    if (hasbound) {
      if (!is.null(gradient)) {
        opt.fun <- "constrOptim"
      }
      else {
        if (meth == "Nelder-Mead") 
          opt.fun <- "constrOptim"
        else if (meth %in% c("L-BFGS-B", "Brent")) 
          opt.fun <- "optim"
        else {
          txt1 <- paste("The method", meth, "cannot be used by constrOptim() nor optim() without gradient and bounds.")
          txt2 <- "Only optimization methods L-BFGS-B, Brent and Nelder-Mead can be used in such case."
          stop(paste(txt1, txt2))
        }
      }
      if (opt.fun == "constrOptim") {
        npar <- length(vstart)
        lower <- as.double(rep_len(lower, npar))
        upper <- as.double(rep_len(upper, npar))
        haslow <- is.finite(lower)
        Mat <- diag(npar)[haslow, ]
        hasupp <- is.finite(upper)
        Mat <- rbind(Mat, -diag(npar)[hasupp, ])
        colnames(Mat) <- names(vstart)
        rownames(Mat) <- paste0("constr", 1:NROW(Mat))
        Bnd <- c(lower[is.finite(lower)], -upper[is.finite(upper)])
        names(Bnd) <- paste0("constr", 1:length(Bnd))
        initconstr <- Mat %*% vstart - Bnd
        if (any(initconstr < 0)) 
          stop("Starting values must be in the feasible region.")
        if (!cens) {
          opttryerror <- try(opt <- constrOptim(theta = vstart, 
                                                f = fnobj, ui = Mat, ci = Bnd, grad = gradient, 
                                                fix.arg = fix.arg, obs = data, ddistnam = ddistname, 
                                                hessian = !is.null(gradient), method = meth, 
                                                ...), silent = TRUE)
        }
        else opttryerror <- try(opt <- constrOptim(theta = vstart, 
                                                   f = fnobjcens, ui = Mat, ci = Bnd, grad = gradient, 
                                                   ddistnam = ddistname, rcens = rcens, lcens = lcens, 
                                                   icens = icens, ncens = ncens, pdistnam = pdistname, 
                                                   fix.arg = fix.arg, obs = data, hessian = !is.null(gradient), 
                                                   method = meth, ...), silent = TRUE)
        if (!inherits(opttryerror, "try-error")) 
          if (length(opt$counts) == 1) 
            opt$counts <- c(opt$counts, NA)
      }
      else {
        if (!cens) 
          opttryerror <- try(opt <- optim(par = vstart, 
                                          fn = fnobj, fix.arg = fix.arg, obs = data, 
                                          gr = gradient, ddistnam = ddistname, hessian = TRUE, 
                                          method = meth, lower = lower, upper = upper, 
                                          ...), silent = TRUE)
        else opttryerror <- try(opt <- optim(par = vstart, 
                                             fn = fnobjcens, fix.arg = fix.arg, gr = gradient, 
                                             rcens = rcens, lcens = lcens, icens = icens, 
                                             ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, 
                                             hessian = TRUE, method = meth, lower = lower, 
                                             upper = upper, ...), silent = TRUE)
      }
    }
    else {
      opt.fun <- "optim"
      if (!cens) 
        opttryerror <- try(opt <- optim(par = vstart, 
                                        fn = fnobj, fix.arg = fix.arg, obs = data, 
                                        gr = gradient, ddistnam = ddistname, hessian = TRUE, 
                                        method = meth, lower = lower, upper = upper, 
                                        ...), silent = TRUE)
      else opttryerror <- try(opt <- optim(par = vstart, 
                                           fn = fnobjcens, fix.arg = fix.arg, gr = gradient, 
                                           rcens = rcens, lcens = lcens, icens = icens, 
                                           ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, 
                                           hessian = TRUE, method = meth, lower = lower, 
                                           upper = upper, ...), silent = TRUE)
    }
    options(warn = owarn)
    if (inherits(opttryerror, "try-error")) {
      warnings("The function optim encountered an error and stopped.")
      if (getOption("show.error.messages")) 
        print(attr(opttryerror, "condition"))
      return(list(estimate = rep(NA, length(vstart)), convergence = 100, 
                  loglik = NA, hessian = NA, optim.function = opt.fun, 
                  fix.arg = fix.arg, optim.method = meth, fix.arg.fun = fix.arg.fun, 
                  counts = c(NA, NA)))
    }
    if (opt$convergence > 0) {
      warnings("The function optim failed to converge, with the error code ", 
               opt$convergence)
    }
    if (is.null(names(opt$par))) 
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, 
                loglik = -opt$value, hessian = opt$hessian, optim.function = opt.fun, 
                fix.arg = fix.arg, optim.method = meth, fix.arg.fun = fix.arg.fun, 
                weights = weights, counts = opt$counts, optim.message = opt$message)
  }
  else {
    options(warn = ifelse(silent, -1, 0))
    if (!cens) 
      opttryerror <- try(opt <- custom.optim(fn = fnobj, 
                                             fix.arg = fix.arg, obs = data, ddistnam = ddistname, 
                                             par = vstart, ...), silent = TRUE)
    else opttryerror <- try(opt <- custom.optim(fn = fnobjcens, 
                                                fix.arg = fix.arg, rcens = rcens, lcens = lcens, 
                                                icens = icens, ncens = ncens, ddistnam = ddistname, 
                                                pdistnam = pdistname, par = vstart, ...), silent = TRUE)
    options(warn = owarn)
    if (inherits(opttryerror, "try-error")) {
      warnings("The customized optimization function encountered an error and stopped.")
      if (getOption("show.error.messages")) 
        print(attr(opttryerror, "condition"))
      return(list(estimate = rep(NA, length(vstart)), convergence = 100, 
                  loglik = NA, hessian = NA, optim.function = custom.optim, 
                  fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                  counts = c(NA, NA)))
    }
    if (opt$convergence > 0) {
      warnings("The customized optimization function failed to converge, with the error code ", 
               opt$convergence)
    }
    argdot <- list(...)
    method.cust <- argdot[argdot == "method"]
    if (length(method.cust) == 0) {
      method.cust <- NULL
    }
    if (is.null(names(opt$par))) 
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, 
                loglik = -opt$value, hessian = opt$hessian, optim.function = custom.optim, 
                fix.arg = fix.arg, method = method.cust, fix.arg.fun = fix.arg.fun, 
                weights = weights, counts = opt$counts, optim.message = opt$message)
  }
  return(res)
}