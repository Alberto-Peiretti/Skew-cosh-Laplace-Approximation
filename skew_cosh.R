# Skew-cosh: theta = center + L z, with L L' = Sigma.
# delta(z) = (log_post(center+Lz) - log_post(center-Lz))/2.
# The correction is cosh(delta) inside the chosen region and 1 outside.
# score is the full log-posterior gradient at the center.
# posterior_remainder_bound(r) >= sup_{||z||<=r} |delta(z) - score' L z|.
# prior_remainder_bound(r) uses the same formula for the prior alone.
# cutoff='paper': n ||theta-center||^2 <= radius_sq; 'mahalanobis': ||z||^2 <= radius_sq.

## 1) Numerical functions ==============================================

log_add <- function(a, b) {
  m <- pmax(a, b)
  ans <- m + log1p(exp(-abs(a-b)))
  ans[is.infinite(m) & m < 0] <- -Inf
  ans
}

make_choice_tree <- function(logweights) {
  n <- length(logweights); base <- 2^ceiling(log2(n))
  mass <- rep(-Inf,2*base)
  mass[base+seq_len(n)-1L] <- logweights
  if (base>1) for (i in seq.int(base-1L,1L)) mass[i] <- log_add(mass[2*i],mass[2*i+1L])
  list(mass=mass,base=base)
}

log_difference <- function(a, b) {
  ans <- a + log(-expm1(b-a))
  ans[is.infinite(b) & b < 0] <- a[is.infinite(b) & b < 0]
  ans
}

log_sum <- function(x) Reduce(log_add,x,init=-Inf)

envelope_log_masses <- function(nodes) vapply(nodes,function(x) x$logweight,numeric(1))

update_shells <- function(model,nodes) {
  model$shell_nodes <- nodes
  model$shell_logweights <- envelope_log_masses(nodes)
  model$shell_choice_tree <- make_choice_tree(model$shell_logweights)
  model$root_log_bound <- model$log_bound_globale <- log_sum(model$shell_logweights)
  model
}

log_cosh <- function(x) {
  a <- abs(x)
  out <- a + log1p(exp(-2*a)) - log(2)
  small <- a < 1
  out[small] <- log1p(2*sinh(x[small]/2)^2)
  out
}

softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

inside_region <- function(model,z) {
  if (identical(model$cutoff$type,"paper"))
    model$n*sum(as.numeric(model$L %*% z)^2) <= model$cutoff$radius_sq else
    sum(z*z) <= model$radius_sq
}

## 2) Global bound and model construction ===============================

cubic_coefficients <- function(eta,cap,bins) {

  logA <- function(t) {
    out <- numeric(length(t)); j <- t > 0
    out[j] <- t[j]+2*log(-expm1(-t[j]))-2*log(t[j])
    out
  }
  lc <- log_cosh(eta)
  amp <- log(abs(tanh(eta/2)))-log(4)
  ap <- rep(0,length(eta)); dp <- log_add(lc,0)
  ql <- qu <- exp(amp-dp)
  for (j in seq_len(bins)) {
    t <- cap*(j/bins)
    an <- logA(t); dn <- log_add(lc,log_cosh(t))
    ql <- pmin(ql,exp(amp+ap-dn))
    qu <- pmax(qu,exp(amp+an-dp))
    ap <- an; dp <- dn
  }
  if (any(!is.finite(c(ql,qu)))) return(NULL)
  pad <- 128*.Machine$double.eps*(1+qu)
  ql <- pmax(0,ql-pad); qu <- qu+pad
  list(c3=-sign(eta)*(ql+qu),remainder=qu-ql)
}

largest_eigenvalue_bound <- function(A) {
  A <- (A+t(A))/2
  eigen(A,symmetric=TRUE,only.values=TRUE)$values[1L]+
    128*.Machine$double.eps*(1+sum(abs(A)))
}

cubic_tensor <- function(B,weights) {
  d <- ncol(B)
  T <- array(0,c(d,d,d))
  for (j in seq_len(d)) T[j,,] <- crossprod(B,B*(weights*B[,j]))
  T
}

harmonic_decomposition <- function(T,scale=sum(abs(T))) {
  d <- dim(T)[1L]
  v <- vapply(seq_len(d),function(j) sum(diag(matrix(T[j,,],d,d))),numeric(1))
  H <- matrix(0,d,d*d)
  I <- diag(d)
  for (j in seq_len(d)) H[j,] <- as.vector(matrix(T[j,,],d,d)-
    (v[j]*I+outer(v,I[j,])+outer(I[j,],v))/(d+2))
  norm <- svd(H,nu=0,nv=0)$d[1L]+128*.Machine$double.eps*(1+scale+sum(abs(H)))
  list(v=v,norm=norm)
}

trace_envelope <- function(B,s,eta,norm_b,r,bins,
                               coefficient_function=cubic_coefficients) {

  d <- ncol(B)
  co <- coefficient_function(eta,r*norm_b,bins)
  if (is.null(co)) stop("Envelope coefficients are not numerically resolvable.")
  T <- cubic_tensor(B,co$c3)
  harmonic <- harmonic_decomposition(T)
  KR <- crossprod(B,B*(co$remainder*norm_b))
  re <- max(0,largest_eigenvalue_bound(KR))
  list(linear=-harmonic$v/(2*(d+2)),cubic_remainder=(harmonic$norm+re)/6,
       score_norm=sqrt(sum(s^2)),radius=r,
       numerical_scale=1+sum(abs(eta)+r*norm_b))
}

normal_interval_mass <- function(lo, hi) {
  upper <- lo >= 0
  big <- small <- numeric(length(lo))
  big[upper] <- pnorm(lo[upper], lower.tail=FALSE, log.p=TRUE)
  small[upper] <- pnorm(hi[upper], lower.tail=FALSE, log.p=TRUE)
  big[!upper] <- pnorm(hi[!upper], log.p=TRUE)
  small[!upper] <- pnorm(lo[!upper], log.p=TRUE)
  lm <- log_difference(big, small)
  if (any(!is.finite(lm))) {
    stop("Interval mass is not numerically resolvable; no region was discarded.")
  }
  list(logmass=lm, anchor=small, upper=upper)
}

make_skew_cosh <- function(X=NULL, y=NULL, center, Sigma, radius_sq, log_post,
                                 score=NULL, quadratic_prior=FALSE,
                                 offset=0,
                                 rotate=TRUE,
                                 partition=c("radial", "boxes"),
                                 search=c("cosh", "exp"),
                                 delta_evaluation=c("analytic", "log_post"),
                                 bound_padding=1e-10,
                                 coefficient_bins=128L,
                                 proposal_variance="auto", proposal_bins=2048L,
                                 prior=NULL, prior_remainder_bound=NULL, log_prior=NULL,
                                 cutoff=c("mahalanobis","paper"),
                                 posterior_remainder_bound=NULL,n_obs=NULL) {
  cutoff <- match.arg(cutoff)
  partition <- match.arg(partition)
  search <- match.arg(search)
  delta_evaluation <- match.arg(delta_evaluation)
  generic_likelihood <- !is.null(posterior_remainder_bound)
  Sigma <- as.matrix(Sigma)
  parameter_names <- names(center)
  center <- as.numeric(center); names(center) <- parameter_names
  d <- length(center)
  prior_terms <- if (is.null(prior) || generic_likelihood) NULL else local({
    d <- length(center)
    expand <- function(x) rep(x,length.out=d)
    loc <- expand(prior$location); scale <- expand(prior$scale); df <- expand(prior$df)
    family <- prior$family; x0 <- center-loc
    log_density <- switch(family,
      normal=function(x) sum(dnorm(x,loc,scale,log=TRUE)),
      student_t=function(x) sum(dt((x-loc)/scale,df,log=TRUE)-log(scale)),
      laplace=function(x) -sum(abs(x-loc)/scale+log(2*scale)))
    gradient <- switch(family,
      normal=function(x) -(x-loc)/scale^2,
      student_t=function(x) -(df+1)*(x-loc)/(df*scale^2+(x-loc)^2),
      laplace=function(x) -sign(x-loc)/scale)
    third_bound <- NULL
    if (family=="normal") remainder <- function(h) 0
    if (family=="laplace") remainder <- function(h) {
      sum(ifelse(x0==0,0,pmax(0,h-abs(x0))/scale))
    }
    if (family=="student_t") {
      c <- df*scale^2; a <- (df+1)/2
      third <- function(x) {
        u <- x/sqrt(c)
        4*a/c^(3/2)*u*(3-u^2)/(1+u^2)^3
      }
      grad_abs <- function(x) 2*a*abs(x)/(c+x^2)
      third_bound <- function(h) {
        lo <- x0-h; hi <- x0+h
        m3 <- pmax(abs(third(lo)),abs(third(hi)))
        for (u in c(-sqrt(3+2*sqrt(2)),-sqrt(3-2*sqrt(2)),
                     sqrt(3-2*sqrt(2)),sqrt(3+2*sqrt(2)))) {
          point <- u*sqrt(c); hit <- lo<=point & hi>=point
          m3[hit] <- pmax(m3[hit],abs(third(point))[hit])
        }
        m3
      }
      remainder <- function(h) {
        lo <- x0-h; hi <- x0+h
        m3 <- third_bound(h)
        mg <- pmax(grad_abs(lo),grad_abs(hi))
        hit <- (lo<=sqrt(c) & hi>=sqrt(c)) | (lo<=-sqrt(c) & hi>=-sqrt(c))
        mg[hit] <- a[hit]/sqrt(c[hit])
        error <- pmin(m3*h^3/6,(mg+grad_abs(x0))*h)
        error[x0==0] <- 0
        sum(error)
      }
    }
    list(log_density=log_density,gradient=gradient,remainder=remainder,third_bound=third_bound)
  })
  if (generic_likelihood) {
    if (is.null(n_obs)) {
      if (!is.null(X)) n_obs <- nrow(as.matrix(X)) else if (!is.null(y)) n_obs <- length(y)
    }
    n <- if (is.null(n_obs)) NA_integer_ else as.integer(n_obs)
    X <- matrix(0,1L,d); y <- 0; offset <- 0
  } else {
    X <- as.matrix(X); n <- nrow(X)
  }
  if (is.null(score)) {
    score <- as.numeric(crossprod(X,y-plogis(X %*% center+offset)))+prior_terms$gradient(center)
  }
  score <- as.numeric(score)
  cutoff_radius_sq <- radius_sq
  if (cutoff=="paper") {
    # The target region is ellipsoidal in z; the bound uses an enclosing ball.
    eig_min <- min(eigen(Sigma,symmetric=TRUE,only.values=TRUE)$values)
    lower <- eig_min-256*.Machine$double.eps*sum(abs(Sigma))
    if (!is.finite(lower) || lower<=0) stop("The covariance is too ill-conditioned for a safe enclosing radius.")
    radius_sq <- cutoff_radius_sq/(n*lower)
    if (!is.finite(radius_sq) || radius_sq<=0) stop("The enclosing radius is not numerically resolvable.")
  }
  radius <- sqrt(radius_sq)
  auto_variance <- identical(proposal_variance,"auto")
  if (generic_likelihood) rotate <- FALSE
  L <- t(chol(Sigma))
  B <- X %*% L
  eta <- as.numeric(X %*% center + offset)
  p <- plogis(eta)
  norm_b <- sqrt(rowSums(B^2))
  c3_at_center <- p*(1-p)*(1-2*p)

  Q <- diag(d)
  if (isTRUE(rotate) && d > 1L && partition == "boxes") {
    v <- as.numeric(crossprod(B, c3_at_center*norm_b^2))
    H <- crossprod(B, B*(abs(c3_at_center)*norm_b))
    if (sqrt(sum(v^2)) > 1e-12) {
      Q0 <- qr.Q(qr(cbind(v, diag(d))))[,seq_len(d),drop=FALSE]
      rest <- Q0[,-1L,drop=FALSE]
      e <- eigen(crossprod(rest, H %*% rest), symmetric=TRUE)
      Q <- cbind(Q0[,1L], rest %*% e$vectors)
    } else {
      Q <- eigen(H, symmetric=TRUE)$vectors
    }
    B <- B %*% Q; L <- L %*% Q
  }
  s <- as.numeric(crossprod(L, score))
  cutoff_info <- list(type=cutoff,radius_sq=cutoff_radius_sq)
  norm_b <- sqrt(rowSums(B^2))
  cap <- radius*norm_b
  bound_scale <- 1+sum(abs(eta)+cap)

  F <- function(t) p*t - 0.5*(softplus(eta+t)-softplus(eta-t))
  Fprime <- function(t) p - 0.5*(plogis(eta+t)+plogis(eta-t))
  delta_base <- function(z) sum(s*z) + sum(F(as.numeric(B %*% z)))
  delta_original <- function(z) {
    shift <- as.numeric(L %*% z)
    lp <- log_post(center+shift); lm <- log_post(center-shift)
    0.5*(lp-lm)
  }
  if (generic_likelihood) {
    raw_bound <- posterior_remainder_bound
    delta_fast <- delta_original
    delta_evaluation <- "log_post"
  } else {
    score_prior <- score-as.numeric(crossprod(X,y-p))
    prior_linear <- as.numeric(crossprod(L,score_prior))
    if (!is.null(prior_terms)) {
      log_prior <- prior_terms$log_density
      row_norm <- sqrt(rowSums(L^2))
      score_error <- sqrt(sum(as.numeric(crossprod(L,prior_terms$gradient(center)-score_prior))^2))
      prior_cubic <- Inf
      if (!is.null(prior_terms$third_bound)) {
        m3 <- prior_terms$third_bound(radius*row_norm)
        prior_cubic <- max(0,largest_eigenvalue_bound(crossprod(L,L*(m3*row_norm))))/6
      }
      prior_remainder_bound <- function(r) {
        cubic <- if (is.finite(prior_cubic)) prior_cubic*r^3 else Inf
        min(prior_terms$remainder(r*row_norm),cubic)+r*score_error
      }
    }
    raw_bound <- if (quadratic_prior) 0 else prior_remainder_bound
    if (!is.null(log_prior)) {
      prior_remainder <- function(z) {
        h <- as.numeric(L %*% z)
        lp <- log_prior(center+h); lm <- log_prior(center-h)
        0.5*(lp-lm)-sum(prior_linear*z)
      }
      delta_fast <- function(z) delta_base(z)+prior_remainder(z)
    } else if (quadratic_prior) {
      delta_fast <- delta_base
    } else {
      delta_fast <- delta_original
      delta_evaluation <- "log_post"
    }
  }
  prior_bound <- local({
    bound <- raw_bound
    padding <- bound_padding
    if (is.numeric(bound)) {
      constant <- bound
      bound <- function(radius) constant
    }
    function(radius) {
      value <- bound(radius)
      if (!is.finite(value) || value<0) stop("Il bound deve essere finito e non negativo.")
      if (value==0) 0 else value+padding*(1+value)
    }
  })
  delta_gradient <- function(z) {
    g <- s+as.numeric(crossprod(B,Fprime(as.numeric(B %*% z))))
    if (!is.null(prior_terms)) {
      h <- as.numeric(L %*% z)
      g <- g+as.numeric(crossprod(L,(prior_terms$gradient(center+h)+
        prior_terms$gradient(center-h))/2))-prior_linear
    }
    g
  }

  # K3 = B' diag(m3_i ||b_i||) B, where m3_i bounds |a'''| over the interval.
  lambda <- 0
  if (!generic_likelihood) {
    logistic_c3 <- function(x) { u <- plogis(x); abs(u*(1-u)*(1-2*u)) }
    ell <- eta-cap; upp <- eta+cap
    c3max <- pmax(logistic_c3(ell),logistic_c3(upp))
    critical <- log(2+sqrt(3))
    crosses <- (ell <= critical & upp >= critical) |
               (ell <= -critical & upp >= -critical)
    c3max[crosses] <- sqrt(3)/18
    K3 <- crossprod(B,B*(c3max*norm_b))
    lambda <- max(0,eigen(K3,symmetric=TRUE,only.values=TRUE)$values[1L])+
      64*.Machine$double.eps*(1+sum(abs(K3)))
  }
  score_norm <- sqrt(sum(s^2))
  cubic_log_bound <- function(r) {
    D <- score_norm*r+lambda*r^3/6+prior_bound(r)
    D <- D+bound_padding*(bound_scale+D)
    M <- if (search == "cosh") log_cosh(D) else D
    M+bound_padding*(1+D)
  }
  global_log_bound <- cubic_log_bound(radius)
  radial_bound <- function(lo,hi) {
    if (lo >= radius) return(0)
    min(global_log_bound,cubic_log_bound(min(hi,radius)))
  }

  radial_grid <- seq(0,radius,length.out=proposal_bins+1L)
  grid_lo <- radial_grid[-length(radial_grid)]
  grid_hi <- radial_grid[-1L]
  grid_M <- vapply(grid_hi,function(r) radial_bound(0,r),numeric(1))
  proposal_objective <- function(c) {
    -0.5*d*log1p(-2*c)+max(grid_M-c*grid_lo^2,-c*radius_sq)
  }
  if (auto_variance) {
    op <- optimize(proposal_objective,c(0,0.5-1e-8),tol=1e-8)
    pc <- if (op$objective < proposal_objective(0)) op$minimum else 0
    proposal_variance <- 1/(1-2*pc)
  }
  proposal_sd <- sqrt(proposal_variance)
  pc <- 0.5*(1-1/proposal_variance)
  proposal_constant <- 0.5*d*log(proposal_variance)
  tilted_radial <- function(lo,hi) {
    if (pc == 0) return(radial_bound(lo,hi))
    if (lo >= radius) return(proposal_constant-pc*lo^2)
    inside_hi <- min(hi,radius)
    ids <- which(grid_hi >= lo & grid_lo <= inside_hi)
    M <- min(max(grid_M[ids]-pc*pmax(grid_lo[ids],lo)^2),
             radial_bound(lo,inside_hi)-pc*lo^2)
    if (hi >= radius) M <- max(M,-pc*radius_sq)
    proposal_constant+M+bound_padding*(1+abs(M))
  }
  base_radial_bound <- function(lo,hi) tilted_radial(proposal_sd*lo,proposal_sd*hi)

  root_log_bound <- base_radial_bound(0,Inf)
  if (length(root_log_bound)!=1L || !is.finite(root_log_bound))
    stop("Invalid global correction bound.")
  list(center=center,Sigma=Sigma,L=L,d=d,n=n,radius=radius,radius_sq=radius_sq,
    cutoff=cutoff_info,log_post=log_post,engine="global_rejection",search=search,
    proposal_sd=proposal_sd,proposal_variance=proposal_variance,
    proposal_c=pc,proposal_constant=proposal_constant,
    envelope_geometry=list(B=B,s=s,eta=eta,norm_b=norm_b),
    coefficient_bins=coefficient_bins,bound_padding=bound_padding,
    base_radial_bound=base_radial_bound,
    delta=if (delta_evaluation=="analytic") delta_fast else delta_original,
    delta_gradient=delta_gradient,prior_bound=prior_bound,
    root_log_bound=root_log_bound,log_bound_globale=root_log_bound)
}

## 3) Rejection sampling and flipping ===================================

branch_log_envelope <- function(br,tt,w) {
  rho2 <- sum(w*w)
  value <- br$constant+br$bp*tt+.5*br$ap*tt^2+.5*br$ar*rho2
  if (!identical(br$kind,"acg")) return(value+sum(br$gperp*w))
  k <- length(w); u <- if (rho2>0) w/sqrt(rho2) else c(1,rep(0,k-1L))
  quad <- sum(br$omega*as.numeric(crossprod(br$V,u))^2)
  value-.5*k*log(quad)+log(2)-softplus(-2*sum(br$gperp*w))
}

draw_directional_cell <- function(node) {
  lws <- vapply(node$branches,function(x) x$logweight,numeric(1))
  br <- node$branches[[which.max(lws-log(rexp(length(lws))))]]
  ni <- br$normal; lp <- log_add(ni$anchor,log(runif(1))+ni$logmass)
  tt <- br$mp+sqrt(br$vp)*qnorm(lp,lower.tail=!ni$upper,log.p=TRUE)
  ri <- br$radial; lp <- log_add(ri$anchor,log(runif(1))+ri$lm)
  rr <- if (isTRUE(ri$stable)) stable_shell_quantile(lp,ri) else if (ri$ncp==0) qchisq(lp,ri$d,lower.tail=!ri$upper,log.p=TRUE) else
    qchisq(lp,ri$d,ncp=ri$ncp,lower.tail=!ri$upper,log.p=TRUE)
  rho <- sqrt(ri$variance*rr)
  if (identical(br$kind,"acg")) {
    u <- as.numeric(br$V %*% (rnorm(node$d-1L)/sqrt(br$omega)))
    nu <- sqrt(sum(u*u))
    u <- u/nu
    if (runif(1)>plogis(2*rho*sum(br$gperp*u))) u <- -u
    angular <- u
  } else {
    gn <- sqrt(sum(br$gperp^2))
    u <- if (gn>0) br$gperp/gn else c(1,rep(0,node$d-2L))
    angular <- draw_vmf(u,rho*gn)
  }
  v <- c(tt,rho*angular)
  if (any(!is.finite(v)) || tt<node$left || tt>node$right ||
      rho<node$rho_lo || rho>node$rho_hi) stop("Sector inversion failed numerically.")

  logenv <- log_sum(vapply(node$branches,function(x)
    branch_log_envelope(x,tt,v[-1L]),numeric(1)))
  list(z=as.numeric(node$Q %*% v),r2=tt^2+rho^2,logenv=logenv)
}

draw_shell_candidate <- function(model) {
  tree <- model$shell_choice_tree
  if (is.null(tree)) {
    j <- which.max(model$shell_logweights-log(rexp(length(model$shell_nodes))))
  } else {
    i <- 1L
    while (i<tree$base) {
      left <- 2*i; right <- left+1L
      small <- if (tree$mass[left]<=tree$mass[right]) left else right
      big <- if (small==left) right else left
      i <- if (log(runif(1))<tree$mass[small]-tree$mass[i]) small else big
    }
    j <- i-tree$base+1L
  }
  node <- model$shell_nodes[[j]]; delta <- NA_real_; logP <- 0
  if (node$kind=="directional") {
    draw <- draw_directional_cell(node); z <- draw$z

    if (draw$r2<node$lo^2 || draw$r2>=node$hi^2) {
      return(list(z=z,delta=delta,log_accept=-Inf))
    }
    if (inside_region(model,z)) {
      delta <- model$delta(z); logP <- log_cosh(delta)
    }
    la <- logP-draw$logenv
  } else {
    lp <- log_add(node$anchor,log(runif(1))+node$lm)
    rr <- if (isTRUE(node$stable)) stable_shell_quantile(lp,node) else if (node$ncp == 0) qchisq(lp,node$d,lower.tail=!node$upper,log.p=TRUE) else
      qchisq(lp,node$d,ncp=node$ncp,lower.tail=!node$upper,log.p=TRUE)
    r <- sqrt(node$variance*rr)
    if (!is.finite(r) || r < node$lo || r > node$hi) stop("Shell inversion failed numerically.")
    gn <- sqrt(sum(node$g^2))
    if (gn == 0) {
      direction <- rep(0,node$d); direction[1L] <- 1
    } else direction <- node$g/gn
    angular <- draw_vmf(direction,r*gn)
    z <- r*angular
    if (gn > 0 && runif(1) < 0.5) z <- -z
    if (inside_region(model,z)) { delta <- model$delta(z); logP <- log_cosh(delta) }
    la <- switch(node$kind,trace=logP-node$err-log_cosh(sum(node$g*z)),
      gaussian=logP+.5*model$d*log(node$variance)-.5*(1-1/node$variance)*sum(z*z)-node$err,
      exterior=logP)
  }
  if (length(la)!=1L || !is.finite(la) || la>0) stop("Directional envelope violated; no valid output.")
  list(z=z,delta=delta,log_accept=la)
}

draw_candidate <- function(model) {
  if (model$engine %in% c("shell_rejection","directional_shell_rejection"))
    return(draw_shell_candidate(model))
  z <- model$proposal_sd*rnorm(model$d); r2 <- sum(z*z)
  delta <- NA_real_; logP <- 0
  if (inside_region(model,z)) {
    delta <- model$delta(z)
    if (length(delta)!=1L || !is.finite(delta)) stop("Nonfinite delta.")
    logP <- if (model$search=="cosh") log_cosh(delta) else delta
  }
  la <- logP+model$proposal_constant-model$proposal_c*r2-model$root_log_bound
  if (length(la)!=1L || !is.finite(la) || la>0) stop("The proposal envelope was violated.")
  list(z=z,delta=delta,log_accept=la)
}

sample_skew_cosh_rejection <- function(N_keep,model,max_evals_per_draw=50000L,
                                      max_seconds_per_draw=Inf,verbose=FALSE) {
  start <- proc.time()[["elapsed"]]
  sim_bar <- sim_skew <- matrix(NA_real_,N_keep,model$d)
  weights <- weights_at_sim_skew <- numeric(N_keep)
  for (k in seq_len(N_keep)) {
    t0 <- proc.time()[["elapsed"]]; ne <- 0L
    repeat {
      if (ne>=max_evals_per_draw || proc.time()[["elapsed"]]-t0>max_seconds_per_draw)
        stop(sprintf("Rejection budget reached at draw %d.",k))
      x <- draw_candidate(model)
      ne <- ne+1L
      if (log(runif(1))<=x$log_accept) break
    }
    if (model$search=="exp" && runif(1)<0.5) {
      x$z <- -x$z
      if (!is.na(x$delta)) x$delta <- -x$delta
    }
    if (is.na(x$delta)) x$delta <- model$delta(x$z)
    if (length(x$delta)!=1L || !is.finite(x$delta)) stop("Nonfinite flipping delta.")
    shift <- as.numeric(model$L %*% x$z); w <- plogis(2*x$delta)
    retain <- runif(1)<=w
    sim_bar[k,] <- model$center+shift
    sim_skew[k,] <- model$center+if (retain) shift else -shift
    weights[k] <- w
    weights_at_sim_skew[k] <- if (retain) w else 1-w
  }
  elapsed <- proc.time()[["elapsed"]]-start
  colnames(sim_bar) <- colnames(sim_skew) <- names(model$center)
  if (verbose) cat(sprintf("Sampling: time %.2f s\n",elapsed))
  # weights refers to sim_bar; weights_at_sim_skew refers to the flipped draws.
  list(sim_bar=sim_bar,sim_skew=sim_skew,weights=weights,
    weights_at_sim_skew=weights_at_sim_skew,elapsed_sec=elapsed)
}

## 4) Radial shells =====================================================

poisson_chisq_logtail <- function(q,d,ncp,upper=FALSE,rtol=1e-13,max_terms=100000L) {

  if (q<=0) return(if (upper) 0 else -Inf)
  if (is.infinite(q)) return(if (upper) -Inf else 0)
  if (ncp==0) return(pchisq(q,d,lower.tail=!upper,log.p=TRUE))
  total <- -Inf; first <- 0L; mu <- ncp/2
  repeat {
    last <- min(first+63L,max_terms-1L); jj <- first:last
    vals <- dpois(jj,mu,log=TRUE)+pchisq(q,d+2*jj,lower.tail=!upper,log.p=TRUE)
    mm <- max(vals)
    block <- if (is.finite(mm)) mm+log(sum(exp(vals-mm))) else -Inf
    total <- log_add(total,block)
    rest <- ppois(last,mu,lower.tail=FALSE,log.p=TRUE)

    if (!upper) rest <- rest+pchisq(q,d+2*(last+1L),log.p=TRUE)
    if (is.finite(total) && rest<=total+log(rtol)) return(min(0,total))
    if (last>=max_terms-1L) stop("Noncentral chi-square series exceeded its evaluation budget.")
    first <- last+1L
  }
}

noncentral_shell <- function(lo,hi,d,ncp=0,variance=1,stable=TRUE) {
  x <- lo^2/variance; y <- hi^2/variance
  upper <- x >= d+ncp
  stable <- isTRUE(stable) && ncp>=80
  prob <- function(q) {
    if (stable) return(poisson_chisq_logtail(q,d,ncp,upper))
    if (ncp == 0) pchisq(q,d,lower.tail=!upper,log.p=TRUE) else
      pchisq(q,d,ncp=ncp,lower.tail=!upper,log.p=TRUE)
  }
  big <- prob(if (upper) x else y)
  small <- prob(if (upper) y else x)
  lm <- log_difference(big,small)
  if (!is.finite(lm)) stop("Shell mass is not numerically resolvable.")
  list(lo=lo,hi=hi,d=d,ncp=ncp,variance=variance,lm=lm,anchor=small,upper=upper,stable=stable)
}

stable_shell_quantile <- function(lp,ri) {
  lo <- ri$lo^2/ri$variance; hi <- ri$hi^2/ri$variance
  ff <- function(x) {
    ans <- poisson_chisq_logtail(x,ri$d,ri$ncp,ri$upper)-lp
    if (is.infinite(ans)) return(sign(ans)*1e300)
    ans
  }
  root <- uniroot(ff,c(lo,hi),tol=128*.Machine$double.eps*(1+hi),
                   maxiter=200L,check.conv=TRUE)$root
  if (!is.finite(root) || root<lo || root>hi ||
      abs(ff(root))>1e-9*(1+abs(lp))) stop("Stable chi-square inversion failed numerically.")
  root
}

draw_vmf <- function(direction,kappa) {

  d <- length(direction)
  if (kappa == 0) {
    z <- rnorm(d); nz <- sqrt(sum(z*z))
    return(z/nz)
  }
  if (d == 1L) {
    z <- direction*if (runif(1) <= plogis(2*kappa)) 1 else -1
    return(z)
  }
  dd <- d-1
  b <- dd/(2*kappa+sqrt(4*kappa^2+dd^2))
  if (!is.finite(b) || b <= 0) stop("vMF concentration is not numerically resolvable.")
  for (j in seq_len(10000L)) {
    zz <- rbeta(1,dd/2,dd/2)
    den <- 1-zz+b*zz
    w <- (1-zz-b*zz)/den
    logacc <- 2*kappa*b*(1-2*zz)/((1+b)*den)+dd*(log1p(b)-log(2)-log(den))
    if (log(runif(1)) <= logacc) break
    if (j == 10000L) stop("Internal vMF budget reached; no provisional sample returned.")
  }
  v <- rnorm(dd); nv <- sqrt(sum(v*v))
  u <- c(w,sqrt((1-w)*(1+w))*v/nv)
  e1 <- c(1,rep(0,dd)); h <- e1-direction; nh <- sqrt(sum(h*h))
  if (nh > 0) { h <- h/nh; u <- u-2*h*sum(h*u) }
  u
}

make_skew_cosh_shells <- function(model,n_shells=64L,n_bound_radii=12L,verbose=FALSE) {
  start <- proc.time()[["elapsed"]]
  gg <- model$envelope_geometry; d <- model$d; R <- model$radius
  radii <- R*sqrt(seq_len(n_bound_radii)/n_bound_radii)
  anchors <- vector("list",n_bound_radii)
  for (j in seq_len(n_bound_radii)) {
    anchors[[j]] <- trace_envelope(gg$B,gg$s,gg$eta,gg$norm_b,radii[j],model$coefficient_bins)
  }
  grid <- R*sqrt(seq(0,1,length.out=n_shells+1L))
  nodes <- vector("list",n_shells+1L)
  for (j in seq_len(n_shells)) {
    a <- grid[j]; b <- grid[j+1L]
    ii <- which(radii >= b)[1L]
    te <- anchors[[ii]]
    g <- te$linear*b^2
    err <- te$score_norm*b+te$cubic_remainder*b^3
    err <- err+model$bound_padding*(te$numerical_scale+1+abs(err)+sum(abs(g)))+model$prior_bound(b)
    node0 <- noncentral_shell(a,b,d,variance=model$proposal_variance)
    node0$g <- rep(0,d); node0$kind <- "gaussian"
    node0$err <- min(model$root_log_bound,model$base_radial_bound(a/model$proposal_sd,b/model$proposal_sd))
    node0$logweight <- node0$err+node0$lm

    node1 <- tryCatch(noncentral_shell(a,b,d,ncp=sum(g*g)),error=function(e) NULL)
    if (is.null(node1)) {
      nodes[[j]] <- node0
    } else {
      node1$g <- g; node1$kind <- "trace"; node1$err <- err
      node1$logweight <- err+0.5*sum(g*g)+node1$lm
      nodes[[j]] <- if (node1$logweight < node0$logweight) node1 else node0
    }
  }
  tail <- noncentral_shell(R,Inf,d)
  tail$g <- rep(0,d); tail$kind <- "exterior"; tail$err <- 0; tail$logweight <- tail$lm
  nodes[[n_shells+1L]] <- tail
  out <- update_shells(model,nodes)
  out$engine <- "shell_rejection"
  print_refinement("Shells",model$root_log_bound,out$root_log_bound,start,verbose)
  invisible(out)
}

make_refinement_control <- function(nodes,threshold,patience) {
  lo <- vapply(nodes,function(x) x$lo,numeric(1))
  hi <- vapply(nodes,function(x) x$hi,numeric(1))
  ord <- order(lo,hi); n <- length(ord)
  first <- c(TRUE,lo[ord[-1L]]!=lo[ord[-n]] | hi[ord[-1L]]!=hi[ord[-n]])
  groups <- integer(n); groups[ord] <- cumsum(first); ng <- sum(first)
  list(groups=groups,stale=integer(ng),stopped=rep(FALSE,ng),
    threshold=threshold,patience=patience)
}

remaining_gain_bound <- function(log_remaining,log_total) {
  if (is.infinite(log_remaining) && log_remaining<0) return(0)
  x <- log_remaining-log_total
  if (!is.finite(x) || x>=0) return(Inf)

  if (x < -log(2)) -log1p(-exp(x)) else -log(-expm1(x))
}

update_refinement_control <- function(control,group,gain,log_remaining=Inf,log_total=0) {
  if (control$threshold==0) return(control)
  potential <- remaining_gain_bound(log_remaining,log_total)
  control$stale[group] <- if (gain<control$threshold) control$stale[group]+1L else 0L
  margin <- 64*.Machine$double.eps*(1+abs(potential))
  small <- is.finite(potential) && potential+margin<control$threshold
  control$stopped[group] <- control$stale[group]>=control$patience && small
  control
}

update_refinement_from_nodes <- function(control,group,gain,nodes,attempted,log_total) {
  pending <- !attempted & control$groups==group &
    vapply(nodes,function(x) is.finite(x$hi),logical(1))
  log_remaining <- if (any(pending)) log_sum(
    envelope_log_masses(nodes[pending])) else -Inf
  update_refinement_control(control,group,gain,log_remaining,log_total)
}

print_refinement <- function(label,before,after,start,verbose) {
  if (verbose) cat(sprintf("%s: log C improvement %.6g; time %.2f s\n",
    label,before-after,proc.time()[["elapsed"]]-start))
}

refine_skew_cosh_shells <- function(model,n_splits=8L,verbose=FALSE,
                                    coefficient_function=cubic_coefficients,
                                    min_log_bound_improvement=0,stagnation_patience=1L) {
  start <- proc.time()[["elapsed"]]
  nodes <- model$shell_nodes; gg <- model$envelope_geometry; d <- model$d
  control <- make_refinement_control(nodes,min_log_bound_improvement,stagnation_patience)
  current_B <- model$root_log_bound
  attempted <- rep(FALSE,length(nodes))
  for (step in seq_len(n_splits)) {
    eligible <- which(!attempted & !control$stopped[control$groups] &
      vapply(nodes,function(x) is.finite(x$hi),logical(1)))
    if (!length(eligible)) break
    j <- eligible[which.max(envelope_log_masses(nodes[eligible]))]
    group <- control$groups[j]
    parent <- nodes[[j]]; middle <- sqrt(0.5*(parent$lo^2+parent$hi^2))
    if (!(middle > parent$lo && middle < parent$hi)) stop("Shell cannot be subdivided further in double precision.")
    ends <- c(parent$lo,middle,parent$hi); children <- vector("list",2L)
    for (i in 1:2) {
      a <- ends[i]; b <- ends[i+1L]

      old <- noncentral_shell(a,b,d,ncp=parent$ncp,variance=parent$variance)
      old$g <- parent$g; old$err <- parent$err; old$kind <- parent$kind
      old$logweight <- parent$logweight-parent$lm+old$lm
      te <- trace_envelope(gg$B,gg$s,gg$eta,gg$norm_b,b,
                              model$coefficient_bins,coefficient_function)
      g <- te$linear*b^2
      err <- te$score_norm*b+te$cubic_remainder*b^3
      err <- err+model$bound_padding*(te$numerical_scale+1+abs(err)+sum(abs(g)))+model$prior_bound(b)
      fresh <- tryCatch(noncentral_shell(a,b,d,ncp=sum(g*g)),error=function(e) NULL)
      if (!is.null(fresh)) {
        fresh$g <- g; fresh$kind <- "trace"; fresh$err <- err
        fresh$logweight <- err+0.5*sum(g*g)+fresh$lm
        if (fresh$logweight < old$logweight) old <- fresh
      }
      children[[i]] <- old
    }
    candidate_nodes <- c(nodes[-j],children)
    B <- log_sum(envelope_log_masses(candidate_nodes))
    previous_B <- current_B
    if (B<=current_B) {
      nodes <- candidate_nodes
      attempted <- c(attempted[-j],rep(FALSE,length(children)))
      control$groups <- c(control$groups[-j],rep(group,length(children)))
      current_B <- B
    } else attempted[j] <- TRUE
    gain <- previous_B-current_B
    control <- update_refinement_from_nodes(control,group,gain,nodes,attempted,current_B)
  }
  out <- update_shells(model,nodes)
  print_refinement("Radial",model$root_log_bound,out$root_log_bound,start,verbose)
  invisible(out)
}

## 5) Directional refinement ============================================

integrated_cubic_coefficients <- function(eta,cap,bins) {

  lc <- log_cosh(eta); amp <- log(abs(tanh(eta/2)))-log(4)
  ap <- rep(0,length(eta)); dp <- log_add(lc,0)
  lower <- upper <- 2*exp(amp-dp)
  il <- iu <- numeric(length(eta))
  for (j in seq_len(bins)) {
    t <- cap*j/bins; an <- numeric(length(t)); nz <- t > 0
    an[nz] <- t[nz]+2*log(-expm1(-t[nz]))-2*log(t[nz])
    dn <- log_add(lc,log_cosh(t))
    ql <- exp(amp+ap-dn); qu <- exp(amp+an-dp)
    volume <- (j^3-(j-1)^3)/bins^3
    il <- il+2*ql*volume; iu <- iu+2*qu*volume
    lower <- pmin(lower,il/(j/bins)^3); upper <- pmax(upper,iu/(j/bins)^3)
    ap <- an; dp <- dn
  }
  pad <- 512*.Machine$double.eps*(1+upper)*bins
  lower <- pmax(0,lower-pad); upper <- upper+pad
  if (any(!is.finite(c(lower,upper)))) stop("Integrated coefficients are not numerically resolvable.")
  list(c3=-sign(eta)*(lower+upper)/2,remainder=(upper-lower)/2)
}

quadratic_remainder_coefficients <- function(eta,cap,bins) {

  n <- length(eta); lc <- log_cosh(eta)
  amp <- log(abs(tanh(eta/2)))-log(4)
  ap <- rep(0,n); dp <- log_add(lc,0); il <- iu <- numeric(n)
  pl <- pu <- 2*exp(amp-dp)
  lower <- upper <- weight <- matrix(0,bins,n)
  for (j in seq_len(bins)) {
    t <- cap*j/bins; an <- numeric(n); nz <- t>0
    an[nz] <- t[nz]+2*log(-expm1(-t[nz]))-2*log(t[nz])
    dn <- log_add(lc,log_cosh(t))
    ql <- exp(amp+ap-dn); qu <- exp(amp+an-dp)
    volume <- (j^3-(j-1)^3)/bins^3
    il <- il+2*ql*volume; iu <- iu+2*qu*volume
    nl <- il/(j/bins)^3; nu <- iu/(j/bins)^3
    lower[j,] <- if (j==1L) nl else pmin(pl,nl)
    upper[j,] <- if (j==1L) nu else pmax(pu,nu)
    weight[j,] <- t
    pl <- nl; pu <- nu; ap <- an; dp <- dn
  }
  pad <- 512*.Machine$double.eps*(1+apply(upper,2L,max))*bins
  lower <- sweep(lower,2L,pad,"-"); lower[lower<0] <- 0
  upper <- sweep(upper,2L,pad,"+")
  lo <- apply(lower,2L,min); hi <- apply(upper,2L,max)

  for (j in seq_len(48L)) {
    cc <- (lo+hi)/2
    ep <- apply(weight*sweep(upper,2L,cc,"-"),2L,max)
    em <- apply(weight*(-sweep(lower,2L,cc,"-")),2L,max)
    move <- ep>em; lo[move] <- cc[move]; hi[!move] <- cc[!move]
  }
  cc <- (lo+hi)/2
  ep <- apply(weight*sweep(upper,2L,cc,"-"),2L,max)
  em <- apply(weight*(-sweep(lower,2L,cc,"-")),2L,max)
  re <- pmax(ep,em)+pad*(1+cap)
  if (any(!is.finite(c(cc,re))) || any(re<0)) stop("Quadratic remainder is not numerically resolvable.")
  list(c3=-sign(eta)*cc,remainder=re)
}

direction_basis <- function(u) {
  u <- u/sqrt(sum(u*u)); Q <- diag(length(u)); h <- Q[,1L]-u
  nh <- sqrt(sum(h*h))
  if (nh > 0) { h <- h/nh; Q <- Q-2*tcrossprod(h) }
  Q
}

skew_cosh_direction <- function(model,radius=model$radius,max_iter=80L) {
  gg <- model$envelope_geometry; B <- gg$B; eta <- gg$eta; p <- plogis(eta)
  w <- p*(1-p); el <- eta-radius*gg$norm_b; eu <- eta+radius*gg$norm_b
  wl <- plogis(el)*(1-plogis(el)); wu <- plogis(eu)*(1-plogis(eu))
  wmax <- pmax(wl,wu); wmax[el<=0 & eu>=0] <- .25
  variation <- pmax(wmax-w,w-pmin(wl,wu))
  info <- crossprod(B,B*variation)
  eig <- eigen(info,symmetric=TRUE)
  trace <- as.numeric(crossprod(B,w*(1-2*p)*gg$norm_b^2))
  starts <- list(eig$vectors[,1L])
  if (sqrt(sum(trace^2)) > 1e-12) starts <- c(starts,list(trace/sqrt(sum(trace^2))))
  best <- -Inf; direction <- starts[[1L]]
  for (u in starts) {
    value <- model$delta(radius*u)
    if (value<0) { u <- -u; value <- -value }
    for (it in seq_len(max_iter)) {
      grad <- radius*model$delta_gradient(radius*u)
      tangent <- grad-u*sum(u*grad); normg <- sqrt(sum(tangent^2))
      if (normg < 1e-7) break
      step <- 1/max(1,normg); moved <- FALSE
      for (ls in seq_len(20L)) {
        v <- u+step*tangent; v <- v/sqrt(sum(v*v))
        nv <- model$delta(radius*v)
        if (is.finite(nv) && nv > value+1e-4*step*normg^2) {
          u <- v; value <- nv; moved <- TRUE; break
        }
        step <- step/2
      }
      if (!moved) break
    }
    if (value>best) { best <- value; direction <- u }
  }
  list(direction=direction)
}

directional_geometry <- function(model,r,Q,coefficient_function=integrated_cubic_coefficients,quadratic_remainder=FALSE) {
  gg <- model$envelope_geometry; A <- gg$B %*% Q; d <- model$d; k <- d-1L
  co <- coefficient_function(gg$eta,r*gg$norm_b,model$coefficient_bins)
  T <- cubic_tensor(A,co$c3)
  D <- array(T[-1L,-1L,-1L,drop=FALSE],c(k,k,k))
  harmonic <- harmonic_decomposition(D,sum(abs(T)))
  re <- max(0,largest_eigenvalue_bound(crossprod(A,A*(co$remainder*gg$norm_b))))
  ce <- eigen(matrix(T[1L,-1L,-1L],k,k),symmetric=TRUE,only.values=TRUE)$values
  ep <- 128*.Machine$double.eps*(1+sum(abs(T)))
  out <- list(a=T[1L,1L,1L],b=as.numeric(T[1L,1L,-1L]),v=harmonic$v,
       cmin=min(ce)-ep,cmax=max(ce)+ep,hn=harmonic$norm,re=re,
       s=as.numeric(crossprod(Q,gg$s)),Q=Q,prior_error=model$prior_bound(r),
       padding=model$bound_padding*(1+sum(abs(T))+sum(abs(gg$eta))+r*sum(gg$norm_b)))
  if (isTRUE(quadratic_remainder)) {
    coq <- quadratic_remainder_coefficients(gg$eta,r*gg$norm_b,model$coefficient_bins)
    Tq <- cubic_tensor(A,coq$c3)
    Dq <- array(Tq[-1L,-1L,-1L,drop=FALSE],c(k,k,k))
    harmonic_q <- harmonic_decomposition(Dq,sum(abs(Tq)))

    M <- crossprod(A,A*(coq$remainder/3))
    M <- (M+t(M))/2+diag(128*.Machine$double.eps*(1+sum(abs(M))),d)
    out$quadratic <- list(a=Tq[1L,1L,1L],b=as.numeric(Tq[1L,1L,-1L]),
      C=matrix(Tq[1L,-1L,-1L],k,k),v=harmonic_q$v,hn=harmonic_q$norm,M=M,s=out$s,prior_error=out$prior_error,
      padding=model$bound_padding*(1+sum(abs(Tq))+sum(abs(M))+sum(abs(gg$eta))+r*sum(gg$norm_b)))
  }
  out
}

gaussian_branch <- function(cell,constant,bp,gperp,ap=0,ar=0,stable_radial=TRUE) {
  vp <- 1/(1-ap); vr <- 1/(1-ar); mp <- vp*bp; mr <- vr*gperp
  ni <- normal_interval_mass((cell$left-mp)/sqrt(vp),(cell$right-mp)/sqrt(vp))
  ri <- noncentral_shell(cell$rho_lo,cell$rho_hi,cell$d-1L,
                             ncp=sum(mr^2)/vr,variance=vr,stable=stable_radial)
  norm <- constant+.5*log(vp)+.5*vp*bp^2+.5*(cell$d-1L)*log(vr)+.5*vr*sum(gperp^2)
  lw <- norm+ni$logmass+ri$lm
  if (!is.finite(lw)) stop("Directional branch mass is not numerically resolvable.")
  list(constant=constant,bp=bp,gperp=gperp,ap=ap,ar=ar,
       vp=vp,mp=mp,normal=ni,radial=ri,logweight=lw)
}

best_gaussian_branch <- function(cell,cc,gt,g,ap0,ar0) {
  tc <- (cell$left+cell$right)/2
  hw <- (cell$right-cell$left)/2
  best <- NULL
  for (cap in c(0,.5,.9,.99)) {
    ap <- min(ap0,cap); ar <- min(ar0,cap)
    constant <- cc+.5*(ap0-ap)*hw^2+.5*(ar0-ar)*cell$rho_hi^2-gt*tc+.5*ap*tc^2-log(2)
    x <- tryCatch(gaussian_branch(cell,constant,gt-ap*tc,g,ap,ar),error=function(e) NULL)
    if (!is.null(x) && (is.null(best) || x$logweight<best$logweight)) best <- x
  }
  best
}

quadratic_envelope_branches <- function(cell,cg) {
  if (is.null(cg)) return(NULL)
  l <- cell$left; h <- cell$right; tc <- (l+h)/2; hw <- (h-l)/2
  k <- cell$d-1L; rl2 <- cell$rho_lo^2; rh2 <- cell$rho_hi^2; rc2 <- (rl2+rh2)/2
  M <- cg$M; mp <- as.numeric(M[1L,-1L]); Mpp <- matrix(M[-1L,-1L],k,k)
  gp <- cg$s[-1L]-.5*tc^2*cg$b-cg$v*rc2/(2*(k+2))
  t2diff <- max(abs(l^2-tc^2),abs(h^2-tc^2),if (l<=0 && h>=0) tc^2 else 0)
  harmonic <- cg$hn*cell$rho_hi^3/6
  err <- .5*t2diff*sqrt(sum(cg$b^2))*cell$rho_hi+
    sqrt(sum(cg$v^2))*(rh2-rl2)*cell$rho_hi/(4*(k+2))+
    harmonic+abs(cg$a)*hw^3/6+hw*sqrt(sum(mp^2))*cell$rho_hi
  err <- err+cg$padding*(1+cell$hi^3+abs(err))+cg$prior_error
  branches <- vector("list",2L)
  for (i in 1:2) {
    sg <- c(1,-1)[i]
    ap0 <- -sg*cg$a*tc+M[1L,1L]

    ar0 <- max(largest_eigenvalue_bound(Mpp-sg*l*cg$C),largest_eigenvalue_bound(Mpp-sg*h*cg$C))
    gt <- sg*(cg$s[1L]-.5*cg$a*tc^2)+M[1L,1L]*tc
    cc <- sg*(cg$s[1L]*tc-cg$a*tc^3/6)+.5*M[1L,1L]*tc^2+err
    best <- best_gaussian_branch(cell,cc,gt,sg*gp+tc*mp,ap0,ar0)
    if (is.null(best)) return(NULL)
    branches[[i]] <- best
  }
  branches
}

make_directional_cell <- function(parent,left,right,cg) {
  d <- nrow(cg$Q); k <- d-1L; a <- parent$lo; b <- parent$hi
  near <- max(left,min(0,right)); far <- max(abs(left),abs(right))
  cell <- list(lo=a,hi=b,left=left,right=right,d=d,Q=cg$Q,
    rho_lo=sqrt(max(0,a*a-far*far)),rho_hi=sqrt(max(0,b*b-near*near)),
    kind="directional",geometry=cg)

  candidates <- list()
  if (parent$kind=="directional") {
    inherited <- tryCatch(lapply(parent$branches,function(x)
      if (identical(x$kind,"acg")) acg_branch_mass(x,cell) else
        gaussian_branch(cell,x$constant,x$bp,x$gperp,x$ap,x$ar)),error=function(e) NULL)
  } else if (parent$kind=="trace") {
    g <- as.numeric(crossprod(cg$Q,parent$g))
    inherited <- tryCatch(lapply(c(1,-1),function(sg)
      gaussian_branch(cell,parent$err-log(2),sg*g[1L],sg*g[-1L])),error=function(e) NULL)
  } else {
    h <- 1-1/parent$variance
    inherited <- tryCatch(list(gaussian_branch(cell,parent$err-.5*d*log(parent$variance),
                                0,rep(0,k),h,h)),error=function(e) NULL)
  }
  if (!is.null(inherited)) candidates <- list(inherited)
  tc <- (left+right)/2; hw <- (right-left)/2
  rl2 <- cell$rho_lo^2; rh2 <- cell$rho_hi^2; rc2 <- (rl2+rh2)/2
  gp <- cg$s[-1L]-.5*tc^2*cg$b-cg$v*rc2/(2*(k+2))
  t2diff <- max(abs(left^2-tc^2),abs(right^2-tc^2),if (left<=0 && right>=0) tc^2 else 0)
  err <- .5*t2diff*sqrt(sum(cg$b^2))*cell$rho_hi+
    sqrt(sum(cg$v^2))*(rh2-rl2)*cell$rho_hi/(4*(k+2))+
    (cg$hn*cell$rho_hi^3+cg$re*b^3+abs(cg$a)*hw^3)/6
  err <- err+cg$padding*(1+b^3+abs(err))+cg$prior_error
  branches <- vector("list",2L)
  for (i in 1:2) {
    sg <- c(1,-1)[i]
    ap0 <- -sg*cg$a*tc
    ar0 <- max(-sg*left*cg$cmin,-sg*left*cg$cmax,-sg*right*cg$cmin,-sg*right*cg$cmax)
    gt <- sg*(cg$s[1L]-.5*cg$a*tc^2)
    cc <- sg*(cg$s[1L]*tc-cg$a*tc^3/6)+err
    best <- best_gaussian_branch(cell,cc,gt,sg*gp,ap0,ar0)
    if (is.null(best)) { branches <- NULL; break }
    branches[[i]] <- best
  }
  if (!is.null(branches)) {
    candidates <- c(candidates,list(branches))
  }
  qb <- quadratic_envelope_branches(cell,cg$quadratic)
  if (!is.null(qb)) {
    candidates <- c(candidates,list(qb))
  }
  if (!length(candidates)) return(NULL)
  masses <- vapply(candidates,function(cc) log_sum(envelope_log_masses(cc)),numeric(1))
  cell$branches <- candidates[[which.min(masses)]]; cell$logweight <- min(masses)
  cell
}

make_skew_cosh_directional <- function(model,n_radial_splits=16L,max_steps=32L,
                                        angular_bins=32L,integrated_bound=TRUE,
                                        directional=TRUE,verbose=TRUE,
                                        quadratic_remainder=FALSE,
                                        min_log_bound_improvement=0,stagnation_patience=1L) {
  start <- proc.time()[["elapsed"]]; before <- model$root_log_bound
  if (model$engine=="global_rejection") model <- make_skew_cosh_shells(model)
  cf <- if (isTRUE(integrated_bound)) integrated_cubic_coefficients else cubic_coefficients
  if (n_radial_splits>0) model <- refine_skew_cosh_shells(model,n_radial_splits,FALSE,cf,
    min_log_bound_improvement=min_log_bound_improvement,stagnation_patience=stagnation_patience)

  nodes <- model$shell_nodes; gg <- model$envelope_geometry
  if (isTRUE(integrated_bound)) for (j in seq_along(nodes)) {
    old <- nodes[[j]]
    if (!is.finite(old$hi)) next
    te <- trace_envelope(gg$B,gg$s,gg$eta,gg$norm_b,old$hi,model$coefficient_bins,cf)
    g <- te$linear*old$hi^2; err <- te$score_norm*old$hi+te$cubic_remainder*old$hi^3
    err <- err+model$bound_padding*(te$numerical_scale+1+abs(err)+sum(abs(g)))+model$prior_bound(old$hi)
    x <- tryCatch(noncentral_shell(old$lo,old$hi,model$d,ncp=sum(g*g)),error=function(e) NULL)
    if (!is.null(x)) {
      x$g <- g; x$err <- err; x$kind <- "trace"
      x$logweight <- err+.5*sum(g*g)+x$lm
      if (x$logweight<old$logweight) nodes[[j]] <- x
    }
  }
  control <- make_refinement_control(nodes,min_log_bound_improvement,stagnation_patience)
  attempted <- rep(FALSE,length(nodes))
  if (isTRUE(directional) && max_steps>0) for (step in seq_len(max_steps)) {
    eligible <- which(!attempted & !control$stopped[control$groups] &
      vapply(nodes,function(x) is.finite(x$hi),logical(1)))
    if (!length(eligible)) break
    j <- eligible[which.max(envelope_log_masses(nodes[eligible]))]
    group <- control$groups[j]
    parent <- nodes[[j]]
    if (parent$kind=="directional") {
      cg <- parent$geometry; grid <- seq(parent$left,parent$right,length.out=3L)
    } else {
      dr <- skew_cosh_direction(model,parent$hi)
      cg <- directional_geometry(model,parent$hi,direction_basis(dr$direction),cf,quadratic_remainder)
      grid <- seq(-parent$hi,parent$hi,length.out=angular_bins+1L)
    }
    children <- lapply(seq_len(length(grid)-1L),function(i) make_directional_cell(parent,grid[i],grid[i+1L],cg))

    cm <- if (any(vapply(children,is.null,logical(1)))) Inf else
      log_sum(envelope_log_masses(children))
    current_before <- log_sum(envelope_log_masses(nodes))
    candidate_B <- if (is.finite(cm)) log_sum(c(envelope_log_masses(nodes[-j]),cm)) else Inf
    accepted <- is.finite(cm) && cm<parent$logweight-1e-8*(1+abs(parent$logweight)) &&
      candidate_B<=current_before
    if (accepted) {
      nodes <- c(nodes[-j],children); attempted <- c(attempted[-j],rep(FALSE,length(children)))
      control$groups <- c(control$groups[-j],rep(group,length(children)))
    } else attempted[j] <- TRUE
    current_B <- log_sum(envelope_log_masses(nodes))
    gain <- current_before-current_B
    control <- update_refinement_from_nodes(control,group,gain,nodes,attempted,current_B)
  }
  out <- update_shells(model,nodes)
  out$engine <- "directional_shell_rejection"
  print_refinement("Directional",before,out$root_log_bound,start,verbose)
  invisible(out)
}

## 6) Local and ACG refinement ==========================================

local_remainder_coefficients <- function(eta,anchor,lower,upper,bins) {
  # Lipschitz padding also covers points between grid nodes.
  n <- length(eta); width <- (upper-lower)/bins
  tt <- outer(width,(seq_len(bins)-.5))+lower
  p <- plogis(eta)
  F <- function(t) p*t-.5*(softplus(eta+t)-softplus(eta-t))
  fa <- F(anchor)
  fpa <- p-.5*(plogis(eta+anchor)+plogis(eta-anchor))
  fppa <- .5*(plogis(eta-anchor)*plogis(-eta+anchor)-
                 plogis(eta+anchor)*plogis(-eta-anchor))
  hh <- tt-anchor; small <- abs(hh)<1e-3; den <- hh^2; den[small] <- 1
  DD <- 2*(F(tt)-fa-fpa*hh)/den
  er <- 512*.Machine$double.eps*(1+abs(eta)+abs(tt)+abs(anchor))/den
  second <- matrix(fppa,n,bins)
  DD[small] <- second[small]
  er[small] <- abs(hh[small])/(18*sqrt(3))

  xt <- tt-rowMeans(tt); yy <- DD-rowMeans(DD)
  ss <- rowSums(xt^2); co <- numeric(n); nz <- ss>0
  co[nz] <- -3*rowSums(xt*yy)[nz]/ss[nz]

  RR <- DD+co*(tt+2*anchor)/3
  pad <- width*(1/(18*sqrt(3))+abs(co)/3)/2
  low <- apply(RR-er,1L,min)-pad
  high <- apply(RR+er,1L,max)+pad
  if (any(!is.finite(c(co,low,high,fa,fpa)))) stop("Local remainder is not numerically resolvable.")
  list(c3=co,lower=low,upper=high,value=fa+co*anchor^3/6,
       derivative=fpa+co*anchor^2/2)
}

local_envelope_candidate <- function(model,node,bins=128L,geometry_only=FALSE) {
  gg <- model$envelope_geometry; Q <- node$Q; A <- gg$B %*% Q
  d <- model$d; k <- d-1L; tc <- (node$left+node$right)/2
  hw <- (node$right-node$left)/2; rh <- node$rho_hi; rl <- node$rho_lo
  alpha <- A[,1L]; beta <- sqrt(rowSums(A[,-1L,drop=FALSE]^2))
  support <- function(sg) {
    tt <- pmax(node$left,pmin(node$right,node$hi*sg*alpha/pmax(gg$norm_b,1e-300)))
    sg*alpha*tt+beta*sqrt(pmax(0,node$hi^2-tt^2))
  }
  sp <- 128*.Machine$double.eps*(1+node$hi*gg$norm_b)
  lower <- -support(-1)-sp; upper <- support(1)+sp
  co <- local_remainder_coefficients(gg$eta,alpha*tc,lower,upper,bins)
  T <- cubic_tensor(A,co$c3)
  D <- array(T[-1L,-1L,-1L,drop=FALSE],c(k,k,k))
  harmonic <- harmonic_decomposition(D,sum(abs(T)))
  vv <- harmonic$v; hn <- harmonic$norm
  aa <- T[1L,1L,1L]; bv <- as.numeric(T[1L,1L,-1L]); C <- matrix(T[1L,-1L,-1L],k,k)
  ss <- as.numeric(crossprod(Q,gg$s)+crossprod(A,co$derivative))
  c0 <- sum(co$value-co$derivative*alpha*tc)
  gp <- ss[-1L]-.5*tc^2*bv-vv*(rh^2+rl^2)/(4*(k+2))
  tr <- max(abs(node$left^2-tc^2),abs(node$right^2-tc^2),
            if (node$left<=0 && node$right>=0) tc^2 else 0)
  harmonic <- hn*rh^3/6
  E0 <- .5*tr*sqrt(sum(bv^2))*rh+
    sqrt(sum(vv^2))*(rh^2-rl^2)*rh/(4*(k+2))+harmonic+abs(aa)*hw^3/6

  E0 <- E0+model$bound_padding*(1+sum(abs(T))+sum(abs(gg$eta))+
         node$hi*sum(gg$norm_b)+sum(abs(ss))+abs(c0))*(1+node$hi^3+abs(E0))+model$prior_bound(node$hi)
  matrices <- vector("list",2L)
  ap0 <- ar0 <- cross <- numeric(2L)
  for (i in 1:2) {
    sg <- c(1,-1)[i]; curvature <- if (sg==1) co$upper else -co$lower
    M <- crossprod(A,A*curvature)
    M <- (M+t(M))/2+diag(128*.Machine$double.eps*(1+sum(abs(M))),d)
    matrices[[i]] <- M
    cross[i] <- hw*sqrt(sum(M[1L,-1L]^2))*rh
    ap0[i] <- -sg*aa*tc+M[1L,1L]
    ar0[i] <- max(largest_eigenvalue_bound(matrix(M[-1L,-1L],k,k)-sg*node$left*C),
                  largest_eigenvalue_bound(matrix(M[-1L,-1L],k,k)-sg*node$right*C))
  }
  geometry <- list(a=aa,b=bv,C=C,v=vv,hn=hn,s=ss,c0=c0,M=matrices,
    tc=tc,gp=gp,E0=E0,cross=cross,ap0=ap0,ar0=ar0)
  if (geometry_only) return(geometry)
  local_envelope_from_geometry(node,geometry)
}

acg_branch_mass <- function(br,cell) {
  k <- cell$d-1L
  br$vp <- 1/(1-br$ap); vr <- 1/(1-br$ar); br$mp <- br$vp*br$bp
  br$normal <- normal_interval_mass((cell$left-br$mp)/sqrt(br$vp),
                                 (cell$right-br$mp)/sqrt(br$vp))
  br$radial <- noncentral_shell(cell$rho_lo,cell$rho_hi,k,0,vr)
  br$logweight <- br$constant+.5*log(br$vp)+.5*br$vp*br$bp^2+
    .5*k*log(vr)+br$normal$logmass+br$radial$lm-.5*sum(log(br$omega))
  if (!is.finite(br$logweight)) stop("ACG branch mass is not numerically resolvable.")
  br
}

acg_branch <- function(node,cc,gt,g,ap0,H,xi) {
  # Bingham-ACG bound of Kent, Ganeiber and Mardia.
  k <- length(g); tc <- (node$left+node$right)/2; hw <- (node$right-node$left)/2

  bc <- if (xi==0) .5 else tanh(xi)/(2*xi)
  ac <- log_cosh(xi)-bc*xi^2
  S <- H+2*bc*tcrossprod(g); S <- (S+t(S))/2
  es <- eigen(S,symmetric=TRUE)
  pad <- 512*.Machine$double.eps*(1+sum(abs(S)))
  lambda <- es$values[1L]+pad
  eigA <- .5*node$rho_lo^2*pmax(0,lambda-es$values)
  if (any(!is.finite(eigA))) return(NULL)
  b <- k
  if (max(eigA)>1e-12) {
    ff <- function(b) sum(1/(b+2*eigA))-1

    b <- tryCatch(uniroot(ff,c(1e-12,k),tol=1e-10,check.conv=TRUE)$root,
                  error=function(e) min(1,k))
  }
  omega <- 1+2*eigA/b
  logCb <- -(k-b)/2+.5*k*log(k/b)
  best <- NULL
  for (cap in c(0,.5,.9,.99)) {
    ap <- min(ap0,cap); ar <- min(lambda,cap)
    constant <- cc+.5*(ap0-ap)*hw^2+.5*(lambda-ar)*node$rho_hi^2-
      gt*tc+.5*ap*tc^2+ac+logCb-log(2)
    br <- list(kind="acg",constant=constant,bp=gt-ap*tc,gperp=g,ap=ap,ar=ar,
      V=es$vectors,omega=omega)
    x <- tryCatch(acg_branch_mass(br,node),error=function(e) NULL)
    if (!is.null(x) && (is.null(best) || x$logweight<best$logweight)) best <- x
  }
  best
}

angular_envelope_candidate <- function(model,node,bins=128L) {
  gg <- node$local_geometry
  if (is.null(gg)) {
    ca <- local_envelope_candidate(model,node,bins,geometry_only=TRUE)
    if (is.null(ca)) return(NULL)
    gg <- ca
  }
  tc <- gg$tc; hw <- (node$right-node$left)/2; k <- node$d-1L
  ec <- eigen(gg$C,symmetric=TRUE)
  absC <- tcrossprod(sweep(ec$vectors,2L,abs(ec$values),"*"),ec$vectors)
  absC <- (absC+t(absC))/2
  branches <- vector("list",2L)
  for (i in 1:2) {
    sg <- c(1,-1)[i]

    H <- matrix(gg$M[[i]][-1L,-1L],k,k)-sg*tc*gg$C+hw*absC
    H <- (H+t(H))/2+diag(512*.Machine$double.eps*(1+sum(abs(H))+hw*sum(abs(gg$C))),k)
    cc <- sg*(gg$c0+gg$s[1L]*tc-gg$a*tc^3/6)+gg$E0+gg$cross[i]
    gt <- sg*(gg$s[1L]-.5*gg$a*tc^2); g <- sg*gg$gp
    base <- sqrt(sum(g*g))*node$rho_hi/sqrt(k)
    xis <- if (base==0) 0 else unique(c(0,.5*base,base,2*base,4*base,base*sqrt(k)))
    best <- NULL
    for (xi in xis) {
      x <- tryCatch(acg_branch(node,cc,gt,g,gg$ap0[i],H,xi),error=function(e) NULL)
      if (!is.null(x) && (is.null(best) || x$logweight<best$logweight)) best <- x
    }
    if (is.null(best)) return(NULL)
    branches[[i]] <- best
  }
  list(branches=branches,logweight=log_add(branches[[1L]]$logweight,branches[[2L]]$logweight),
    geometry=gg)
}

local_envelope_from_geometry <- function(node,gg) {
  tc <- gg$tc
  branches <- vector("list",2L)
  for (i in 1:2) {
    sg <- c(1,-1)[i]
    gt <- sg*(gg$s[1L]-.5*gg$a*tc^2)
    cc <- sg*(gg$c0+gg$s[1L]*tc-gg$a*tc^3/6)+gg$E0+gg$cross[i]
    best <- best_gaussian_branch(node,cc,gt,sg*gg$gp,gg$ap0[i],gg$ar0[i])
    if (is.null(best)) return(NULL)
    branches[[i]] <- best
  }
  list(branches=branches,logweight=log_add(branches[[1L]]$logweight,branches[[2L]]$logweight),
    geometry=gg)
}

hybrid_cache_key <- function(node,bins) list(bins=as.integer(bins),
  limits=c(node$lo,node$hi,node$left,node$right,node$rho_lo,node$rho_hi))

hybrid_envelope_candidate <- function(model,node,bins=128L) {
  reuse <- !is.null(node$local_geometry) &&
    identical(node$local_geometry_key,hybrid_cache_key(node,bins))
  local <- tryCatch(if (reuse) local_envelope_from_geometry(node,node$local_geometry) else
    local_envelope_candidate(model,node,bins),error=function(e) NULL)
  angular_node <- node
  if (!is.null(local)) angular_node$local_geometry <- local$geometry else if (!reuse)
    angular_node$local_geometry <- NULL
  angular <- tryCatch(angular_envelope_candidate(model,angular_node,bins),error=function(e) NULL)
  choices <- list(local,angular)
  lw <- vapply(choices,function(x) if (is.null(x)) Inf else x$logweight,numeric(1))
  good <- which(is.finite(lw))
  if (length(good)) choices[[good[which.min(lw[good])]]] else NULL
}

refine_skew_cosh_directional_angular <- function(model,max_cells=Inf,n_radial_steps=16L,
                                                angular_bins=128L,local_bins=128L,verbose=TRUE,
                                                min_log_bound_improvement=0,stagnation_patience=1L) {
  refine_skew_cosh_directional_local(model,max_cells,n_radial_steps,angular_bins,
    local_bins,verbose,bound_method="angular",
    min_log_bound_improvement=min_log_bound_improvement,stagnation_patience=stagnation_patience)
}

refine_skew_cosh_directional_local <- function(model,max_cells=Inf,n_radial_steps=16L,
                                                angular_bins=128L,local_bins=128L,verbose=TRUE,
                                                bound_method=c("local","angular"),
                                                min_log_bound_improvement=0,stagnation_patience=1L) {
  bound_method <- match.arg(bound_method); angular_method <- bound_method=="angular"
  label <- if (angular_method) "ACG" else "Local"
  start <- proc.time()[["elapsed"]]; before <- model$root_log_bound
  nodes <- model$shell_nodes; attempted <- if (angular_method)
    vapply(nodes,function(x) x$kind=="directional" &&
      identical(x$hybrid_checked,hybrid_cache_key(x,local_bins)),logical(1)) else rep(FALSE,length(nodes))
  control <- make_refinement_control(nodes,min_log_bound_improvement,stagnation_patience)
  used <- radial_used <- 0L
  refresh <- function(x) {
    if (angular_method) {
      ca <- hybrid_envelope_candidate(model,x,local_bins)
      x$hybrid_checked <- hybrid_cache_key(x,local_bins)
    } else ca <- tryCatch(local_envelope_candidate(model,x,local_bins),error=function(e) NULL)
    if (!is.null(ca) && ca$logweight<x$logweight-1e-10*(1+abs(x$logweight))) {
      x$branches <- ca$branches; x$logweight <- ca$logweight
      x$local_geometry <- ca$geometry
      x$local_geometry_key <- hybrid_cache_key(x,local_bins)
    }
    x
  }
  repeat {
    budget_eligible <- vapply(nodes,function(x) {
      if (!is.finite(x$hi)) return(FALSE)
      if (x$kind=="directional") return(used<max_cells)
      radial_used<n_radial_steps && used+angular_bins<=max_cells
    },logical(1))
    eligible <- which(!attempted & !control$stopped[control$groups] & budget_eligible)
    if (!length(eligible)) break
    j <- eligible[which.max(envelope_log_masses(nodes[eligible]))]
    group <- control$groups[j]
    parent <- nodes[[j]]
    current_before <- log_sum(envelope_log_masses(nodes))
    if (parent$kind=="directional") {
      child <- refresh(parent); used <- used+1L
      cm <- child$logweight
      candidate_B <- log_sum(c(envelope_log_masses(nodes[-j]),cm))
      accepted <- child$logweight<parent$logweight && candidate_B<=current_before

      if (accepted) nodes[[j]] <- child else if (angular_method) {
        parent$hybrid_checked <- child$hybrid_checked
        nodes[[j]] <- parent
      }
      attempted[j] <- TRUE
    } else {
      radial_used <- radial_used+1L
      dr <- skew_cosh_direction(model,parent$hi)
      cg <- directional_geometry(model,parent$hi,direction_basis(dr$direction),integrated_cubic_coefficients,TRUE)
      grid <- seq(-parent$hi,parent$hi,length.out=angular_bins+1L)
      children <- lapply(seq_len(angular_bins),function(i) {
        x <- make_directional_cell(parent,grid[i],grid[i+1L],cg)
        if (is.null(x)) return(NULL)
        refresh(x)
      })
      used <- used+angular_bins
      cm <- if (any(vapply(children,is.null,logical(1)))) Inf else
        log_sum(envelope_log_masses(children))
      candidate_B <- if (is.finite(cm)) log_sum(c(envelope_log_masses(nodes[-j]),cm)) else Inf
      accepted <- is.finite(cm) && cm<parent$logweight-1e-8*(1+abs(parent$logweight)) &&
        candidate_B<=current_before
      if (accepted) {
        nodes <- c(nodes[-j],children); attempted <- c(attempted[-j],rep(TRUE,length(children)))
        control$groups <- c(control$groups[-j],rep(group,length(children)))
      } else attempted[j] <- TRUE
    }
    current <- log_sum(envelope_log_masses(nodes))
    gain <- current_before-current
    control <- update_refinement_from_nodes(control,group,gain,nodes,attempted,current)
  }
  out <- update_shells(model,nodes)
  print_refinement(label,before,out$root_log_bound,start,verbose)
  invisible(out)
}
