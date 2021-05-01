#' @rdname r2beta
#' @export r2beta.gamm4
#' @export
r2beta.gamm4 <- function(model, partial=TRUE, method='sgv',
                         data = NULL,
                         formula,
                         partial.terms=NULL,
                         ...) {

  model <- model$mer  ## do everything with the 'mer' component

  if (is.null(data)) {
    data <- try(eval.parent(getCall(model)$data))
    if (inherits(data,"try-error") || is.null(data)) {
      data = model.frame(model)
    }
  }

  if(is.null(data) && partial) {
    stop('Please specify the dataframe used to fit the model.')
  }

  # Get model matrices
  X = lme4::getME(model, 'X')
  n <- nrow(X)

  if(nrow(X) != nrow(data)) {
    stop('data do not have same number of rows as model data. Try removing rows with missing values for model variables')
  }

  # The Kenward Roger Approach
  if (toupper(method) == 'KR') {

    # Calculate F statistic using the approach of Kenward & Roger
    krfun <- function(comp_model, effect='Model') {
      mc = pbkrtest::KRmodcomp(model, comp_model)$stats
      ss = with(mc, ndf*Fstat/ddf)
      return(data.frame(Effect=effect, 'F' = mc$Fstat,
                    v1 = mc$ndf, v2 = mc$ddf,
                    pval = mc$p.value,
                    ncp = mc$ndf*mc$Fstat,
                    Rsq = ss / (1+ss)))
    }

    null_model <- gamm4(update(formula, . ~ 1), ..., data=data)
    R2 <- krfun(null_model$mer)

    # For partial R2 statistics:
    if(partial) {

      if (is.null(partial.terms)) {
        warning("attempting to reconstruct terms by splitting on '+'")
        partial.terms <- trimws(strsplit(deparse1(formula[[3]]),"\\+")[[1]])
        partial.terms <- setdiff(partial.terms, "1") ## remove intercept
      }

      R2_partial <- list()
      for (t0 in partial.terms) {
        reduced_form <- update(formula, as.formula(sprintf(". ~ . - %s",t0)))
        reduced_model <- gamm4(reduced_form, data=data, ...)
        R2_partial <- c(R2_partial, list(krfun(reduced_model$mer,t0)))
      }
      R2 <- do.call("rbind", c(list(R2), R2_partial))
    }

  } else {
    ## not KR
    # Get grouping information from the model (not needed for KR)
    clust.id = names(model@flist)[ length(model@flist) ]

    if(!clust.id %in% names(data)){

      ids = strsplit(clust.id, ':')[[1]]
      data[[clust.id]] = interaction(data[[ids[1]]], data[[ids[2]]])
      data = data = droplevels(data)

    }

    obsperclust = as.numeric(table(data[ , clust.id ]))
    mobs = mean(obsperclust)
    nclusts = length(obsperclust)

    if (toupper(method) %in% c('SGV', 'NSJ')) {

      beta = lme4::fixef(model)
      p <- length(beta)

      if(p==1) stop('Model must have at least one fixed effect')

      # Get random effects design matrix
      Z = lme4::getME(model, 'Z')

      # Get variance component estimates
      s2e  = lme4::getME(model, 'sigma')^2
      lam  = lme4::getME(model, 'Lambda')
      lamt = lme4::getME(model, 'Lambdat')

      G = s2e * (lam %*% lamt)

      # Compute estimated covariance matrix
      SigHat = Z %*% ( G%*%Matrix::t(Z) )

      # Add the residual component
      Matrix::diag(SigHat) = Matrix::diag(SigHat) + s2e/model@resp$weights

      if(toupper(method)=='NSJ'){

        # NSJ approach takes the mean of the trace of the model covariance
        SigHat = mean(Matrix::diag(SigHat))

      }

      if(toupper(method)=='SGV'){

        # SGV approach takes standardized
        # determinant of the model covariance
        SigHat=calc_sgv(nblocks=nclusts,
                        blksizes=obsperclust,
                        vmat=SigHat)

      }

      # C matrix defines the Wald Test for Fixed Effects
      C = list(); nms = c('Model', names(beta)[-1])

      # Define the model Wald statistic for all fixed effects

      C[['Model']] = cbind(rep(0, p-1),diag(p-1))

      # For partial R2 statistics:
      if (partial && p>1) {

        labs = attr(stats::terms(model), 'term.labels')
        if (identical(labs, "X")) {
          ## FIXME: can gamm4 models have non-X elements??
          ## HACK for gamm$mer objects; true 'assign' value has been
          ## lost due to replacement of orig formula by ~X ...
          assign <- seq(ncol(X))
          names(assign) <- gsub("^X","",colnames(X))
          nms <- c("Model", names(assign)[-1])
          nTerms <- ncol(X)
        } else {
          asgn = attr(X, 'assign')
          nmrs = 1:length(asgn)
          assign = split(nmrs, asgn)
          nTerms = length(assign)
          nms = c('Model', labs)
          names(assign) = c('(Intercept)',labs)
        }

        # add the partial contrast matrices to C
        for(i in 2:(nTerms)) {
          C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = assign[[i]])
        }

      }

      # Compute the specified R2
      r2=lapply(C, FUN=cmp_R2, x=X, SigHat=SigHat, beta=beta, method=method,
                obsperclust=obsperclust, nclusts=nclusts)

      # initialize a dataframe to hold results
      R2 = data.frame(Effect = names(r2))

      # place results in the dataframe
      for(i in names(r2[[1]])){
        R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
      }

    }
  }

  # Calculate confidence limits with the non-central beta quantile function.
  # Arrange the resulting dataframe from highest to lowest R^2 values.

  R2 = within(R2, {
    lower.CL = stats::qbeta(0.025, R2$v1/2, R2$v2/2, R2$ncp)
    upper.CL = stats::qbeta(0.975, R2$v1/2, R2$v2/2, R2$ncp)
  } )
  R2 = R2[order(-R2$Rsq),]

  class(R2) <- c('R2', 'data.frame')
  return(R2)

}
