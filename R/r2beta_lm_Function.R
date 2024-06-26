
#' @export

r2beta.lm <- function(model, partial=TRUE, method='sgv',
                      data = NULL, ...) {

  check_dots(...)

  beta = stats::coef(model)
  p = length(beta)
  X = stats::model.matrix(model)
  SigHat = summary(model)$sigma^2

  # C matrix defines the Wald Test for Fixed Effects
  C = list(); nms = c('Model', names(beta)[-1])

    # Define the model Wald statistic for all fixed effects
    C[['Model']] = cbind(rep(0, p-1),diag(p-1))

    # For partial R2 statistics:
    if (partial) {

      ## copied from lmer (DNR!)
      labs = attr(stats::terms(model), 'term.labels')
      asgn = attr(X, 'assign')
      nmrs = 1:length(asgn)
      assign = split(nmrs, asgn)
      nTerms = length(assign)
      nms = c('Model', labs)
      names(assign) = c('(Intercept)',labs)

      # add the partial contrast matrices to C
      for(i in 2:(nTerms)) {
        C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = assign[[i]])
      }
    }

    # Compute the specified R2
    r2=lapply(C, FUN=cmp_R2, x=X, SigHat=SigHat,
              beta=beta, method='lm')

    # initialize a dataframe to hold results
    R2 = data.frame(Effect = names(r2))

    # place results in the dataframe
    for(i in names(r2[[1]])){
      R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
    }

    R2 = within(R2, {
      lower.CL = stats::qbeta(0.025, R2$v1/2, R2$v2/2, R2$ncp)
      upper.CL = stats::qbeta(0.975, R2$v1/2, R2$v2/2, R2$ncp)
    } )
    R2 = R2[order(-R2$Rsq),]

    class(R2) <- c('R2', 'data.frame')

    return(R2)

}
