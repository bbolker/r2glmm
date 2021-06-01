if (
    require("gamm4") && require("splines") && require("mgcv")
) {

  library(r2glmm)
  ## example with estimated zero among-block variance
  set.seed(101)
  dat2 <- gamSim(1,n=200,dist="normal", verbose=FALSE)
  dat2$fac <- rep(factor(1:20),length.out=nrow(dat2))
  b3 <- gamm4(y~s(x0)+s(x1)+s(x2),family=gaussian,
              data=dat2,random=~(1|fac))
  class(b3) <- "gamm4"
  r2_3 <- r2beta(b3)
  b4 <- lm(y~ns(x0,3)+ns(x1,3) + ns(x2,8), data=dat2)
  r2_4 <- r2beta(b4)

  ## similar enough that I'm not horrified
  test_that("gamm4 and raw splines similar", {
    expect_equal(r2_3$Rsq, r2_4$Rsq, tolerance=0.2)
  })

  r2_3kr <- r2beta.gamm4(b3, method="kr",
                         formula=y~s(x0)+s(x1)+s(x2),
                         partial.terms=c("s(x0)","s(x1)","s(x2)"),
                         data=dat2, random=~(1|fac))

  test_that("KR on gamm4", {
    expect_equal(r2_3kr$Rsq,
        c(0.662728341545902, 0.523745379980698,
          0.0400949752570711, 0.0382369181231837),
        tolerance=1e-4)

    ## automatic split of formula to get partial terms
    expect_warning(r2_3krB <- r2beta.gamm4(b3, method="kr",
                            formula=y~s(x0)+s(x1)+s(x2),
                            partial.terms=NULL,
                            data=dat2, random=~(1|fac)),
                   "attempting to reconstruct terms by splitting")
    expect_equal(r2_3krB, r2_3kr)

  })

  r2_3kr2 <- r2beta.gamm4(b3, method="kr",
                         formula=y~s(x0)+s(x1)+s(x2),
                         partial.terms=c("(s(x0)+s(x1))"),
                         data=dat2, random=~(1|fac))


}
