

library(r2glmm)
context("Computing r2beta")

library(nlme)
library(lme4)
data(Orthodont)

# Linear mixed models
mermod = lmer(distance ~ age*Sex + (1|Subject), data = Orthodont)
lmemod = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)

lme.r2 = r2beta(model = lmemod, method = 'sgv', partial = TRUE)
mer.r2 = r2beta(model = mermod, method = 'sgv', partial = TRUE)

r2.diff = lme.r2$Rsq - mer.r2$Rsq

test_that("r2beta works", {
  expect_equal(lme.r2$Rsq, mer.r2$Rsq, tolerance=0.01)
})


# Example from documentation of `r2beta`

dis = data.frame(discoveries)
dis$year = 1:nrow(dis)

glmod1 = glm(
  discoveries ~ year + I(year^2),
  family = 'poisson',
  data = dis
)

glmod2 = glm(
  as.formula("discoveries ~ year + I(year^2)"),
  family = 'poisson',
  data = dis
)

glmod3 = glm(
  "discoveries ~ year + I(year^2)",
  family = 'poisson',
  data = dis
)

test_that("r2beta allows programmatic formulas", {
  expect_equal(r2glmm::r2beta(glmod1), r2glmm::r2beta(glmod2))
  expect_equal(r2glmm::r2beta(glmod1), r2glmm::r2beta(glmod3))
})


