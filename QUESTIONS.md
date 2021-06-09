The documentation says that Kenward Roger is only available for `lme` models ([here](https://github.com/bcjaeger/r2glmm/blob/master/man/r2beta.Rd#L19-L20)); this is a simple typo but led me down the wrong track for quite a while!

## standardized generalized variance

- Where is the standardized generalized variance approach documented? It's not easily located in Jaeger 2016.
- Handling of crossed models (cf [this issue](https://github.com/bcjaeger/r2glmm/issues/9#issuecomment-600783392))

My models have multiple random effect terms/grouping variables (roughly speaking, `flor_realms` in the example below is the largest-scale/smallest number of levels (n=6, mean 103.3 obs/block); `biomes` is smaller-scale but crossed with `flor_realms` (14, 44.3); and `biome_FR` (55, 11.2) represents the interaction between these two sampling levels.  In the paper, I see (below eq 10)

> \hat\Sigma is the estimated covariance matrix for the pseudo otucomes having a block diagonal structure with subject-specific estimated covariance matrices \hat\Sigma_i = Z_i \hat\Sigma_b Z'_i + \hat W_i{-1}, i=1, ..., n constituting the blocks

However, my Sigma is not block-diagonal because it combines information from three different sampling levels ...  `r2beta.lmerMod` appears to choose the *last* grouping variable automatically to do the block size/samples per block calculations (which if I recall correctly is always the grouping variables with the *fewest* levels given lme4's internal rearrangements).	`

The block size info seems to be used in two places: 

(1) in `calc_sgv()`. I see that `?calc_sgv` says "calculates the standardized generalized variance for a block diagonal matrix", but I don't see SGV referred to anywhere in Jaeger et al 2017 (maybe there's a later reference?). As far as I can tell the definition of SGV is more generally

```{r sgv}
exp(determinant(SigHat,log=TRUE)$modulus/nrow(SigHat))
```

Could perhaps use this expression instead of `calc_sgv()` ?

More generally, I'm not sure why the computation is using SGV. I see that `r2beta.lmerMod` is internally converting SigHat (the full covariance matrix) to its SGV before calling `cmp_R2`. I'm not sure why ... ?

Instead of what's done in `r2beta`, `cmp_R2` will (apparently) take SigHat either as a scalar or as a matrix ... the latter seems more general and closer in spirit to what the paper suggests ... ?

`cmp_R2` (i.e., allowing for arbitrary contrasts) has a different setup/interface than `r2glmm(..., partial=TRUE)`.  Seems hard to use K-R ddf in a `cmp_R2` context ?

## issues with gamm4 models

`gamm4` models are pretty terrible to work with for a variety of reasons:

- `call` object is mangled to uselessness
- the formulas and terms are not very useful; `drop1` doesn't really work
- I've resorted to setting up an interface where the user has to specify the partial terms manually, although I have a hack to guess it
- hard to know what the null model should be (include smooth terms or not?)
