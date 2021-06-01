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
Could perhaps use this instead of `calc_sgv()` ?

More generally, I'm not sure why the computation is using SGV. I see that `r2beta.lmerMod` is internally converting SigHat (the full covariance matrix) to its SGV before calling `cmp_R2`. I'm not sure why this is done?
`cmp_R2` will take SigHat as a scalar or as a matrix ... the latter seems more general and closer in spirit to what the paper suggests ... ?

(2) The block size is also used in ddf computations. Presumably if I have the more general situation I could obtain the Kenward-Roger ddf for the model.

These are the results from trying all three levels of clustering all three methods of collapsing Sigma-hat ("blk_sgv" is the current method; det_sgv uses the determinant calculation above; none leaves Sigma-hat as a matrix. If I were to implement the K-R ddf approach as well that would give me another factorial level.

       clust.id collapse_sighat      F v1    v2    ncp    Rsq
cc        biome         blk_sgv  96.63 11 619.0 1063.0 0.6320
cc1       biome         det_sgv 151.96 11 619.0 1671.5 0.7298
cc2       biome            none  17.52 11 619.0  192.7 0.2374
cc3 flor_realms         blk_sgv 140.54 11 619.0 1545.9 0.7141
cc4 flor_realms         det_sgv 151.96 11 619.0 1671.5 0.7298
cc5 flor_realms            none  17.52 11 619.0  192.7 0.2374
cc6    biome_FR         blk_sgv  37.19 11 540.6  409.1 0.4308
cc7    biome_FR         det_sgv 151.96 11 540.6 1671.5 0.7556
cc8    biome_FR            none  17.52 11 540.6  192.7 0.2628


  Naively, I would think that using K-R ddf would be best (since it is what is used in the original paper, and doesn't depend on pretending Sigma-hat is block-diagonal), and that leaving Sigma-hat as a matrix is also best (although it disappoints me a bit that this leads to considerably lower R^2 values ...)

## issues with gamm4 models

For my own information

- `gamm4` doesn't return objects with class `gamm4` (although it should); it returns a list containing a `mer` object and a(n) `mgcv` object
- `cmp_R2` (i.e., allowing for arbitrary contrasts) has a different setup/interface than `r2glmm(..., partial=TRUE)`.  Seems hard to use K-R ddf in a `cmp_R2` context ? Maybe I should just give up and do model comparison?
- model updating is tricky for `gamm4` objects:
    - `update()` and `drop1()` don't work well
