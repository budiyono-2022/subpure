# subpure

Method-purified subscore added value for bifactor models with wording factors.

Questionnaires are often built from several theoretical facets and
balance positively and negatively worded items. Whether the facets deserve
separate scores is usually judged with PRMSE/value-added ratios or bifactor
indices that ignore the wording factor. `subpure` fits a bifactor model with a
wording factor, defines each facet's target true score **without** wording
variance, and reports how well the facet subscore, the total score and an
augmented score recover that target — with Monte Carlo confidence intervals, a
one-sided test of added value, between-group contrasts, and the method-blind
values for comparison.

```r
# install.packages("remotes")
remotes::install_github("budiyono-2022/subpure")
library(subpure)
data(rmbangka)
items  <- names(rmbangka)[7:42]
facets <- split(items, substr(items, 1, 1))
method <- c("A2","A3","A6","A7","A10","B12","B14","B16","B19",
            "C21","C22","C25","C26","C29","C30","D32","D34","D36")
pure_fit(rmbangka, facets, method)            # five competing structures
pure_value(rmbangka, facets, method, R = 1000) # purified added value
pure_value(rmbangka, facets, method, group = "school_type")
```

Main functions: `pure_syntax()`, `pure_fit()`, `pure_value()`,
`simulate_pure()`, `pure_population()`. The bundled `rmbangka` data
(N = 4,197) are from Pratama et al. (2026, *Journal of Open Psychology Data*,
doi:10.5334/jopd.163), released under CC BY 4.0.

`inst/simulation/` contains the simulation script and `inst/analysis/` the
script reproducing the empirical example in the accompanying article.
