# North American Breeding Bird Survey status and trend estimates to inform a wide-range of conservation needs, using a flexible Bayesian hierarchical generalized additive model

The status and trend estimates derived from the North American Breeding Bird Survey (BBS), are critical sources of information for bird conservation. However, the estimates are partly dependent on the statistical model used. Therefore, multiple models are useful because not all of the varied uses of these estimates (e.g. inferences about long-term change, annual fluctuations, population cycles, recovery of once declining populations) are supported equally well by a single statistical model. Here we describe Bayesian hierarchical generalized additive models (GAM) for the BBS, which share information on the pattern of population change across a species' range. We demonstrate the models and their benefits using data from a selection of species; and we run a full cross-validation of the GAMs against two other models to compare predictive fit. The GAMs have better predictive fit than the standard model for all species studied here, and comparable predictive fit to an alternative first difference model. In addition, one version of the GAM described here (GAMYE) estimates a population trajectory that can be decomposed into a smooth component and the annual fluctuations around that smooth. This decomposition allows trend estimates based only on the smooth component, which are more stable between years and are therefore particularly useful for trend-based status assessments, such as those by the International Union for the Conservation of Nature. It also allows for the easy customization of the model to incorporate covariates that influence the smooth component separately from those that influence annual fluctuations (e.g., climate cycles vs annual precipitation). For these reasons and more, this GAMYE model is a particularly useful model for the BBS-based status and trend estimates.

This repository contains all of the code, models (JAGS), and facilities to download the raw BBS data.

This work was conducted primarily on a Windows computer that had 148 GB of RAM and 64 processor cores, so users should adjust the calls to the foreach and doParallel packages to suite their local set-up.

There are 3 main scripts that provided the analysis, 1 that produces the figures for the paper, and 1 that compares some alternative priors for the hyperparameters of the GAM.

1.  script_kfold.R applies the four statistical models used to the BBS data for a selection of species, and then conducts a 15-fold crossvalidation for each model by species combination.

2.  comparison_script.R applies an additional model to the crossvalidation results to account for the imbalanced sampling (in time and space) of the BBS.

3.  comparison_summary.R further summarizes the results of the previous two scripts.

4.  prior_comparisons_GAM.r applies the GAM using two alternative priors on the variance of the beta hyperparameters, demonstrating that the population inferences of the GAM are not dependent on this prior (added based on a suggestion by reviewers).

5.  publication_figures.R draws on the output from the previous 4 scripts to produce the final figures for the publication, including the supplements.

The paper will be published in Condor, and a preprint of the paper is available at: doi: <https://doi.org/10.1101/2020.03.26.010215>
