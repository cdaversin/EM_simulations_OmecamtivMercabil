# UPV Matlab scripts
MATLAB files to run cellular electro-mechanical simulations with Omecamtiv Mecarbil

## [Main codes](#main)
  1 - `fitting_3ks_om.m` --> script to obtain parameters of the OM model \
    2 - `basic_simulation.m` --> script to run a single simulation to obtain electro-mechanical signals and compute biomarkers. Selection of BCL, OM concentration and failing or non-failing conditions.

## Subfunctions
  `ObjFuncT.m` --> Objective function
    `model_Tor_Land_OM.m` --> healthy myocyte model
    `model_Tor_Land_HF_OM.m` --> myocyte model with heart failure conditions
    `hill_eq.m` --> Hill equation
## Biomarker subfunctions
  `APD.m` --> action potential metrics
    `Tx` --> active tension metrics

## Data
  `drug_points_shen21.mat` --> Experimental data points extracted from Shen et al. 2021 (DOI: 10.1161/JAHA.121.020860), Figure C
    `om_paramsTest1.mat` --> OM parameters obtained with ([main codes](#main) - 1)
