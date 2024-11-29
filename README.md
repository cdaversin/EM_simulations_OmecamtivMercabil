
Software, source code and data used in "Insights from Electromechanical Simulations to Assess Omecamtiv Mecarbil Efficacy in Heart Failure", by Maria Teresa Mora, Ilse van Herck, Cécile Daversin-Catty, Henrik Finsberg, Jordi Llopis-Lorente,
Javier Saiz, Hermenegild Arevalo, Samuel Wall, and Beatriz Trenor.

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


# Simula Cardiac Electro-Mechanics Solver

`simcardems` is a FEniCS-based cardiac electro-mechanics solver and is developed as a part of the [SimCardio Test project](https://www.simcardiotest.eu/wordpress/). The solver depends on [`pulse`](https://github.com/ComputationalPhysiology/pulse) and [`cbcbeat`](https://github.com/ComputationalPhysiology/cbcbeat).


## Installation

See [Installation instructions](https://computationalphysiology.github.io/simcardems/install.html)

## Getting started

See [the demos](https://computationalphysiology.github.io/simcardems/simple_demo.html)

## Documentation

Documentation is hosted at http://computationalphysiology.github.io/simcardems.

## Automated test

Tests are provided in the folder [tests](https://github.com/ComputationalPhysiology/simcardems/tree/main/tests). You can run the tests with pytest

```
python3 -m pytest tests -vv
```

## Contributing
See [the contributing section](https://computationalphysiology.github.io/simcardems/CONTRIBUTING.html)

## Citing
If you use `simcardems` in your research, please cite it as follows
```
@article{Finsberg2023,
  doi = {10.21105/joss.04753},
  url = {https://doi.org/10.21105/joss.04753},
  year = {2023},
  publisher = {The Open Journal},
  volume = {8},
  number = {81},
  pages = {4753},
  author = {Henrik Nicolay Topnes Finsberg and Ilsbeth Gerarda Maria van Herck and Cécile Daversin-Catty and Hermenegild Arevalo and Samuel Wall},
  title = {simcardems: A FEniCS-based cardiac electro-mechanics solver},
  journal = {Journal of Open Source Software}
 }
```

## Known issues

- Issue with h5py, see https://github.com/ComputationalPhysiology/pulse#known-issues


## Authors
- Henrik Finsberg (henriknf@simula.no)
- Ilsbeth van Herck (ilse@simula.no)
- Cécile Daversin-Catty (cecile@simula.no)
