# adsorption_isotherm_fitting
[![DOI](https://zenodo.org/badge/255999657.svg)](https://zenodo.org/badge/latestdoi/255999657)

Fit temperature-dependent isotherms to equilibrium data.

# Getting Started

[Documentation](https://adsorption-isotherm-fitting.readthedocs.io/en/latest/)


Installation
============

Docker
------
```bash
docker pull dejac001/isotherm-fitting
docker run -ti -v $PWD:/home/pyomo/shared/ # run interactively inside container (ubuntu-based)
```
Singularity
-----------
```bash
module load singularity
singularity pull docker://dejac001/isotherm-fitting
mv isotherm-fitting_latest.sif /path/to/shared/directory/isotherm-fitting_latest.sif
singularity exec -B $PWD:/home/pyomo/shared /path/to/shared/directory/isotherm-fitting_latest.sif python3
```

Examples
========

H2S/CH4 on MFI Zeolite
----------------------

Fit temperature-dependent Langmuir unary and binary isotherms for H2S/CH4 mixture isotherms.
Data for fitting taken from [Shah et al, 2015](https://doi.org/10.1021/acs.langmuir.5b03015).

Example script [here](examples/h2s_ch4.py)

<p align="center">
    <img
        src="examples/h2s_ch4_example.png"
        width="640"
    />
</p>

# See Also
[PyGAPS](https://github.com/pauliacomi/pyGAPS)
[PyIAST](https://github.com/CorySimon/pyIAST)


