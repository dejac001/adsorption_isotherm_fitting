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
