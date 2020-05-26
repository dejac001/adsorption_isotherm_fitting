Installation
============

Docker
------

.. code-block:: bash

    docker pull dejac001/isotherm-fitting-users:0.0.3
    docker run -ti -v $PWD:/home/pyomo/shared/ dejac001/isotherm-fitting-users:0.0.3 # run interactively inside container (ubuntu-based)

Singularity
-----------

.. code-block:: bash

    module load singularity
    singularity pull docker://dejac001/isotherm-fitting-users:0.0.3
    mv isotherm-fitting-users_0.0.3.sif /path/to/shared/directory/isotherm-fitting-users_0.0.3.sif
    singularity exec -B $PWD:/home/pyomo/shared /path/to/shared/directory/isotherm-fitting-users_0.0.3.sif python3 path/to/input/file.py
