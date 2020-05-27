Installation
============

Docker
------

.. code-block:: bash

    docker pull dejac001/isotherm-fitting-users:0.0.4
    docker run -ti -v $PWD:/home/pyomo/shared/ dejac001/isotherm-fitting-users:0.0.4 # run interactively inside container (ubuntu-based)

Singularity
-----------

.. code-block:: bash

    module load singularity
    singularity pull docker://dejac001/isotherm-fitting-users:0.0.4
    mv isotherm-fitting-users_0.0.4.sif /path/to/shared/directory/isotherm-fitting-users_0.0.4.sif
    singularity exec -B $PWD:/home/pyomo/shared /path/to/shared/directory/isotherm-fitting-users_0.0.4.sif python3 path/to/input/file.py

Scipy Only
----------

.. code-block:: bash

    pip3 install Pyomo chem-util matplotlib pandas numpy realgas>=1.0.2
    python3 -m pip3 install https://github.com/dejac001/adsorption_isotherm_fitting/archive/v0.0.4.tar.gz