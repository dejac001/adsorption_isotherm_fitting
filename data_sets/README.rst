

CO2 and N2 on BEA
-----------------
Experiment from :cite:`Pham2014`

.. csv-table:: N2 adsorption on BEA, scanned from images in paper
    :widths: 20 20 20 20
    :file: N2_BEA.csv

.. csv-table:: CO2 adsorption on BEA, scanned from images in paper
    :widths: 20 20 20 20
    :file: CO2_BEA.csv


H2S and CH4 on MFI
------------------

Molecular simulation from :cite:`Shah2015`,


.. csv-table:: H2S adsorption on MFI, taken from tables in SI
    :widths: 20 20 20 20 20
    :file: H2S_MFI.csv

.. csv-table:: CH4 adsorption on MFI, taken from tables in SI
    :widths: 20 20 20 20 20
    :file: H2S_MFI.csv

.. csv-table:: H2S/CH4 Binary adsorption on MFI, taken from tables in SI
    :widths: 20 20 20 20 20 20 20 20 20 20 20 20
    :file: CH4_H2S_MFI_binary.csv

Fugacity coefficients and fugacities
************************************

Calculated from :code`RealGas` python package using virial equation of state.
The critical temperatures and densities were obtained from the TraPPE website.
The accentric factors and critical compressibilities were obtained from DIPPR :cite:`DIPPR`.
The critical pressure was calculated from all the other critical properties above.
The :math:`k_ij` parameter was set to 0.


.. csv-table:: H2S/CH4 Binary adsorption on MFI, taken from tables in SI and added fugacities
    :widths: 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20
    :file: CH4_H2S_MFI_binary_with_fugacity.csv


.. bibliography:: ../source/references.bib
