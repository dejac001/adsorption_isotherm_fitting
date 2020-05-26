Unary Isotherms
===============

Langmuir
--------
The temperature-dependent unary Langmuir isotherm is expressed as

.. math::
    q_i = \frac{q_{\text{m},i}k_if_i}{1 + k_i f_i}
    :label: eq_lang_unary

where :math:`f_i`, is the fugacity of component *i*,
can be calculated assuming ideal gas

.. math::
    f_i^\text{IG} = y_i P

or, using the RealGas_ package to calculate :math:`\phi_i`
from :math:`y_i,P,T` data,

.. math::
    f_i = \phi_i y_i P

An Arrhenius relationship for :math:`k_i` is assumed as

.. math::
    k_i = k_{i,\infty}\exp\left(\frac{-\Delta H_i}{RT}\right)

Introducing the dimensionless parameters

.. math::
    \theta_i = \frac{q_i}{q_\text{ref}}
    :label: eq_theta

.. math::
    f_i^\star = \frac{f_i}{f_\text{ref}}
    :label: eq_fis_unary

.. math::
    T^\star = \frac{T}{T_\text{ref}}
    :label: eq_Ts

The variables to be fit in dimensionless form are

.. math::
    H_i^\star = \frac{\Delta H_i}{R T_\text{ref}}
    :label: H_i_star

.. math::
    q_{\text{m},i}^\star = \frac{q_{\text{m},i}}{q_\text{ref}}
    :label: q_mi_star

.. math::
    A_i = \ln\left(k_{i,\infty} f_\text{ref}\right)
    :label: A_i

So that Equation :eq:`eq_lang_unary` becomes

.. math::
    \theta_i = \frac{q_{\text{m},i}^\star\exp\left(A_i - \frac{H_i^\star}{T^\star}\right)f_i^\star}{1 + \exp\left(A_i - \frac{H_i^\star}{T^\star}\right)f_i^\star}
    :label: eq_lang_unary_dimensionless


Modules
-------

.. automodule:: isotherm_models.unaryisotherm
    :members:



.. _RealGas: https://github.com/dejac001/RealGas
