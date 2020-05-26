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

.. _RealGas: https://github.com/dejac001/RealGas

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

.. autoclass:: isotherm_models.unaryisotherm.UnaryIsotherm
    :members:

.. autoclass:: isotherm_models.unaryisotherm.LangmuirUnary
    :members:


Binary Isotherms
================

Binary Langmuir
---------------

.. math::
    q_i = \frac{q_{\text{m},i}k_i\hat{f}_i}{1 + k_i \hat{f}_i + k_j \hat{f}_j}
    :label: eq_lang_binary

Arrhenius relationships are used for :math:`k_i`  and :math:`k_j`,
and dimensionless variables are used as
illustrated in :class:`isotherm_models.unaryisotherm.LangmuirUnary`

.. note::
    This isotherm is not equivalent to the conventional extended langmuir isotherm,
    because *both* :math:`k_i` and :math:`k_j` are fit simultaneously to binary data.

For completeness, the relationships are repeated for the binary case below

.. math::
    \begin{align}
        k_i &= k_{i,\infty}\exp\left(\frac{-\Delta H_i}{RT}\right)\\
        k_j &= k_{j,\infty}\exp\left(\frac{-\Delta H_j}{RT}\right)\\
    \end{align}

The dimensionless parameters :math:`\theta_i` and
:math:`T^\star` are calculated as the unary case,
as shown in Equations :eq:`eq_theta` and :eq:`eq_Ts`, respectively.
The other dimensionless parameters are

.. math::
    \hat{f}_i^\star = \frac{\hat{f}_i}{f_\text{ref}}
    :label: eq_fis_binary

.. math::
    \hat{f}_j^\star = \frac{\hat{f}_j}{f_\text{ref}}
    :label: eq_fjs_binary


The dimensionless variables to be fit
include
:math:`H_i^\star`,
:math:`q_{\text{m},i}^\star`,
:math:`A_i`,
:math:`H_j^\star`,
:math:`q_{\text{m},j}^\star`,
and
:math:`A_j`.
The former three (math:`H_i^\star`,
:math:`q_{\text{m},i}^\star`,
and :math:`A_i`) have the same
expression as the unary case,
as shown in Equations :eq:`H_i_star`, :eq:`q_mi_star`, and :eq:`A_i`, respectively.
The latter two are expressed as

.. math::
    H_j^\star = \frac{\Delta H_j}{R T_\text{ref}}
    :label: H_j_star

.. math::
    A_j = \ln\left(k_{j,\infty} f_\text{ref}\right)
    :label: A_j

So that Equation :eq:`eq_lang_binary` becomes

.. math::
    \theta_i = \frac{
            q_{\text{m},i}^\star\exp\left(A_i - \frac{H_i^\star}{T^\star}\right)\hat{f}_i^\star
        }{
            1 + \exp\left(A_i - \frac{H_i^\star}{T^\star}\right)\hat{f}_i^\star
              + \exp\left(A_j - \frac{H_j^\star}{T^\star}\right)\hat{f}_j^\star
        }
    :label: eq_lang_binary_dimensionless


Modules
-------

.. autoclass:: isotherm_models.binaryisotherm.BinaryIsotherm
    :members:

.. autoclass:: isotherm_models.binaryisotherm.BinaryLangmuir
    :members:
