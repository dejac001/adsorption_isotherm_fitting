
Real Adsorbed Solution Theory
=============================

The *Gibbs adsorption isotherm* is

.. math::
    -a \mathrm{d}\Pi + \sum_i x_i \mathrm{d}\mu_i = 0
    :label: gibbs_iso

where :math:`a` is the surface area per mole of adsorbate,
:math:`\Pi` is the spreading pressure,
:math:`x_i` is the adsorbed mole fraction of component *i*,
and :math:`\mu_i` is the adsorbed-phase chemical potential of component *i*.

For change in equilibrium conditions,

.. math::
    \mathrm{d}\mu_i = \mathrm{d}\mu_i^\text{g} = RT \mathrm{d} \ln{\hat{f}_i^\text{g}}
    :label: mu

where :math:`\hat{f}_i^\text{g}` is the fugacity of component *i* in the gas phase.
And the substituting surface area of the adsorbent is


.. math::
    A = a \sum_i q_i
    :label: A

where :math:`q_i` is the loading of component *i*.
Substituting Equations :eq:`mu` and :eq:`A` into Equation :eq:`gibbs_iso` yields

.. math::
    \frac{A}{RT} \mathrm{d}\Pi = \sum_i q_i \mathrm{d} \ln{\hat{f}_i^\text{g}}
    :label: fugacity_iso

If we have a good description of the multicomponent isotherms,

.. math::
    q_i = F(\{\hat{f}_k\})

where :math:`F` is an isotherm function,
Equation :eq:`fugacity_iso` can be simplified to

.. math::
    \frac{A\Pi}{RT} = \sum_i \int_0^{\hat{f}_i^\text{g}}\frac{q_i}{f_i^\prime}\mathrm{d}f_i^\prime

where :math:`f_i^\prime` is a dummy variable for integration.