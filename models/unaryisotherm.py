import pyomo.environ as pyo
import matplotlib.pyplot as plt
from .default_solver import default_solver
from chem_util.chem_constants import gas_constant as R
import numpy as np


class UnaryIsotherm(pyo.ConcreteModel):
    r"""Base class for Unary Isotherms

    The following dimensionless variables are used in computations

    .. math::
        \theta_i = \frac{q_i}{q_\text{ref}}
        :label: eq_theta

    .. math::
        p_i^\star = \frac{p_i}{p_\text{ref}}
        :label: eq_pis

    .. math::
        T^\star = \frac{T}{T_\text{ref}}
        :label: eq_Ts

    :param p_i: pressures (or fugacities) of component *i*
    :type p_i: list
    :param q_i: loadings of component *i*
    :type q_i: list
    :param T: temperatures in [K], defaults to None
    :type T: list, optional
    :param points: state points at which a pressure and temperature are provided
    :type points: list, derived from input
    :param q_calc: calculated loding at each state point
    :type q_calc: pyo.Var, derived from input
    :param objective: objective function to be minimized for isotherm fitting, calculated from :meth:`models.unaryisotherm.UnaryIsotherm.objective_rule`
    :type objective: pyo.Objective, derived from input
    :param q_ref: reference loading, defaults to maximum loading in :attr:`q_i`
    :type q_ref: float, optional
    :param p_ref: reference pressure, defaults to maximum pressure in :attr:`p_i`
    :type p_ref: float, optional
    :param T_ref: reference temperature, defaults to maximum temperature in :attr:`T`
    :type T_ref: float, optional
    :param q_calc: calculated loading in units
    :type q_calc: pyo.Expression, derived
    :param R: gas constant, set to 8.314 J/mol/K [=] m**3*Pa/mol/K
    :type R: float, hard-coded
    """
    def __init__(self, p_i, q_i, T, q_ref=None, p_ref=None, T_ref=None, **kwargs):
        """

        :param kwargs: passed to :code:`pyo.ConcreteModel`
        """
        pyo.ConcreteModel.__init__(self, **kwargs)

        assert len(p_i) == len(q_i), 'Inconsistent input data'
        assert len(p_i) == len(T), 'Inconsistent number of temperatures'

        self.R = R
        self.points = pyo.Set(initialize=list(range(len(p_i))), ordered=True)
        self.p_i = p_i
        self.q_i = q_i
        self.T = T
        self.q_ref = q_ref
        self.p_ref = p_ref
        self.T_ref = T_ref

        # override defaults
        if self.q_ref is None:
            self.q_ref = float(max(q_i))
        if self.p_ref is None:
            self.p_ref = float(max(p_i))
        if self.T_ref is None:
            self.T_ref = float(max(T))

        self.p_i_star = [i / self.p_ref for i in self.p_i]
        self.theta = [i/self.q_ref for i in self.q_i]
        self.T_star = [i/self.T_ref for i in self.T]

        self.theta_calc = pyo.Var(self.points, initialize=1., bounds=(0., None))

        self.q_calc = pyo.Expression(self.points, rule=UnaryIsotherm.q_calc_expr)
        self.R2 = pyo.Expression(expr=self.R2_rule())
        self.objective = pyo.Objective(expr=self.objective_rule(), sense=pyo.minimize)

    def isotherm_eq_rule(self, point):
        """Constraint for dimensionless expression"""
        return self.theta_calc[point] == self.dimensionless_isotherm_expression(point)

    def display_fit_quality(self, **kwargs):
        self.R2.display(**kwargs)
        self.objective.display(**kwargs)

    def q_calc_expr(self, point):
        return self.theta_calc[point]*self.q_ref

    def dimensionless_isotherm_expression(self, point):
        raise NotImplementedError

    def R2_rule(self):
        """Calculate coefficient of determination squared, :math:`R^2` """
        # q_mean = pyo.summation(self.theta) / len(self.points)
        q_mean = sum(self.theta[i] for i in self.points) / len(self.points)
        ss_tot = sum((self.theta[i] - q_mean)*(self.theta[i] - q_mean) for i in self.points)
        ss_res = self.objective_rule()
        return 1 - ss_res/ss_tot

    def get_R2(self):
        return pyo.value(self.R2)

    def get_objective(self):
        return pyo.value(self.objective)

    def objective_rule(self):
        r"""Sum of squared errors between calculated loading and predicted loading

        .. math::
            \sum_i \left(\theta_i^{\text{raw}}-\theta_i^{\text{calc}}\right)^2

        where *raw* denotes the raw data obtained by experiment or molecular simulation
        and *calc* denotes the data calculated from the isotherm function

        """
        return sum(
            (self.theta[i]-self.theta_calc[i])*(self.theta[i]-self.theta_calc[i]) for i in self.points
        )

    def solve(self, solver=None, **kwargs):
        """Solve constraints subject to objective function

        :param solver: solver for solving model equations, defaults to pyo.SolverFactory('ipopt')
        :type solver: pyo.SolverFactory, optional
        :param kwargs: for solve argument
        """
        if solver is None:
            solver = default_solver

        solver.solve(self, **kwargs)

    def plot_unary(self, ax=None, fig=None):
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
            ax.set_xlabel('p_i')
            ax.set_ylabel('q_i')

        vals = self.q_calc.extract_values()
        q_calc = [pyo.value(vals[i]) for i in self.points]

        ax.plot(self.p_i, self.q_i, 'o', label='Raw Data units')
        ax.plot(self.p_i, q_calc, 'x', label='Fit units')
        return fig, ax

    def plot_unary_dimensionless(self, ax=None, fig=None):
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
            ax.set_xlabel('p_i')
            ax.set_ylabel('q_i')

        vals = self.theta_calc.extract_values()
        theta_calc = [vals[i] for i in self.points]
        ax.plot(self.p_i_star, self.theta, 'o', label='Raw data dimensionless')
        ax.plot(self.p_i_star, theta_calc, 'x', label='Fit dimensionless')
        return fig, ax

    def plot_comparison_base(self, q_i, q_calc, xlabel, ylabel, fig=None, ax=None, marker='x', **kwargs):
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        if not kwargs:
            kwargs = dict(color='C1')

        ax.plot(q_i, q_calc, marker, **kwargs)

        max_q = max(q_i + q_calc)
        ax.plot([0., max_q], [0., max_q], '--', color='black')
        return fig, ax

    def plot_comparison_units(self, **kwargs):
        xlabel = 'Loading from Raw Data, with units'
        ylabel = 'Predicted Loading from Fit, with units'

        vals = self.q_calc.extract_values()
        q_calc = [pyo.value(vals[i]) for i in self.points]
        return self.plot_comparison_base(self.q_i, q_calc, xlabel, ylabel, **kwargs)

    def plot_comparison_dimensionless(self, **kwargs):
        xlabel = 'Loading from Raw Data, dimensionless'
        ylabel = 'Predicted Loading from Fit, dimensionless'

        vals = self.theta_calc.extract_values()
        theta_calc = [vals[i] for i in self.points]
        return self.plot_comparison_base(self.theta, theta_calc, xlabel, ylabel, **kwargs)

    def display_results(self, **kwargs):
        self.display_fit_quality(**kwargs)


class LangmuirUnary(UnaryIsotherm):
    r"""Temperature-dependent unary Langmuir isotherm, expressed as

    .. math::
        q_i = \frac{M_ik_ip_i}{1 + k_i p_i}
        :label: eq_lang_unary

    where an Arrhenius relationship for :math:`k_i` is assumed as

    .. math::
        k_i = k_{i,\infty}\exp\left(\frac{-\Delta H_i}{RT}\right)


    Introducing the dimensionless parameters in Equation :eq:`eq_theta` :eq:`eq_pis` and :eq:`eq_Ts`,
    the dimensionless variables to be fit are

    .. math::
        \begin{align}
            H_i^\star &= \frac{\Delta H_i}{R T_\text{ref}} \\
            M_i^\star &= \frac{M_i}{q_\text{ref}} \\
            A_i &= \ln\left(k_{i,\infty} p_\text{ref}\right) \\
        \end{align}

    So that Equation :eq:`eq_lang_unary` becomes

    .. math::
        \theta_i = \frac{M_i^\star\exp\left(A_i - \frac{H_i^\star}{T^\star}\right)p_i^\star}{1 + \exp\left(A_i - \frac{H_i^\star}{T^\star}\right)p_i^\star}
        :label: eq_lang_unary_dimensionless


    :param H_i_star: :math:`H_i^\star`, dimensionless heat of adsorption of component *i*
    :type H_i_star: pyo.Var
    :param A_i: :math:`A_i`, dimensionless langmuir constant in logarithmic space
    :type A_i: pyo.Var
    :param M_i_star: :math:`M_i^\star`, dimensionless saturation loading
    :type M_i_star: pyo.Var
    :param M_i: langmuir saturaiton loading
    :type M_i: pyo.Expression
    :param k_i_inf: langmuir adsorption constant independent of temperature
    :type k_i_inf: pyo.Expression
    :param dH_i: heat of adsorption of component *i*
    :type dH_i: pyo.Expression
    """
    def __init__(self, *args, **kwargs):
        UnaryIsotherm.__init__(self, *args, **kwargs)

        # add variables
        self.H_i_star = pyo.Var(initialize=self.initial_guess_H_i_star())
        self.A_i = pyo.Var(initialize=self.initial_guess_A_i())
        self.M_i_star = pyo.Var(initialize=self.initial_guess_M_i_star())

        # add expressions for dimensional quantities
        self.M_i = pyo.Expression(expr=self.M_i_star*self.q_ref)
        self.k_i_inf = pyo.Expression(expr=pyo.exp(self.A_i) / self.p_ref)
        self.dH_i = pyo.Expression(expr=self.R*self.T_ref*self.H_i_star)

        self.isotherm_eq = pyo.Constraint(self.points, rule=LangmuirUnary.isotherm_eq_rule)

    def initial_guess_H_i_star(self):
        r"""Initial guess for :math:`H_i^\star` variable

        This value of 10 corresponds to an absolute value for heat of adsorption of :math:`10RT`
        which is approximately 25 kJ/mol
        """
        return -10.

    def initial_guess_A_i(self):
        r"""Initial guess for :math:`A_i` variable


        .. todo::
            Come up with logical initial guess

        """
        return -1.

    def initial_guess_M_i_star(self):
        r"""Initial guess for :math:`M_i^\star` variable

        If :math:`q_\text{ref}` is chosen to be the saturation loading, :math:`M_i^\star` will be 1.
        Thus, we return 1 as initial guess
        """
        return 1.

    def isotherm_expression(self, point):
        """Isotherm expression in unit quantities, see Equation :eq:`eq_lang_unary`"""
        K = self.k_i_inf*pyo.exp(-self.dH_i/self.R/self.T[point])*self.p_i[point]
        return self.M_i * K / (1. + K)

    def dimensionless_isotherm_expression(self, point):
        """Dimensionless isotherm expression, see Equation :eq:`eq_lang_unary_dimensionless`"""
        K = pyo.exp(self.A_i - self.H_i_star / self.T_star[point]) * self.p_i_star[point]
        return self.M_i_star * K / (1. + K)

    def display_results(self, **kwargs):
        super().display_results(**kwargs)
        self.H_i_star.display(**kwargs)
        self.A_i.display(**kwargs)
        self.M_i_star.display(**kwargs)
        self.M_i.display(**kwargs)
        self.k_i_inf.display(**kwargs)
        self.dH_i.display(**kwargs)
