import pyomo.environ as pyo
from .default_solver import default_solver
from .constants import R
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

    :param p_i: pressures (or concentrations) of component *i*
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
    """
    def __init__(self, p_i, q_i, T, q_ref=None, p_ref=None, T_ref=None):
        pyo.ConcreteModel.__init__(self)

        assert len(p_i) == len(q_i), 'Inconsistent input data'
        assert len(p_i) == len(T), 'Inconsistent number of temperatures'

        self.points = pyo.Set(initialize=list(range(len(p_i))), ordered=True)
        self.p_i = p_i
        self.q_i = q_i
        self.T = T
        self.q_ref = q_ref
        self.p_ref = p_ref
        self.T_ref = T_ref

        # override defaults
        if self.q_ref is None:
            self.q_ref = max(q_i)
        if self.p_ref is None:
            self.p_ref = max(p_i)
        if self.T_ref is None:
            self.T_ref = max(T)

        self.p_star = [i/self.p_ref for i in self.p_i]
        self.theta = [i/self.q_ref for i in self.q_i]
        self.T_star = [i/self.T_ref for i in self.T]

        self.theta_calc = pyo.Var(self.points, initialize=1., bounds=(0., None))

        self.q_calc = pyo.Expression(self.points, rule=UnaryIsotherm.q_calc_expr)
        self.objective = pyo.Objective(expr=self.objective_rule(), sense=pyo.minimize)

    def q_calc_expr(self, point):
        return self.theta_calc[point]*self.q_ref

    def isotherm_eq_rule(self, point):
        raise NotImplementedError

    def get_R2(self):
        """Calculate coefficient of determination squared, :math:`R^2` """
        q_mean = np.mean(self.theta)
        ss_tot = np.sum([(i-q_mean)*(i-q_mean) for i in self.theta])
        ss_res = pyo.value(self.objective)
        return 1 - ss_res/ss_tot

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

    def plot(self, ax=None):
        if ax is None:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('p_i')
            ax.set_ylabel('q_i')

        vals = self.q_calc.extract_values()
        q_calc = [pyo.value(vals[i]) for i in self.points]

        ax.plot(self.p_i, self.q_i, 'o', label='Raw Data units')
        ax.plot(self.p_i, q_calc, 'x', label='Fit units')

    def plot_dimensionless(self, ax=None):
        if ax is None:
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('p_i')
            ax.set_ylabel('q_i')

        vals = self.theta_calc.extract_values()
        theta_calc = [vals[i] for i in self.points]
        ax.plot(self.p_star, self.theta, 'o', label='Raw data dimensionless')
        ax.plot(self.p_star, theta_calc, 'x', label='Fit dimensionless')

    def get_objective(self):
        return pyo.value(self.objective)


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
    :param R: gas constant
    :type R: float
    """
    def __init__(self, *args, **kwargs):
        UnaryIsotherm.__init__(self, *args, **kwargs)

        # add variables
        self.H_i_star = pyo.Var(initialize=10.)
        self.A_i = pyo.Var(initialize=-3.)
        self.M_i_star = pyo.Var(initialize=0.5)

        # add expressions for dimensional quantities
        self.R = R
        self.M_i = pyo.Expression(expr=self.M_i_star*self.q_ref)
        self.k_i_inf = pyo.Expression(expr=pyo.exp(self.A_i) / self.p_ref)
        self.dH_i = pyo.Expression(expr=self.R*self.T_ref*self.H_i_star)

        self.isotherm_eq = pyo.Constraint(self.points, rule=LangmuirUnary.isotherm_eq_rule)

    def isotherm_expression(self, point):
        """Isotherm expression in unit quantities, see Equation :eq:`eq_lang_unary`"""
        K = self.k_i_inf*pyo.exp(-self.dH_i/self.R/self.T[point])*self.p_i[point]
        return self.M_i * K / (1. + K)

    def dimensionless_isotherm_expression(self, point):
        """Dimensionless isotherm expression, see Equation :eq:`eq_lang_unary_dimensionless`"""
        K = pyo.exp(self.A_i - self.H_i_star / self.T_star[point]) * self.p_star[point]
        return self.M_i_star * K / (1. + K)

    def isotherm_eq_rule(self, point):
        """Constraint for dimensionless expression"""
        return self.theta_calc[point] == self.dimensionless_isotherm_expression(point)
