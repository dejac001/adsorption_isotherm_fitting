import pyomo.environ as pyo
import matplotlib.pyplot as plt
from isotherm_models import solver as default_solver
from chem_util.chem_constants import gas_constant as R


def K_expr(k_inf, dH, T, f):
    r"""

    .. math::
        k_\infty \exp{\left(\frac{-\Delta H}{RT}\right)}f

    :param k_inf:
    :param dH:
    :param T:
    :param f:
    :return:
    """
    return k_inf * pyo.exp(dH / -R / T) * f


def K_star_expr(A, H_star, T_star, f_star):
    r"""

    .. math::
        \exp{\left(A - H^\star/T^\star\right)}f^\star

    :param A:
    :param H_star:
    :param T_star:
    :param f_star:
    :return:
    """
    return pyo.exp(A - H_star/T_star) * f_star


class UnaryIsotherm(pyo.ConcreteModel):
    r"""Base class for Unary Isotherms

    :param f_i: fugacities of component *i* (can be calculated assuming ideal gas or real gas)
    :type f_i: list
    :param q_i: loadings of component *i*
    :type q_i: list
    :param T: temperatures in [K], defaults to None
    :type T: list, optional
    :param q_ref: reference loading, defaults to maximum loading in :attr:`q_i`
    :type q_ref: float, optional
    :param f_ref: reference fugacity, defaults to maximum fugacity in :attr:`f_j`
    :type f_ref: float, optional
    :param T_ref: reference temperature, defaults to maximum temperature in :attr:`T`
    :type T_ref: float, optional
    :param points: state points at which a pressure and temperature are provided
    :type points: list, derived from input
    :param f_i_star: dimensionless fugacitities, calculated by Equation :eq:`eq_fis_unary`
    :type f_i_star: list, derived
    :param theta: dimensionless loadings, calculated by Equation :eq:`eq_theta`
    :type theta: list, derived
    :param T_star: dimensionless temperatures, calculated by Equation :eq:`eq_Ts`
    :type T_star: list, derived
    :param theta_calc: calculated dimensionless at each state point
    :type theta_calc: pyo.Var, derived from input
    :param objective: objective function to be minimized for isotherm fitting, calculated from :meth:`isotherm_models.unaryisotherm.UnaryIsotherm.objective_rule`
    :type objective: pyo.Objective, derived from input
    :param R2: coefficient of determination, see :meth:`isotherm_models.unaryisotherm.UnaryIsotherm.R2_rule`
    :type R2: pyo.Expression, derived
    :param q_calc: calculated loading in units
    :type q_calc: pyo.Expression, derived
    """

    def __init__(self, f_i, q_i, T, q_ref=None, f_ref=None, T_ref=None, **kwargs):
        """

        :param kwargs: passed to :code:`pyo.ConcreteModel`
        """
        pyo.ConcreteModel.__init__(self, **kwargs)

        assert len(f_i) == len(q_i), 'Inconsistent input data'
        assert len(f_i) == len(T), 'Inconsistent number of temperatures'

        self.f_i = f_i
        self.q_i = q_i
        self.T = T
        self.q_ref = q_ref
        self.f_ref = f_ref
        self.T_ref = T_ref
        self.points = pyo.Set(initialize=list(range(len(f_i))), ordered=True)

        # override defaults
        if self.q_ref is None:
            self.q_ref = float(max(q_i))
        if self.f_ref is None:
            self.f_ref = float(max(f_i))
        if self.T_ref is None:
            self.T_ref = float(max(T))

        self.f_i_star = [i / self.f_ref for i in self.f_i]
        self.theta = [i / self.q_ref for i in self.q_i]
        self.T_star = [i / self.T_ref for i in self.T]

        self.theta_calc = pyo.Var(self.points, initialize=1., bounds=(0., None))
        self.x_data_dimensionless = [(i, j) for i, j in zip(self.f_i_star, self.T_star)]

        self.q_calc = pyo.Expression(self.points, rule=UnaryIsotherm.q_calc_expr)
        self.R2 = pyo.Expression(expr=self.R2_rule())
        self.objective = pyo.Objective(expr=self.objective_rule(), sense=pyo.minimize)
        self.has_isotherm_variables = False

    def isotherm_eq_rule(self, point):
        """Constraint for dimensionless expression"""
        return self.theta_calc[point] == self.dimensionless_isotherm_expression(point)

    def display_fit_quality(self, **kwargs):
        self.R2.display(**kwargs)
        self.objective.display(**kwargs)

    def q_calc_expr(self, point):
        return self.theta_calc[point] * self.q_ref

    def dimensionless_isotherm_expression(self, point):
        raise NotImplementedError

    def R2_rule(self):
        """Calculate coefficient of determination squared, :math:`R^2` """
        # q_mean = pyo.summation(self.theta) / len(self.points)
        q_mean = sum(self.theta[i] for i in self.points) / len(self.points)
        ss_tot = sum((self.theta[i] - q_mean) * (self.theta[i] - q_mean) for i in self.points)
        ss_res = self.objective_rule()
        return 1 - ss_res / ss_tot

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
            (self.theta[i] - self.theta_calc[i]) * (self.theta[i] - self.theta_calc[i]) for i in self.points
        )

    def solve_scipy(self, loss='soft_l1', max_nfev=5000, bounds=None, **kwargs):
        from scipy.optimize import curve_fit
        if bounds is None:
            bounds = (-500, 500)
        return curve_fit(
            self.eval_dimensionless, self.x_data_dimensionless, self.theta,
            p0=self.initial_guess_vector(),
            loss=loss, max_nfev=max_nfev, bounds=bounds
        )

    def initial_guess_vector(self):
        raise NotImplementedError

    def eval(self, *args, **kwargs):
        raise NotImplementedError

    def eval_dimensionless(self, *args, **kwargs):
        raise NotImplementedError

    def solve(self, solver=default_solver, **kwargs):
        """Solve constraints subject to objective function

        :param solver: solver for solving model equations, defaults to pyo.SolverFactory('ipopt')
        :type solver: pyo.SolverFactory, optional
        :param kwargs: for solve argument
        """
        if not self.has_isotherm_variables:
            raise Exception('Need to add isotherm variables before solving; '
                            'otherwise will get perfect fit with no params')

        solver.solve(self, **kwargs)

    def plot_unary(self, ax=None, fig=None):
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
            ax.set_xlabel('hat_f_j')
            ax.set_ylabel('q_i')

        vals = self.q_calc.extract_values()
        q_calc = [pyo.value(vals[i]) for i in self.points]

        ax.plot(self.f_i, self.q_i, 'o', label='Raw Data units')
        ax.plot(self.f_i, q_calc, 'x', label='Fit units')
        return fig, ax

    def plot_unary_dimensionless(self, ax=None, fig=None):
        if fig is None:
            fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111)
            ax.set_xlabel('hat_f_j')
            ax.set_ylabel('q_i')

        vals = self.theta_calc.extract_values()
        theta_calc = [vals[i] for i in self.points]
        ax.plot(self.f_i_star, self.theta, 'o', label='Raw data dimensionless')
        ax.plot(self.f_i_star, theta_calc, 'x', label='Fit dimensionless')
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
    r"""Langmuir isotherm for unary mixture

    Isotherm is Equation :eq:`eq_lang_unary`.
    Dimensionless isotherm is Equation :eq:`eq_lang_unary_dimensionless`.
    Dimensionless variables to be fit are :math:`H_i^\star`, :math:`A_i`, and :math:`q_{\text{m},i}^\star`,
    as defined in Equations :eq:`H_i_star`, :eq:`q_mi_star`, and :eq:`A_i`, respectively.

    :param H_i_star: :math:`H_i^\star`, dimensionless heat of adsorption of component *i*
    :type H_i_star: pyo.Var
    :param A_i: :math:`A_i`, dimensionless langmuir constant in logarithmic space
    :type A_i: pyo.Var
    :param q_mi_star: :math:`q_{\text{m}i}^\star`, dimensionless saturation loading
    :type q_mi_star: pyo.Var
    :param q_mi: langmuir saturaiton loading
    :type q_mi: pyo.Expression
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
        self.q_mi_star = pyo.Var(initialize=self.initial_guess_q_mi_star())

        # add expressions for dimensional quantities
        self.q_mi = pyo.Expression(expr=self.q_mi_star * self.q_ref)
        self.k_i_inf = pyo.Expression(expr=pyo.exp(self.A_i) / self.f_ref)
        self.dH_i = pyo.Expression(expr=R * self.T_ref * self.H_i_star)

        self.isotherm_eq = pyo.Constraint(self.points, rule=LangmuirUnary.isotherm_eq_rule)
        self.has_isotherm_variables = True

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

    def initial_guess_q_mi_star(self):
        r"""Initial guess for :math:`q_mi^\star` variable

        If :math:`q_\text{ref}` is chosen to be the saturation loading, :math:`q_mi^\star` will be 1.
        Thus, we return 1 as initial guess
        """
        return 1.

    def initial_guess_vector(self):
        """:code:`p0` in scipy curve fit; initial guess for *dimensionless* parameters

        .. note::
            order here must be the same as last args in :meth:`.LangmuirUnary.eval_dimensionless`

        """
        return [
            self.initial_guess_q_mi_star(),
            self.initial_guess_A_i(),
            self.initial_guess_H_i_star()
        ]

    def eval(self, f_i, T, q_mi, k_i_inf, dH_i):
        """evaluate using generic types (any type)"""
        K = K_expr(k_i_inf, dH_i, T, f_i)
        return q_mi * K / (1. + K)

    def eval_pyomo(self, f_i, T):
        """evaluate using pyomo types (any type)"""
        return self.eval(
            f_i, T, self.q_mi, self.k_i_inf, self.dH_i
        )

    def eval_dimensionless(self, f_i_star, T_star, q_mi_star, A_i, H_i_star):
        K = K_star_expr(A=A_i, H_star=H_i_star, T_star=T_star, f_star=f_i_star)
        return q_mi_star * K / (1. + K)

    def eval_dimensionless_pyomo(self, f_i_star, T_star):
        return self.eval_dimensionless(
            f_i_star, T_star, self.q_mi_star, self.A_i, self.H_i_star
        )

    def isotherm_expression(self, point):
        """Isotherm expression in unit quantities, see Equation :eq:`eq_lang_unary`"""
        return self.eval_pyomo(self.f_i[point], self.T[point])

    def dimensionless_isotherm_expression(self, point):
        """Dimensionless isotherm expression, see Equation :eq:`eq_lang_unary_dimensionless`"""
        return self.eval_dimensionless_pyomo(self.f_i_star[point], self.T_star[point])

    def display_results(self, **kwargs):
        super().display_results(**kwargs)
        self.H_i_star.display(**kwargs)
        self.A_i.display(**kwargs)
        self.q_mi_star.display(**kwargs)
        self.q_mi.display(**kwargs)
        self.k_i_inf.display(**kwargs)
        self.dH_i.display(**kwargs)
