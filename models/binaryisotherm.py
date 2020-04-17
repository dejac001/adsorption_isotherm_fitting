import matplotlib.pyplot as plt
import pyomo.environ as pyo
from models.unaryisotherm import UnaryIsotherm, LangmuirUnary


class BinaryIsotherm(UnaryIsotherm):
    r"""Base class for Binary Isotherms, inherits from UnaryIsotherm

    The following additional dimensionless variables are used in computations

    .. math::
        p_j^\star = \frac{p_j}{p_\text{ref}}
        :label: eq_pjs

    :param p_j: pressures (or fugacities) of component *j*
    :type p_j: list
    """
    def __init__(self, p_i, p_j, q_i, T, p_ref=None, **kwargs):
        if p_ref is None:
            p_ref = max(p_i + p_j)
        UnaryIsotherm.__init__(self, p_i, q_i, T, p_ref=p_ref, **kwargs)
        self.p_j = p_j
        assert len(self.p_j) == len(self.p_i), 'Inconsistent input data'
        self.p_j_star = [i/self.p_ref for i in self.p_j]

        self.unary_points = [i for i in self.points if self.p_j[i] < 1e-12]

    def plot_adsorption_surface(self):
        """plot surface of adsorption data

        .. todo::
            implement this
            helpful for debugging

        """
        pass

    def plot_unary(self, ax=None):
        assert len(self.unary_points) > 0, "No unary points found, doesnt make sense to plot it"
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_xlabel('p_i')
            ax.set_ylabel('q_i')

        vals = self.q_calc.extract_values()
        q_calc = [pyo.value(vals[i]) for i in self.unary_points]
        p_i = [self.p_i[i] for i in self.unary_points]
        q_i = [self.q_i[i] for i in self.unary_points]

        ax.plot(p_i, q_i, 'o', label='Raw Data units')
        ax.plot(p_i, q_calc, 'x', label='Fit units')


class BinaryLangmuir(BinaryIsotherm, LangmuirUnary):
    r"""Temperature-dependent extended unary Langmuir isotherm, expressed as

    .. math::
        q_i = \frac{M_ik_ip_i}{1 + k_i p_i + k_j p_j}
        :label: eq_lang_binary

    Arrhenius relationships are used for :math:`k_i`  and :math:`k_j`,
    and dimensionless variables are used as
    illustrated in :class:`models.unaryisotherm.LangmuirUnary`

    .. note::
        This isotherm is not equivalent to the conventional extended langmuir isotherm,
        because *both* :math:`k_i` and :math:`k_j` are fit simultaneously to binary data.

    For completeness, the relationships are repeated for the binary case below

    .. math::
        \begin{align}
            k_i &= k_{i,\infty}\exp\left(\frac{-\Delta H_i}{RT}\right)\\
            k_j &= k_{j,\infty}\exp\left(\frac{-\Delta H_j}{RT}\right)\\
        \end{align}

    The dimensionless variables to be fit are

    .. math::
        \begin{align}
            H_i^\star &= \frac{\Delta H_i}{R T_\text{ref}} \\
            H_j^\star &= \frac{\Delta H_j}{R T_\text{ref}} \\
            M_i^\star &= \frac{M_i}{q_\text{ref}} \\
            A_i &= \ln\left(k_{i,\infty} p_\text{ref}\right) \\
            A_j &= \ln\left(k_{j,\infty} p_\text{ref}\right) \\
        \end{align}

    So that Equation :eq:`eq_lang_binary` becomes

    .. math::
        \theta_i = \frac{M_i^\star\exp\left(A_i - \frac{H_i^\star}{T^\star}\right)p_i^\star}{1 + \exp\left(A_i - \frac{H_i^\star}{T^\star}\right)p_i^\star + \exp\left(A_j - \frac{H_j^\star}{T^\star}\right)p_j^\star}
        :label: eq_lang_binary_dimensionless


    :param H_i_star: :math:`H_i^\star`, dimensionless heat of adsorption of component *i*
    :type H_i_star: pyo.Var
    :param A_i: :math:`A_i`, dimensionless langmuir constant in logarithmic space
    :type A_i: pyo.Var
    :param H_j_star: :math:`H_j^\star`, dimensionless heat of adsorption of component *j*
    :type H_j_star: pyo.Var
    :param A_j: :math:`A_j`, dimensionless langmuir constant in logarithmic space
    :type A_j: pyo.Var
    :param M_i_star: :math:`M_i^\star`, dimensionless saturation loading
    :type M_i_star: pyo.Var
    :param M_i: langmuir saturaiton loading
    :type M_i: pyo.Expression
    :param k_i_inf: langmuir adsorption constant independent of temperature
    :type k_i_inf: pyo.Expression
    :param dH_i: heat of adsorption of component *i*
    :type dH_i: pyo.Expression
    :param k_j_inf: langmuir adsorption constant independent of temperature
    :type k_j_inf: pyo.Expression
    :param dH_j: heat of adsorption of component *j*
    :type dH_j: pyo.Expression
    """
    def __init__(self, *args, **kwargs):
        BinaryIsotherm.__init__(self, *args, **kwargs)

        # add variables
        self.H_i_star = pyo.Var(initialize=self.initial_guess_H_i_star())
        self.A_i = pyo.Var(initialize=self.initial_guess_A_i())
        self.H_j_star = pyo.Var(initialize=self.initial_guess_H_i_star())
        self.A_j = pyo.Var(initialize=self.initial_guess_A_i())
        self.M_i_star = pyo.Var(initialize=self.initial_guess_M_i_star(), bounds=(0., None))

        # add expressions for dimensional quantities
        self.M_i = pyo.Expression(expr=self.M_i_star*self.q_ref)
        self.k_i_inf = pyo.Expression(expr=pyo.exp(self.A_i) / self.p_ref)
        self.dH_i = pyo.Expression(expr=self.R*self.T_ref*self.H_i_star)
        self.k_j_inf = pyo.Expression(expr=pyo.exp(self.A_j) / self.p_ref)
        self.dH_j = pyo.Expression(expr=self.R*self.T_ref*self.H_j_star)

        self.isotherm_eq = pyo.Constraint(self.points, rule=BinaryLangmuir.isotherm_eq_rule)

    def isotherm_expression(self, point):
        """Isotherm expression in unit quantities, see Equation :eq:`eq_lang_binary`"""
        H_i = self.k_i_inf*pyo.exp(-self.dH_i/self.R/self.T[point])*self.p_i[point]
        if point in self.unary_points:
            return self.M_i * H_i / (1. + H_i)

        H_j = self.k_j_inf*pyo.exp(-self.dH_j/self.R/self.T[point])*self.p_j[point]
        return self.M_i * H_i / (1. + H_i + H_j)

    def dimensionless_isotherm_expression(self, point):
        """Dimensionless isotherm expression, see Equation :eq:`eq_lang_binary_dimensionless`"""
        K_i = pyo.exp(self.A_i - self.H_i_star / self.T_star[point]) * self.p_i_star[point]
        if point in self.unary_points:
            return self.M_i_star * K_i / (1. + K_i)

        K_j = pyo.exp(self.A_j - self.H_j_star / self.T_star[point]) * self.p_j_star[point]
        return self.M_i_star * K_i / (1. + K_i + K_j)
