import matplotlib.pyplot as plt
import pyomo.environ as pyo
from isotherm_models.unaryisotherm import UnaryIsotherm, LangmuirUnary
from chem_util.chem_constants import gas_constant as R


class BinaryIsotherm(UnaryIsotherm):
    r"""Base class for Binary Isotherms, inherits from UnaryIsotherm

    The following additional dimensionless variables are used in computations

    :param hat_f_j: mixture fugacities of component *i*
    :type hat_f_j: list
    :param hat_f_j: mixture fugacities of component *j*
    :type hat_f_j: list
    :param q_i: loadings of component *i*
    :type q_i: list
    :param T: temperatures in [K], defaults to None
    :type T: list, optional
    :param points: state points at which a pressure and temperature are provided
    :type points: list, derived from input
    :param hat_f_i_star: dimensionless fugacitities of component *i*, calculated by Equation :eq:`eq_fis_binary`
    :type hat_f_i_star: list, derived
    :param hat_f_j_star: dimensionless fugacitities of component *j*, calculated by Equation :eq:`eq_fjs_binary`
    :type hat_f_j_star: list, derived
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
    :param unary_points: points where only *i* is present, derived from where :math:`\hat{f}_j < 1\times10^{-12}`
    :type unary_points: list, derived

    """
    def __init__(self, hat_f_i, hat_f_j, q_i, T, f_ref=None, **kwargs):
        if f_ref is None:
            f_ref = max(hat_f_i + hat_f_j)
        UnaryIsotherm.__init__(self, hat_f_i, q_i, T, f_ref=f_ref, **kwargs)
        self.hat_f_i = self.f_i[:]
        self.hat_f_i_star = self.f_i_star[:]
        del self.f_i  # keep notation clear use hat for mixture
        del self.f_i_star  # keep notation clear use hat for mixture
        self.hat_f_j = hat_f_j
        assert len(self.hat_f_j) == len(self.hat_f_i), 'Inconsistent input data'
        self.hat_f_j_star = [i / self.f_ref for i in self.hat_f_j]

        self.unary_points = [i for i in self.points if self.hat_f_j[i] < 1e-12]

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
            ax.set_xlabel('hat_f_j')
            ax.set_ylabel('q_i')

        vals = self.q_calc.extract_values()
        q_calc = [pyo.value(vals[i]) for i in self.unary_points]
        p_i = [self.hat_f_i[i] for i in self.unary_points]
        q_i = [self.q_i[i] for i in self.unary_points]

        ax.plot(p_i, q_i, 'o', label='Raw Data units')
        ax.plot(p_i, q_calc, 'x', label='Fit units')


class BinaryLangmuir(BinaryIsotherm, LangmuirUnary):
    r"""Temperature-dependent extended unary Langmuir isotherm, expressed as

    Isotherm is Equation :eq:`eq_lang_binary`.
    Dimensionless isotherm is Equation :eq:`eq_lang_binary_dimensionless`.
    Dimensionless variables to be fit are :math:`H_i^\star`, :math:`A_i`, :math:`q_{\text{m},i}^\star`,
    :math:`H_j_star`, and :math:`A_j`,
    as defined in Equations :eq:`H_i_star`, :eq:`q_mi_star`, :eq:`A_i`,
    :eq:`H_j_star`, and :eq:`A_j`, respectively.

    :param H_i_star: :math:`H_i^\star`, dimensionless heat of adsorption of component *i*
    :type H_i_star: pyo.Var
    :param A_i: :math:`A_i`, dimensionless langmuir constant in logarithmic space
    :type A_i: pyo.Var
    :param H_j_star: :math:`H_j^\star`, dimensionless heat of adsorption of component *j*
    :type H_j_star: pyo.Var
    :param A_j: :math:`A_j`, dimensionless langmuir constant in logarithmic space
    :type A_j: pyo.Var
    :param q_mi_star: :math:`q_{\text{m},i}^\star`, dimensionless saturation loading
    :type q_mi_star: pyo.Var
    :param q_mi: langmuir saturation loading
    :type q_mi: pyo.Expression
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
        self.q_mi_star = pyo.Var(initialize=self.initial_guess_M_i_star(), bounds=(0., None))

        # add expressions for dimensional quantities
        self.q_mi = pyo.Expression(expr=self.q_mi_star * self.q_ref)
        self.k_i_inf = pyo.Expression(expr=pyo.exp(self.A_i) / self.f_ref)
        self.dH_i = pyo.Expression(expr=R*self.T_ref*self.H_i_star)
        self.k_j_inf = pyo.Expression(expr=pyo.exp(self.A_j) / self.f_ref)
        self.dH_j = pyo.Expression(expr=R*self.T_ref*self.H_j_star)

        self.isotherm_eq = pyo.Constraint(self.points, rule=BinaryLangmuir.isotherm_eq_rule)

    def isotherm_expression(self, point):
        """Isotherm expression in unit quantities, see Equation :eq:`eq_lang_binary`"""
        H_i = self.k_i_inf*pyo.exp(-self.dH_i/R/self.T[point])*self.hat_f_i[point]
        if point in self.unary_points:
            return self.q_mi * H_i / (1. + H_i)

        H_j = self.k_j_inf*pyo.exp(-self.dH_j/R/self.T[point])*self.hat_f_j[point]
        return self.q_mi * H_i / (1. + H_i + H_j)

    def dimensionless_isotherm_expression(self, point):
        """Dimensionless isotherm expression, see Equation :eq:`eq_lang_binary_dimensionless`"""
        K_i = pyo.exp(self.A_i - self.H_i_star / self.T_star[point]) * self.hat_f_i_star[point]
        if point in self.unary_points:
            return self.q_mi_star * K_i / (1. + K_i)

        K_j = pyo.exp(self.A_j - self.H_j_star / self.T_star[point]) * self.hat_f_j_star[point]
        return self.q_mi_star * K_i / (1. + K_i + K_j)

    def display_results(self, **kwargs):
        super().display_results(**kwargs)
        self.H_j_star.display(**kwargs)
        self.A_j.display(**kwargs)
        self.k_j_inf.display(**kwargs)
        self.dH_j.display(**kwargs)