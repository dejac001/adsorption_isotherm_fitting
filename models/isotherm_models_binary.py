import numpy as np
import pyomo.environ as pyo


class BinaryIsotherms:
    def __init__(self, isodata, sorbate, otherSorbate):
        I = pyo.ConcreteModel()
        I.points = pyo.RangeSet(0, len(isodata['q'][sorbate])-1)  #for rangeset, last point is also included
        I.main_sorbate = pyo.Param(initialize=sorbate)
        I.otherSorbate = pyo.Param(initialize=otherSorbate)
        I.sorbates = pyo.Set(initialize=[sorbate, otherSorbate])

        def _init_pressures(model, sorbate, point):
            return isodata['p'][sorbate][point]
        I.p = pyo.Param(I.sorbates, I.points, initialize=_init_pressures)

        I.loadings = pyo.Param(I.points, initialize={i:isodata['q'][sorbate][i] for i in I.points})
        I.max_loading = pyo.Param(initialize=np.max(isodata['q'][sorbate]))

        def _init_Q(m, point):
            return m.loadings[point] * 0.5
        I.Q_predicted = pyo.Var(I.points, initialize=_init_Q)  # usually get problems here when set bounds

        self.model = I


def langmuir_site(M, k_i, k_j, c_i, c_j):
    return M * c_i / (1 + k_i * c_i + k_j * c_j)


def langmuir_deriv_ii(M, k_i, k_j, c_i, c_j):
    """
    derivative of langmuir_site with respect to c_i
    """
    denominator = (1 + k_i * c_i + k_j * c_j)
    return M * (1 + k_j * c_j) / denominator / denominator


def langmuir_deriv_ij(M, k_i, k_j, c_i, c_j):
    """
    derivative of langmuir_site with respect to c_j
    """
    denominator = (1 + k_i * c_i + k_j * c_j)
    return -M * c_i * k_j / denominator / denominator


class LangmuirSites(BinaryIsotherms):

    def __init__(self, model_with_data, num_sites=1):
        self.model = model_with_data
        self.n_sites = num_sites
        self.model.sites = pyo.RangeSet(num_sites)

        # def _init_k(model, site, sorbate):
        #     if sorbate == 'butane':
        #         factor = -1.
        #     else:
        #         factor = 1.
        #     return 1.*site*factor
        self.model.k_self = pyo.Var(self.model.sites, initialize=0., bounds=(-1e12, 1e12))
        self.model.k_other = pyo.Var(self.model.sites, initialize=0., bounds=(-1e12, 1e12))
        self.model.scale_factor = pyo.Param(initialize=1e3)
        self.R2 = None

        def _init_C(model, sorbate, point):
            return model.p[sorbate,point]*model.scale_factor
        self.model.c = pyo.Param(self.model.sorbates, self.model.points, initialize=_init_C)

        def _M_bounds(m):
            return m.max_loading/1000., 1e12
        self.model.M = pyo.Var(self.model.sites, bounds=_M_bounds,
                               initialize=0.01)

        def _isotherm(model, data_point):
            predicted_loading = sum(
                langmuir_site(model.M[site], model.k_self[site], model.k_other[site],
                              model.c[model.main_sorbate.value, data_point],
                              model.c[model.otherSorbate.value, data_point])
                for site in model.sites
            )
            return model.Q_predicted[data_point] == predicted_loading
        self.model.isotherm_constraint = pyo.Constraint(self.model.points, rule=_isotherm)

        def _sum_residuals(model):
            residuals = (
                model.loadings[i] - model.Q_predicted[i] for i in model.points
            )
            return sum(k*k for k in residuals)
        self.model.sum_residuals = pyo.Expression(rule=_sum_residuals)

        def _chi_squared(model):
            return sum(
                (model.loadings[i] - model.Q_predicted[i]) * (model.loadings[i] - model.Q_predicted[i])
                / model.Q_predicted[i]
                for i in self.model.points
            )
        self.model.chi_squared = pyo.Expression(rule=_chi_squared)

        def _my_objective(model):
            return 0.01 * model.sum_residuals

        self.model.objective = pyo.Objective(rule=_my_objective, sense=pyo.minimize)

    def calculateR2(self):
        q = [i for i in self.model.loadings.values()]
        q_mean = np.mean(q)
        ss_tot = np.sum([(i-q_mean)*(i-q_mean) for i in q])
        ss_res = pyo.value(self.model.sum_residuals())
        self.R2 = 1 - ss_res/ss_tot

    def parameter_string(self):
        string = ''
        for site in self.model.sites:
            unscaled_params = [k*self.model.scale_factor for k in [self.model.M[site], self.model.k_self[site],
                                                                   self.model.k_other[site]]]
            q_n, k_in, k_jn = map(pyo.value, unscaled_params)
            string += ' %e %e %e' % (q_n, k_in, k_jn)
        return string

    def deactivate(self):
        for attr in (self.model.k_self, self.model.k_other, self.model.M, self.model.sites,
                     self.model.isotherm_constraint,
                     self.model.c, self.model.c_index, self.model.scale_factor):
            self.model.del_component(attr)
        # if self.n_sites > 1:
        #     self.model.del_component(self.model.q_site)


def sips_binary(M, k, n, C):
    k_main, k_other = k
    n_main, n_other = n
    c_main, c_other = C
    return (
        M*k_main*c_main**n_main / ( 1 + k_main*c_main**n_main + k_other*c_other**n_other )
    )
