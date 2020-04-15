import pyomo.environ as pyo
import numpy as np


class UnaryIsotherms:
    def __init__(self, conc, q):
        I = pyo.ConcreteModel()

        I.points = pyo.RangeSet(len(conc))

        def _init_conc(model, i):
            return conc[i-1]
        I.conc = pyo.Param(I.points, initialize=_init_conc)

        def _init_q(model, i):
            return q[i-1]
        I.loadings = pyo.Param(I.points, initialize=_init_q)
        I.k1 = pyo.Var()
        I.M1 = pyo.Var(within=pyo.NonNegativeReals)
        self.model = I

    def calculateR2(self):
        q = [i for i in self.model.loadings.values()]
        q_mean = np.mean(q)
        ss_tot = np.sum([(i-q_mean)*(i-q_mean) for i in q])
        ss_res = pyo.value(self.model.objective)
        self.R2 = 1 - ss_res/ss_tot


class LangmuirUnary(UnaryIsotherms):
    def __init__(self, conc, q):
        UnaryIsotherms.__init__(self, conc, q)

        def _my_objective(model):
            residuals = (
                model.loadings[i] - model.k1*model.M1*model.conc[i]/(1 + model.k1*model.conc[i])
                for i in model.points
            )
            return sum(k*k for k in residuals)
        self.model.objective = pyo.Objective(rule=_my_objective, sense=pyo.minimize)


class BET_Unary(UnaryIsotherms):
    def __init__(self, conc, q):
        UnaryIsotherms.__init__(self, conc, q)
        self.model.k2 = pyo.Var()

        def _my_objective(m):
            residuals = (
                m.loadings[i] -
                m.M1*m.k1*m.conc[i]/(1-m.k2*m.conc[i])/(1+m.k1*m.conc[i] - m.k2*m.conc[i])
                for i in m.points
            )
            return sum(k*k for k in residuals)

        self.model.objective = pyo.Objective(rule=_my_objective, sense=pyo.minimize)


class Quadratic_Unary(UnaryIsotherms):
    def __init__(self, conc, q):
        UnaryIsotherms.__init__(self, conc, q)
        self.model.k2 = pyo.Var()

        def _my_objective(m):
            residuals = (
                m.loadings[i] -
                m.M1*(m.k1+2*m.k2*m.conc[i])*m.conc[i]/(1 + m.k1*m.conc[i]+m.k2*m.conc[i]*m.conc[i])
                for i in m.points
            )
            return sum(k*k for k in residuals)

        self.model.objective = pyo.Objective(rule=_my_objective, sense=pyo.minimize)


class DualSite(UnaryIsotherms):
    def __init__(self, conc, q):
        UnaryIsotherms.__init__(self, conc, q)
        self.model.k2 = pyo.Var()
        self.model.M2 = pyo.Var(within=pyo.NonNegativeReals)

        def _my_objective(m):
            residuals = (
                m.loadings[i] -
                (m.M1*m.k1*m.conc[i]/(1 + m.k1*m.conc[i]) + m.M2*m.k2*m.conc[i]/(1+m.k2*m.conc[i]))
                for i in m.points
            )
            return sum(k*k for k in residuals)

        self.model.objective = pyo.Objective(rule=_my_objective, sense=pyo.minimize)


class DualSite_Quadratic(UnaryIsotherms):
    def __init__(self, conc, q):
        UnaryIsotherms.__init__(self, conc, q)
        self.model.k2 = pyo.Var()
        self.model.k3 = pyo.Var()
        self.model.k4 = pyo.Var()
        self.model.M2 = pyo.Var(bounds=(0,None))

        def _my_objective(m):
            residuals = (
                m.loadings[i] -
                m.M1*(m.k1+2*m.k2*m.conc[i])*m.conc[i]/(1 + m.k1*m.conc[i]+m.k2*m.conc[i]*m.conc[i])
                - m.M2*(m.k3+2*m.k4*m.conc[i])*m.conc[i]/(1 + m.k3*m.conc[i]+m.k4*m.conc[i]*m.conc[i])
                for i in m.points
            )
            return sum(k*k for k in residuals)

        self.model.objective = pyo.Objective(rule=_my_objective, sense=pyo.minimize)
