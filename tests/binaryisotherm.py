import unittest
from models.binaryisotherm import BinaryLangmuir
import pandas as pd


class MyTestCase(unittest.TestCase):

    def test_H2S_MFI(self):
        df = pd.read_csv('../data_sets/CH4_H2S_MFI_binary_with_fugacity.csv')
        p_i, p_j, q_i, T = df['fugacity H2S [Pa]'], df['fugacity CH4 [Pa]'], df['Q H2S [mmol/g]'], df['T [K]']
        # clean data to remove state points at which no H2S is adsorbed
        points = [i for i in range(len(q_i)) if q_i[i] > 0.]
        p_i = [p_i[i] for i in points]
        p_j = [p_j[i] for i in points]
        q_i = [q_i[i] for i in points]
        T = [T[i] for i in points]
        model = BinaryLangmuir(p_i, p_j, q_i, T)
        model.solve()
        self.assertLess(model.get_objective(), 1.0)
        self.assertGreater(model.get_R2(), 0.99)

    def test_CH4_MFI(self):
        df = pd.read_csv('../data_sets/CH4_H2S_MFI_binary_with_fugacity.csv')
        p_i, p_j, q_i, T = df['fugacity CH4 [Pa]'], df['fugacity H2S [Pa]'], df['Q CH4 [mmol/g]'], df['T [K]']
        # clean data to remove state points at which no H2S is adsorbed
        points = [i for i in range(len(q_i)) if q_i[i] > 0.]
        p_i = [p_i[i] for i in points]
        p_j = [p_j[i] for i in points]
        q_i = [q_i[i] for i in points]
        T = [T[i] for i in points]
        model = BinaryLangmuir(p_i, p_j, q_i, T)
        model.solve()
        self.assertLess(model.get_objective(), 1.0)
        self.assertGreater(model.get_R2(), 0.99)


if __name__ == '__main__':
    unittest.main()
