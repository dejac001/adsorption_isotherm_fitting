import unittest
from models.unaryisotherm import LangmuirUnary
import pandas as pd


class MyTestCase(unittest.TestCase):

    def test_co2_bea(self):
        data = pd.read_csv('../data_sets/CO2_BEA.csv')
        model = LangmuirUnary(data['P [atm]'], data['Q [mmol/g]'], data['T [K]'])
        model.solve()
        self.assertGreater(model.get_R2(), 0.99)
        self.assertLess(model.get_objective(), 0.01)

    def test_N2_bea(self):
        data = pd.read_csv('../data_sets/N2_BEA.csv')
        model = LangmuirUnary(data['P [atm]'], data['Q [mmol/g]'], data['T [K]'])
        model.solve()
        self.assertGreater(model.get_R2(), 0.99)
        self.assertLess(model.get_objective(), 0.01)

    def test_H2S_MFI(self):
        data = pd.read_csv('../data_sets/H2S_MFI.csv')
        model = LangmuirUnary(data['P [bar]'], data['Q [mmol/g]'], data['T [K]'])
        model.solve()
        self.assertGreater(model.get_R2(), 0.99)
        self.assertLess(model.get_objective(), 0.01)

    def test_CH4_MFI(self):
        data = pd.read_csv('../data_sets/CH4_MFI.csv')
        model = LangmuirUnary(data['P [bar]'], data['Q [mmol/g]'], data['T [K]'])
        model.solve()
        self.assertGreater(model.get_R2(), 0.99)
        self.assertLess(model.get_objective(), 0.001)


if __name__ == '__main__':
    unittest.main()
