import unittest
from models.unaryisotherm import LangmuirUnary
import pandas as pd


class MyTestCase(unittest.TestCase):
    def test_co2_bea(self):
        data = pd.read_csv('../data_sets/co2_bea.csv')
        model = LangmuirUnary(data['P [atm]'], data['Q [mmol/g]'], data['T [K]'])
        model.solve()
        print(model.pprint())
        print(model.display())


if __name__ == '__main__':
    unittest.main()
