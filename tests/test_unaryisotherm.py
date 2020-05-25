from models.unaryisotherm import LangmuirUnary
import pandas as pd


def test_co2_bea():
    data = pd.read_csv('data_sets/CO2_BEA.csv')
    model = LangmuirUnary(data['P [atm]'], data['Q [mmol/g]'], data['T [K]'])
    model.solve()
    assert model.get_R2() > 0.99, 'R2 not above 0.99'
    assert model.get_objective() < 0.01, 'Objective function not below 0.01'


def test_N2_bea():
    data = pd.read_csv('data_sets/N2_BEA.csv')
    model = LangmuirUnary(data['P [atm]'], data['Q [mmol/g]'], data['T [K]'])
    model.solve()
    assert model.get_R2() > 0.99, 'R2 not above 0.99'
    assert model.get_objective() < 0.01, 'Objective function not below 0.01'


def test_H2S_MFI():
    data = pd.read_csv('data_sets/H2S_MFI.csv')
    model = LangmuirUnary(data['P [bar]'], data['Q [mmol/g]'], data['T [K]'])
    model.solve()
    assert model.get_R2() > 0.99, 'R2 not above 0.99'
    assert model.get_objective() < 0.01, 'Objective function not below 0.01'


def test_CH4_MFI():
    data = pd.read_csv('data_sets/CH4_MFI.csv')
    model = LangmuirUnary(data['P [bar]'], data['Q [mmol/g]'], data['T [K]'])
    model.solve()
    assert model.get_R2() > 0.99, 'R2 not above 0.99'
    assert model.get_objective() < 0.01, 'Objective function not below 0.01'

