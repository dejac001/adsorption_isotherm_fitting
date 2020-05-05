"""Compare binary fitting of isotherm to extended langmuir isotherm"""

import pyomo.environ as pyo
import matplotlib.pyplot as plt
import pandas as pd
from models.unaryisotherm import LangmuirUnary
from models.binaryisotherm import BinaryLangmuir


def co2():
    # step 1: get CO2 data
    data = pd.read_csv('../data_sets/CO2_BEA.csv')
    P_i = data['P [atm]']*101325  # convert to Pa -- si units
    # step 2: create CO2 Model
    co2_model = LangmuirUnary(P_i, data['Q [mmol/g]'], data['T [K]'], name='CO2')
    # step 3: solve CO2 Model
    co2_model.solve()
    # step 4: print results to file
    with open(co2_model.name + '.output', 'w') as f:
        co2_model.display_results(ostream=f)
    # step 7: plot results and save
    fig = plt.figure()
    fig, ax = co2_model.plot_unary(fig=fig)
    ax.legend()
    fig.savefig('CO2_example.png')


def n2():
    # step 1: get CO2 data
    data = pd.read_csv('../data_sets/N2_BEA.csv')
    P_i = data['P [atm]']*101325  # convert to Pa -- si units
    # step 2: create CO2 Model
    co2_model = LangmuirUnary(P_i, data['Q [mmol/g]'], data['T [K]'], name='N2')
    # step 3: solve CO2 Model
    co2_model.solve()
    # step 4: print results to file
    with open(co2_model.name + '.output', 'w') as f:
        co2_model.display_results(ostream=f)
    # step 7: plot results and save
    fig = plt.figure()
    fig, ax = co2_model.plot_unary(fig=fig)
    ax.legend()
    fig.savefig('N2.png')


if __name__ == '__main__':
    co2()
    n2()