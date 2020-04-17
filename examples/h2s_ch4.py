"""Compare binary fitting of isotherm to extended langmuir isotherm"""

import pyomo.environ as pyo
import matplotlib.pyplot as plt
import pandas as pd
from models.unaryisotherm import LangmuirUnary
from models.binaryisotherm import BinaryLangmuir


def main():
    # step 1: get H2S data
    df = pd.read_csv('../data_sets/CH4_H2S_MFI_binary_with_fugacity.csv')
    p_i, p_j, q_i, T = df['fugacity H2S [Pa]'], df['fugacity CH4 [Pa]'], df['Q H2S [mmol/g]'], df['T [K]']
    all_points = [i for i in range(len(q_i)) if q_i[i] > 0.]
    unary_points = [i for i in range(len(q_i)) if p_j[i] < 1e-12]

    # step 2: create models for h2s
    h2s_unary = LangmuirUnary(
        [p_i[i] for i in unary_points],
        [q_i[i] for i in unary_points],
        [T[i] for i in unary_points],
    )
    h2s_binary = BinaryLangmuir(
        [p_i[i] for i in all_points],
        [p_j[i] for i in all_points],
        [q_i[i] for i in all_points],
        [T[i] for i in all_points]
    )

    # step 3: get CH4 data
    df = pd.read_csv('../data_sets/CH4_H2S_MFI_binary_with_fugacity.csv')
    p_i, p_j, q_i, T = df['fugacity CH4 [Pa]'], df['fugacity H2S [Pa]'], df['Q CH4 [mmol/g]'], df['T [K]']
    # clean data to remove state points at which no H2S is adsorbed
    all_points = [i for i in range(len(q_i)) if q_i[i] > 0.]
    unary_points = [i for i in range(len(q_i)) if p_j[i] < 1e-12]

    # step 4: create models for CH4
    ch4_unary = LangmuirUnary(
        [p_i[i] for i in unary_points],
        [q_i[i] for i in unary_points],
        [T[i] for i in unary_points],
    )
    ch4_binary = BinaryLangmuir(
        [p_i[i] for i in all_points],
        [p_j[i] for i in all_points],
        [q_i[i] for i in all_points],
        [T[i] for i in all_points]
    )

    # step 5: solve unary models and assess quality of fit
    ch4_unary.solve()
    h2s_unary.solve()
    print('R2 H2S unary: ', h2s_unary.get_R2())
    print('SSE H2S unary: ', h2s_unary.get_objective())
    print('R2 CH4 unary: ', ch4_unary.get_R2())
    print('SSE CH4 unary: ', ch4_unary.get_objective())

    # step 6: plot unary models
    fig = plt.figure()
    fig, ax = h2s_unary.plot_comparison_dimensionless(fig=fig, color='red', marker='o', markerfacecolor='None', label='H2S unary')
    ch4_unary.plot_comparison_dimensionless(fig=fig, ax=ax, color='blue', marker='x', markerfacecolor='None', label='CH4 unary')

    # step 7: initialize binary variables from Langmuir combining rule
    h2s_binary.H_i_star = pyo.value(h2s_unary.H_i_star)
    h2s_binary.A_i = pyo.value(h2s_unary.A_i)
    h2s_binary.M_i_star = pyo.value(h2s_unary.M_i_star)
    h2s_binary.A_j = pyo.value(ch4_unary.A_i)
    h2s_binary.H_j_star = pyo.value(ch4_unary.H_i_star)
    ch4_binary.H_i_star = pyo.value(ch4_unary.H_i_star)
    ch4_binary.A_i = pyo.value(ch4_unary.A_i)
    ch4_binary.M_i_star = pyo.value(ch4_unary.M_i_star)
    ch4_binary.A_j = pyo.value(h2s_unary.A_i)
    ch4_binary.H_j_star = pyo.value(h2s_unary.H_i_star)

    # step 8: solve binary models and assess quality of fit
    h2s_binary.solve()
    ch4_binary.solve()
    print('R2 H2S binary: ', h2s_binary.get_R2())
    print('SSE H2S binary: ', h2s_binary.get_objective())
    print('R2 CH4 binary: ', ch4_binary.get_R2())
    print('SSE CH4 binary: ', ch4_binary.get_objective())

    # step 9: add binary models to plot
    h2s_binary.plot_comparison_dimensionless(fig=fig, ax=ax, color='purple', marker='d', markerfacecolor='None', label='H2S binary')
    ch4_binary.plot_comparison_dimensionless(fig=fig, ax=ax, color='cyan', marker='s', markerfacecolor='None', label='CH4 binary')
    ax.legend()
    fig.savefig('h2s_ch4_example.png')


if __name__ == '__main__':
    main()