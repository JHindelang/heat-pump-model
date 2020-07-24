import os
import pandas as pd
import numpy as np
import model_functions as mf
import CoolProp.CoolProp as CP
import CoolProp.HumidAirProp as HA

# change working directory
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

# directory input and output files
directory = os.path.join(dir_path, 'input_data')
directory_out = os.path.join(dir_path, 'output')

# read file names
files = pd.read_csv('input_filenames_simulation.csv', sep=';')

# assume temperature difference for air source heat pumps
delta_T_l = 5

for i, row_file in files.iterrows():
    filename_samples = os.path.join(directory,  files.filename_operating_points[i])
    # read samples
    samples = pd.read_csv(
        filename_samples,
        sep=';', header=[3],
        dtype={'m_s': np.float64, 'm_l': np.float64, 'q_source': np.float64,
               'q_load': np.float64})
    hp = pd.read_csv(filename_samples, sep=';', header=[0], nrows=2)
    hp_array = hp.loc[0, :].to_numpy()

    # list containing heat pump specifications
    # 0 - hp-type, 1 - compressor type, 2 - refrigerant,
    # 3 - source composition/ relative humidity , 4 - load composition
    hp_specification = hp.loc[0, "HP-Type":"pressure load"].to_numpy()

    # convert sample values
    samples.T_s = samples.T_s + 273.15
    samples.T_l = samples.T_l + 273.15
    samples.q_load = samples.q_load * 10 ** 3
    samples.q_source = samples.q_source * 10 ** 3
    samples.w_el = samples.w_el * 10 ** 3

    # convert source volume flow rate to mass flow rate
    if hp_array[10] == 'l/s':
        if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
            for k, row_samples in samples.iterrows():
                rho = 1 / HA.HAPropsSI(
                    'Vha', 'T', samples.T_s[k], 'P', hp_specification[5], 'R',
                    samples.comment[k])
                m_s = samples.m_s[k] * rho * 10 ** (-3)
                samples.at[k, 'm_s'] = m_s
            hp.set_value('flow type source', 0, 'kg/s')
        else:
            for k, row_samples in samples.iterrows():
                rho = CP.PropsSI('D', 'P', hp_specification[5], 'T',
                                 samples.T_s[k], hp_specification[3])
                m_s = samples.m_s[k] * rho * 10 ** (-3)
                samples.at[k, 'm_s'] = m_s
            hp.set_value('flow type source', 0, 'kg/s')
    elif hp_array[10] == 'm^3/h':
        if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
            for k, row_samples in samples.iterrows():
                rho = 1 / HA.HAPropsSI(
                    'Vha', 'T', samples.T_s[k], 'P', hp_specification[5], 'R',
                    samples.comment[k])
                m_s = samples.m_s[k] * rho / 3600
                samples.at[k, 'm_s'] = m_s  # m_s
            hp.set_value('flow type source', 0, 'kg/s')
        else:
            for k, row_samples in samples.iterrows():
                rho = CP.PropsSI('D', 'P', hp_specification[5], 'T',
                                 samples.T_s[k], hp_specification[3])
                m_s = samples.m_s[k] * rho / 3600
                samples.at[k, 'm_s'] = m_s
            hp.set_value('flow type source', 0, 'kg/s')

    # convert load volume flow rate to mass flow rate
    if hp_array[11] == 'l/s':  # load
        for k, row_samples in samples.iterrows():
            rho = CP.PropsSI('D', 'P', hp_specification[6], 'T',
                             samples.T_l[k], hp_specification[4])
            m_l = samples.m_l[k] * rho * 10 ** (-3)
            samples.at[k, 'm_l'] = m_l
        hp.set_value('flow type load', 0, 'kg/s')
    elif hp_array[11] == 'm^3/h':
        for k, row_samples in samples.iterrows():
            rho = CP.PropsSI('D', 'P', hp_specification[6], 'T',
                             samples.T_l[k], hp_specification[4])
            m_l = samples.m_l[k] * rho / 3600
            samples.at[k, 'm_l'] = m_l
        hp.set_value('flow type load', 0, 'kg/s')

    if hp_array[9] != 0:
        if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
            for k, row_samples in samples.iterrows():
                c_p = HA.HAPropsSI(
                    'cp_ha', 'T', samples.T_s[k], 'P', hp_specification[5],
                    'R', samples.comment[k])
                q_source = samples.m_s[k] * c_p * hp_array[9]
                samples.at[k, 'q_source'] = q_source

        else:
            for k, row_samples in samples.iterrows():
                c_p = CP.PropsSI('C', 'P', hp_specification[5], 'T',
                                 samples.T_s[k], hp_specification[3])
                q_source = samples.m_s[k] * c_p * hp_array[9]
                samples.at[k, 'q_source'] = q_source


    print(filename_samples)
    filename_parameter = os.path.join(directory, files.parameter_file[i])
    parameters = pd.read_csv(
        filename_parameter, header=[0], skiprows=[0], index_col=[0], sep=';')
    params = {'DE': parameters.loc['DE', 'ua_l':'delta_p'].to_numpy()}

    results = samples.copy()
    results['m_ref'] = np.nan
    results['iter'] = np.nan
    results['q_load_calc'] = np.nan
    results['load_dev'] = np.nan

    for k, row_samples in samples.iterrows():

        [m_ref, q_source, comp_p, q_load_calc, iter] = \
            mf.compute_hp_performance(
                params['DE'], samples.T_s[k], samples.T_l[k], samples.q_load[k],
                samples.m_s[k], samples.m_l[k], samples.comment[k],
                hp_specification)
        results.at[k, 'w_el'] = comp_p
        results.at[k, 'cop'] = samples.q_load[k] / comp_p
        results.at[k, 'q_source'] = q_source
        results.at[k, 'm_ref'] = m_ref
        results.at[k, 'iter'] = iter
        results.at[k, 'q_load_calc'] = q_load_calc
        results.at[k, 'load_dev'] = (q_load_calc - samples.q_load[k]) / samples.q_load[k]

    results.T_s = results.T_s - 273.15
    results.T_l = results.T_l - 273.15
    results.q_load = results.q_load / (10 ** 3)
    results.q_source = results.q_source / (10 ** 3)
    results.w_el = results.w_el / (10 ** 3)
    results.q_load_calc = results.q_load_calc / (10 ** 3)

    filename_evaluation = os.path.join(
        directory_out, 'sim_' + files.filename_operating_points[i])
    hp.to_csv(filename_evaluation, header=True, index=False, sep=';')
    results.to_csv(filename_evaluation, index=False, header=True, sep=';', mode='a')
