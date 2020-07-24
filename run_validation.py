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
files = pd.read_csv('input_filenames_validation.csv', sep=';')

for i, row_file in files.iterrows():
    filename_samples = os.path.join(directory, files.filename_samples[i])

    samples = pd.read_csv(
        filename_samples, sep=';', header=[3],
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
                    'Vha', 'T', samples.T_s[k], 'P', hp_specification[5],
                    'R', samples.comment[k])
                m_s = samples.m_s[k] * rho * 10 ** (-3)
                samples.at[k, 'm_s'] = m_s
            hp.at[0, 'flow type source'] = 'kg/s'
        else:
            for k, row_samples in samples.iterrows():
                rho = CP.PropsSI(
                    'D', 'P', hp_specification[5], 'T', samples.T_s[k],
                    hp_specification[3])
                m_s = samples.m_s[k] * rho * 10 ** (-3)
                samples.at[k, 'm_s'] = m_s
            hp.at[0, 'flow type source'] = 'kg/s'
    elif hp_array[10] == 'm^3/h':
        if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
            for k, row_samples in samples.iterrows():
                rho = 1 / HA.HAPropsSI(
                    'Vha', 'T', samples.T_s[k], 'P', hp_specification[5],
                    'R', samples.comment[k])
                m_s = samples.m_s[k] * rho / 3600
                samples.at[k, 'm_s'] = m_s  # m_s
            hp.at[0, 'flow type source'] = 'kg/s'
        else:
            for k, row_samples in samples.iterrows():
                rho = CP.PropsSI('D', 'P', hp_specification[5], 'T',
                                 samples.T_s[k], hp_specification[3])
                m_s = samples.m_s[k] * rho / 3600
                samples.at[k, 'm_s'] = m_s
            hp.at[0, 'flow type source'] = 'kg/s'

    # convert load volume flow rate to mass flow rate
    if hp_array[11] == 'l/s':  # load
        for k, row_samples in samples.iterrows():
            rho = CP.PropsSI('D', 'P', hp_specification[6], 'T',
                             samples.T_l[k], hp_specification[4])
            m_l = samples.m_l[k] * rho * 10 ** (-3)
            samples.at[k, 'm_l'] = m_l
        hp.at[0, 'flow type load'] = 'kg/s'
    elif hp_array[11] == 'm^3/h':
        for k, row_samples in samples.iterrows():
            rho = CP.PropsSI('D', 'P', hp_specification[6], 'T',
                             samples.T_l[k], hp_specification[4])
            m_l = samples.m_l[k] * rho / 3600
            samples.at[k, 'm_l'] = m_l
        hp.at[0, 'flow type load'] = 'kg/s'

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
                c_p = CP.PropsSI('C', 'P', hp_specification[5],
                                 'T', samples.T_s[k], hp_specification[3])
                q_source = samples.m_s[k] * c_p * hp_array[9]
                samples.at[k, 'q_source'] = q_source

    print(filename_samples)
    filename_parameter = os.path.join(directory, files.filename_parameter[i])
    parameters = pd.read_csv(
        filename_parameter, header=[0], skiprows=[0], index_col=[0], sep=';')
    params = {'DE': parameters.loc['DE', 'ua_l':'delta_p'].to_numpy()}

    results = samples.copy()
    results['w_el_comp'] = np.nan
    results['q_source_comp'] = np.nan
    results['m_ref'] = np.nan
    results['iter'] = np.nan
    results['q_load_calc'] = np.nan
    results['load_dev'] = np.nan

    sse_w_el = 0
    sse_q_s = 0
    sse_check = 0

    for k, row in samples.iterrows():
        [m_ref, q_source, comp_p, q_load_calc, iter] = \
            mf.compute_hp_performance(
                params['DE'], samples.T_s[k], samples.T_l[k], samples.q_load[k],
                samples.m_s[k], samples.m_l[k], samples.comment[k], hp_specification)
        results.at[k, 'w_el_comp'] = comp_p
        results.at[k, 'q_source_comp'] = q_source
        results.at[k, 'm_ref'] = m_ref
        results.at[k, 'iter'] = iter
        results.at[k, 'q_load_calc'] = q_load_calc
        results.at[k, 'load_dev'] = (q_load_calc - samples.q_load[k]) / samples.q_load[k]

        print('\n', m_ref, q_source, comp_p, iter)
        print(samples.q_source[k], samples.w_el[k])
        print('rel. deviation heat source: ',
              (q_source/samples.q_source[k] - 1)*100, '%', ' rel. deviation power: ',
              (comp_p/samples.w_el[k] - 1)*100, '%')

        sse_w_el = sse_w_el + (comp_p-samples.w_el[k])**2
        sse_q_s= sse_q_s + (q_source - samples.q_source[k]) ** 2

    rms_w_el = np.sqrt(sse_w_el/k)
    rms_q_s = np.sqrt(sse_q_s/k)
    sse = mf.sqsumerr_check(params['DE'], samples, hp_specification)
    hp.loc[hp.index[0], 'rms_w_el'] = rms_w_el
    hp.loc[hp.index[0], 'rms_q_s'] = rms_q_s
    hp.loc[hp.index[0], 'sse'] = sse

    results.T_s = results.T_s - 273.15
    results.T_l = results.T_l - 273.15
    results.q_load = results.q_load / (10 ** 3)
    results.q_source = results.q_source / (10 ** 3)
    results.w_el = results.w_el / (10 ** 3)
    results.q_source_comp = results.q_source_comp / (10 ** 3)
    results.w_el_comp = results.w_el_comp / (10 ** 3)
    results.q_load_calc = results.q_load_calc / (10 ** 3)

    filename_validation = os.path.join(
        directory_out, 'validation_' + files.filename_samples[i])
    hp.to_csv(filename_validation, header=True, index=False, sep=';')
    results.to_csv(filename_validation, index=False, header=True, sep=';', mode='a')
