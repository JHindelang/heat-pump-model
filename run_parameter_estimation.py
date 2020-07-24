if __name__ == '__main__':
    import os
    import pandas as pd
    import numpy as np
    from scipy import optimize
    import model_functions as mf
    import CoolProp.CoolProp as CP
    import CoolProp.HumidAirProp as HA
    import time

    seed_DE = [5, 15, 25, 31, 37, 49, 55, 86, 95, 98]

    # initial parameters
    class ParameterConversion:
        pass

    # change working directory
    dir_path = os.path.dirname(os.path.realpath(__file__))
    os.chdir(dir_path)

    # directory input and output files
    directory = os.path.join(dir_path, 'input_data')
    directory_out = os.path.join(dir_path, 'output')

    # read file with list of heat pumps
    try:
        files = pd.read_csv('input_filenames_estimation.csv')
    except FileNotFoundError:
        print("File not found")

    for i, row_file in files.iterrows():
        # compose directory of sample file
        filename_samples = os.path.join(directory, files.filename[i])

        # read sample file
        try:
            samples = pd.read_csv(filename_samples, sep=';', header=[3], dtype={'m_s': np.float64, 'm_l': np.float64,
                                                                            'q_source': np.float64, 'q_load': np.float64})
        except FileNotFoundError:
            print('file ', filename_samples, ' not found')
            continue

        print(filename_samples)

        # read file header
        hp = pd.read_csv(filename_samples, sep=';', header=[0], nrows=1)

        # copy to array
        hp_array = hp.loc[0,:].to_numpy()

        # extract heat pump specifications
        # list containing heat pump specifications
        # 0 - hp-type, 1 - compressor type, 2 - refrigerant,
        # 3 - source composition/ relative humidity , 4 - load composition
        hp_specification = hp.loc[0, "HP-Type":"pressure load"].to_numpy()

        # convert sample values
        samples.T_s = samples.T_s + 273.15  # [K]
        samples.T_l = samples.T_l + 273.15  # [K]
        samples.q_load = samples.q_load * 10 ** 3 #[W]
        samples.q_source = samples.q_source * 10 ** 3 #[W]
        samples.w_el = samples.w_el * 10 ** 3 #[W]

        # convert source volume flow rate to mass flow rate
        if hp_array[10] == 'l/s':
            if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
                for k, row_samples in samples.iterrows():
                    rho = 1 / HA.HAPropsSI('Vha', 'T', samples.T_s[k], 'P', hp_specification[5], 'R', samples.comment[k])
                    m_s = samples.m_s[k] * rho * 10 ** (-3)
                    samples.at[k, 'm_s'] = 5.78 #m_s
            else:
                for k, row_samples in samples.iterrows():
                    rho = CP.PropsSI('D', 'P', hp_specification[5], 'T', samples.T_s[k], hp_specification[3])
                    m_s = samples.m_s[k] * rho * 10 **(-3)
                    samples.at[k, 'm_s'] = m_s
        elif hp_array[10] == 'm^3/h':
            if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
                for k, row_samples in samples.iterrows():
                    rho = 1 / HA.HAPropsSI('Vha', 'T', samples.T_s[k], 'P', hp_specification[5], 'R', samples.comment[k])
                    m_s = samples.m_s[k] * rho / 3600
                    samples.at[k, 'm_s'] = m_s
            else:
                for k, row_samples in samples.iterrows():
                    rho = CP.PropsSI('D', 'P', hp_specification[5], 'T', samples.T_s[k], hp_specification[3])
                    m_s = samples.m_s[k] * rho / 3600
                    samples.at[k, 'm_s'] = m_s

        # convert load volume flow rate to mass flow rate
        if hp_array[11] == 'l/s': # load
            for k, row_samples in samples.iterrows():
                rho = CP.PropsSI('D', 'P', hp_specification[6], 'T', samples.T_l[k], hp_specification[4])
                m_l = samples.m_l[k] * rho * 10 **(-3)
                samples.at[k, 'm_l'] = m_l
        elif hp_array[11] == 'm^3/h':
            for k, row_samples in samples.iterrows():
                rho = CP.PropsSI('D', 'P', hp_specification[6], 'T', samples.T_l[k], hp_specification[4])
                m_l = samples.m_l[k] * rho / 3600
                samples.at[k, 'm_l'] = m_l

        # compute source heat transfer rate if required
        if hp_array[9] != 0:
            if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
                for k, row_samples in samples.iterrows():
                    c_p = HA.HAPropsSI('cp_ha', 'T', samples.T_s[k], 'P', hp_specification[5], 'R', samples.comment[k])
                    q_source = samples.m_s[k] * c_p * hp_array[11]
                    samples.at[k, 'q_source'] = q_source

            else:
                for k, row_samples in samples.iterrows():
                    c_p = CP.PropsSI('C', 'P', hp_specification[5], 'T', samples.T_s[k], hp_specification[3])
                    q_source = samples.m_s[k] * c_p * hp_array[11]
                    samples.at[k, 'q_source'] = q_source

        # specify bounds for parameters according to compressor type
        if hp_specification[1] == 'reciprocating':
            lb = [0.001, 0.001, 0, 0, 0.001, 0, 0, 0]
            ub = [1000000, 1000000, 50, 1, 1, 1000, 1, 100000]
            A = np.eye(8)
        elif hp_specification[1] == 'scroll':
            lb = [0.001, 0.001, 0, 0.0005, 0.00001, 0, 0.00001, 0.5]
            ub = [15000, 15000, 20, 10, 1, 1000, 1, 3.5]
            A = np.eye(7)
        elif hp_specification[1] == 'rotary':
            lb = [0, 0, 0, 0, 0, 0, 0]
            ub = [np.inf, np.inf, 20, np.inf, 10, np.inf, 10]
            A = np.eye(8)

        # convert constraints to list
        bounds = list(zip(lb, ub))

        filename_params = os.path.join(directory_out, 'param_' + files.filename[i])
        hp.to_csv(filename_params, header=None, index=False, sep=';')

        # perform optimization
        results = dict()
        n = 0
        for rand in seed_DE:
            print('Start DE')
            start_time = time.time()
            results_DE = optimize.differential_evolution(
                mf.sqsumerr,
                bounds,
                args=(samples, hp_specification),
                mutation=(0.5, 1.5), recombination=0.8, popsize=15,
                seed=rand, updating='deferred', workers=-1, strategy='best1bin'
            )
            time_DE = (time.time() - start_time)/60
            print('execution time: ', time_DE, ' min')
            print('differential evolution completed')
            print(results_DE)

            if hp_specification[1] == 'rotary':
                col_names = ['ua_l', 'ua_s', 'd_T_sh', 'V_PD', 'eta', 'P_loss', 'delta_p', 'fun', 'message', 'nit',
                             'success', 'execution_time']
            if hp_specification[1] == 'scroll':
                col_names = ['ua_l', 'ua_s', 'd_T_sh', 'V_s', 'eta', 'P_loss', 'L', 'nu', 'fun', 'message', 'nit',
                             'success', 'execution_time']
            if hp_specification[1] == 'reciprocating':
                col_names = ['ua_l', 'ua_s', 'd_T_sh', 'V_PD', 'eta', 'P_loss', 'C', 'delta_P', 'fun', 'message', 'nit',
                             'success', 'execution_time']

            df = pd.DataFrame(
                [np.append(results_DE.x,
                 [results_DE.fun, results_DE.message, results_DE.nit,
                  results_DE.success, time_DE])],
                columns=col_names, index=['DE'])

            if n == 0:
                df.to_csv(filename_params, index=True, header=True, sep=';', mode='a')
            else:
                df.to_csv(filename_params, index=True, header=False, sep=';', mode='a')
            n += 1



