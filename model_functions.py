import numpy as np
import pandas as pd
import CoolProp
import CoolProp.CoolProp as CP
import CoolProp.HumidAirProp as HA



def sqsumerr(params, samples, hp_specification):
    # computes the sum of square errors for a given parameter set and given
    # measurements
    sse = 0

    for i, row in samples.iterrows():
        # compute the heat pump performance
        [q_load_calc_i, w_comp_calc_i, q_source_calc_i] = \
            compute_hp_param_estimation(
                params, samples.m_s[i], samples.m_l[i], samples.T_s[i],
                samples.T_l[i], samples.q_load[i], samples.q_source[i],
                samples.comment[i], hp_specification)
        # compute the sse
        sse = sse + ((samples.w_el[i] - w_comp_calc_i) / samples.w_el[i]) ** 2 + \
              ((samples.q_load[i] - q_load_calc_i) / samples.q_load[i]) ** 2 + \
              ((samples.q_source[i] - q_source_calc_i) / samples.q_source[i]) ** 2

    return sse


def sqsumerr_check(params, samples, hp_specification):
    # computes the sum of square errors for a given parameter set and given
    # measurements
    sse = 0

    for i, row in samples.iterrows():
        # compute the heat pump performance
        [m_ref, q_source_calc_i, w_comp_calc_i, q_load_calc_i, iter] \
            = compute_hp_performance(
            params, samples.T_s[i], samples.T_l[i], samples.q_load[i],
            samples.m_s[i], samples.m_l[i], samples.comment[i],
            hp_specification)
        # compute the sse
        sse = sse + ((samples.w_el[i] - w_comp_calc_i) / samples.w_el[i]) ** 2 + \
              ((samples.q_load[i] - q_load_calc_i) / samples.q_load[i]) ** 2 + \
              ((samples.q_source[i] - q_source_calc_i) / samples.q_source[i]) ** 2

    return sse


def thermal_eff(ua, m_dot, c_p):
    # compute thermal efficiency of counter flow heat exchanger for given UA
    # ua, mass flow rate m_dot and heat capacity c_p
    epsilon = 1 - np.exp(-ua / (c_p * m_dot))
    return epsilon


def thermal_eff_ha(ua, m_dot, c_p, C):
    # compute thermal efficiency of heat exchanger using humid air for given
    # UA, mass flow rate m_dot, heat capacity c_p
    # and heaf flow ratio C
    NTU = -ua / (c_p * m_dot)
    epsilon = 1 - np.exp(((NTU**0.22)*np.exp(-C*NTU**0.78 - 1))/C)
    return epsilon


def sat_prop(fluid, temperature):
    # return the saturation pressure, saturated liquid and saturated vapor
    # enthalpy for fluid at temperature from CoolProp
    pressure = CP.PropsSI('P', 'T', temperature, 'Q', 1, fluid)
    h_liq = CP.PropsSI('H', 'P', pressure, 'Q', 0, fluid)
    h_vap = CP.PropsSI('H', 'P', pressure, 'Q', 1, fluid)
    return [pressure, h_liq, h_vap]


def enth_sh(fluid, temperature, pressure):
    # return enthalpy of superheated vapor at compressor inlet for fluid at
    # pressure and temperature from CoolProp
    phase = CP.PhaseSI('P', pressure, 'T', temperature, fluid)
    if phase == 'gas':
        h_sh = CP.PropsSI('H', 'T', temperature, 'P', pressure, fluid)
    else:
        h_sh = CP.PropsSI('H', 'P', pressure, 'Q', 1, fluid)
        # print('Warning: Values set to saturation (h_sh)')

    return h_sh


def volume_suc(fluid, temperature, pressure):
    # return specific volume for fluid at temperature and pressure from
    # CoolProp
    T_min = CP.PropsSI('T', 'P', pressure, 'Q', 1, fluid)
    if temperature <= T_min +0.1 :
        # if temperature below saturation temperature at pressure, return
        # specific volume of saturated vapor
        v_suc = 1 / CP.PropsSI('D', 'P', pressure, 'Q', 1, fluid)
        # print('Warning: Values set to saturation (volume_suc)')
    else:
        v_suc = 1 / CP.PropsSI('D', 'T', temperature, 'P', pressure, fluid)
    return v_suc


def gamma(temperature, pressure, fluid):
    # return isentropic coefficient for fluid at temperature and pressure
    # from CoolProp
    T_min = CP.PropsSI('T', 'P', pressure, 'Q', 1, fluid)
    if temperature <= T_min +0.1:
        # if temperature below saturation temperature at pressure, return
        # isentropic coefficient of saturated vapor
        gamma_calc = CP.PropsSI(
            'ISENTROPIC_EXPANSION_COEFFICIENT', 'P', pressure, 'Q', 1, fluid)
        # print('Warning: Values set to saturation (gamma)')
    else:
        gamma_calc = CP.PropsSI(
            'ISENTROPIC_EXPANSION_COEFFICIENT', 'T', temperature, 'P',
            pressure, fluid)
    return gamma_calc


def compressor_inlet_state(temperature, pressure, fluid):
    # return isentropic coefficient ad specific volume for fluid at temperature
    # and pressure from CoolProp
    phase = CP.PhaseSI('P', pressure, 'T', temperature, fluid)

    if phase == 'gas':
        v_suc = 1 / CP.PropsSI('D', 'T', temperature, 'P', pressure, fluid)
        gamma_calc = CP.PropsSI('ISENTROPIC_EXPANSION_COEFFICIENT', 'T',
                                temperature, 'P', pressure, fluid)
    else:  #if temperature below saturation temperature at pressure, assume
        # saturated state
        v_suc = 1 / CP.PropsSI('D', 'P', pressure,'Q', 1, fluid)
        gamma_calc = CP.PropsSI('ISENTROPIC_EXPANSION_COEFFICIENT', 'P',
                                pressure, 'Q', 1, fluid)
        # print('Warning: Values set to saturation (compressor inlet state)')
    return [v_suc, gamma_calc]

def compressor_inlet_state_reci(enthalpy, pressure, fluid):
    # return isentropic coefficient and specific volume for fluid at enthalpy
    # and pressure from CoolProp
    phase = CP.PhaseSI('P', pressure, 'H', enthalpy, fluid)

    if phase == 'gas':
        v_suc = 1 / CP.PropsSI('D', 'H', enthalpy, 'P', pressure, fluid)
        gamma_calc = CP.PropsSI('ISENTROPIC_EXPANSION_COEFFICIENT', 'H',
                                enthalpy, 'P', pressure, fluid)
    else:  # if enthalpy below saturated vapor enthalpy at pressure, assume
        # saturated state
        v_suc = 1 / CP.PropsSI('D', 'P', pressure,'Q', 1, fluid)
        gamma_calc = CP.PropsSI('ISENTROPIC_EXPANSION_COEFFICIENT', 'P',
                                pressure, 'Q', 1, fluid)
        # print('Warning: Values set to saturation (compressor inlet state reci)')
    return [v_suc, gamma_calc]


def reciprocating_compressor(displacement, clearance, pressure_suc,
                             pressure_dis, gamma_comp, v_suc):
    # compute refrigerant mass flow rate m_ref and isentropic compression work
    # for reciprocating compressor
    pressure_rat = pressure_dis / pressure_suc
    m_ref = (displacement / v_suc) * (1 + clearance - clearance *
                                      (pressure_rat ** (1 / gamma_comp)))
    compr_p_is = (gamma_comp / (gamma_comp - 1)) * m_ref * pressure_suc * v_suc * (
                (pressure_rat ** ((gamma_comp - 1) / gamma_comp) - 1))

    return [m_ref, compr_p_is]


def scroll_compressor(pocket_volume, leakage_fact, pressure_suc, pressure_dis,
                      gamma_comp, v_suc, volume_ratio):
    # compute refrigerant mass flow rate m_ref and technical work for scroll
    # compressor
    pressure_rat = pressure_dis / pressure_suc
    m_ref = pocket_volume / v_suc
    m_leakage = leakage_fact * pressure_rat
    if m_leakage >= m_ref:
        # if leakage mass flow rate is to high, set m_ref approximately to zero
        m_ref = 0.00000001
    else:
        m_ref = m_ref-m_leakage

    compr_p_is = (pressure_suc * pocket_volume)/(gamma_comp - 1) * \
                 (volume_ratio ** (gamma_comp - 1) - gamma_comp) + \
                 (pressure_suc * pocket_volume * pressure_rat) / volume_ratio

    return [m_ref, compr_p_is]


def rotary_compressor(volume, pressure_suc, pressure_dis, gamma_comp, v_suc):
    # compute refrigerant mass flow rate m_ref and isentropic compression work
    # for rotary compressor
    pressure_rat = pressure_dis / pressure_suc
    m_ref = volume / v_suc
    compr_p_is = (gamma_comp / (gamma_comp - 1)) * m_ref * pressure_suc * v_suc * (
                (pressure_rat ** ((gamma_comp - 1) / gamma_comp) - 1))
    return [m_ref, compr_p_is]


def compute_hp_param_estimation(params, m_s, m_l, T_s, T_l,
                                q_load, q_source, humidity, hp_specification):
    refrigerant = hp_specification[2]

    T_min_ref = CP.PropsSI('TMIN', refrigerant)
    T_max_ref = CP.PropsSI('TCRIT', refrigerant)
    p_min_ref = CP.PropsSI('PMIN', refrigerant)
    p_max_ref = CP.PropsSI('PCRIT', refrigerant)
    source_composition = hp_specification[3]
    load_composition = hp_specification[4]
    # constants
    cp_l = CP.PropsSI('CP0MASS', 'T', T_l, 'P', hp_specification[6],
                      load_composition)
    if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
        # humid air properties
        cp_s = HA.HAPropsSI('cp_ha', 'T', T_s, 'P', hp_specification[5], 'R',
                            humidity)
        epsilon_s = thermal_eff(params[1], m_s, cp_s)
    else:
        cp_s = CP.PropsSI('C', 'T', T_s, 'P', hp_specification[5],
                          source_composition)
        epsilon_s = thermal_eff(params[1], m_s, cp_s)

    epsilon_l = thermal_eff(params[0], m_l, cp_l)

    [m_ref, q_source, q_heating, comp_p] = heat_pump_cycle(
        params, T_s, T_l, q_load, q_source, m_s, m_l, hp_specification,
        [cp_l, cp_s], [epsilon_l, epsilon_s],[T_min_ref, T_max_ref, p_min_ref,
                                              p_max_ref])

    return [q_heating, comp_p, q_source]


def compute_hp_performance(params, T_s, T_l, q_load, m_s, m_l, humidity,
                           hp_specification):
    # constants
    iter = 0
    maxiter = 1000
    tol = 5 * 10 ** (-2)
    q_load_calc = 0
    refrigerant = hp_specification[2]
    source_composition = hp_specification[3]
    load_composition = hp_specification[4]
    T_min_ref = CP.PropsSI('TMIN', refrigerant)
    T_max_ref = CP.PropsSI('TCRIT', refrigerant)
    p_min_ref = CP.PropsSI('PMIN', refrigerant)
    p_max_ref = CP.PropsSI('PCRIT', refrigerant)

    cp_l = CP.PropsSI('CP0MASS', 'T', T_l, 'P', hp_specification[6],
                      load_composition)
    epsilon_l = thermal_eff(params[0], m_l, cp_l)

    if hp_specification[0] == 'a-a' or hp_specification[0] == 'a-w':
        # humid air properties
        cp_s = HA.HAPropsSI('cp_ha', 'T', T_s, 'P', hp_specification[5], 'R',
                            humidity)
        epsilon_s = thermal_eff(params[1], m_s, cp_s)
    else:
        cp_s = CP.PropsSI('C', 'T', T_s, 'P', hp_specification[5],
                          source_composition)
        epsilon_s = thermal_eff(params[1], m_s, cp_s)
    # initial guess for q_source
    q_source = q_load * 0.8

    while abs((q_load - q_load_calc) / q_load) > tol and iter < maxiter:
        [m_ref, q_source, q_load_calc, comp_p] = heat_pump_cycle(params, T_s,
         T_l, q_load, q_source, m_s, m_l, hp_specification, [cp_l, cp_s],
         [epsilon_l, epsilon_s], [T_min_ref, T_max_ref, p_min_ref, p_max_ref])
        iter = iter + 1

    return [m_ref, q_source, comp_p, q_load_calc, iter]


def heat_pump_cycle(params, T_s, T_l, q_load,q_source, m_s, m_l,
                    hp_specification, cp, epsilon, thermo_lim):

    refrigerant = hp_specification[2]

    # condenser and evaporator temperatures
    t_l_in = T_l - q_load / (m_l * cp[0])
    t_cond = t_l_in + q_load / (epsilon[0] * cp[0] * m_l)
    t_evap = T_s - q_source / (epsilon[1] * cp[1] * m_s)
    # check temperature range
    if t_cond < thermo_lim[0]:
        t_cond = thermo_lim[0]
    elif t_cond > 0.98 * thermo_lim[1]:
        t_cond = 0.98 * thermo_lim[1]

    if t_evap < thermo_lim[0]:
        t_evap = thermo_lim[0]
    elif t_evap > 0.98 * thermo_lim[1]:
        t_evap = 0.98 * thermo_lim[1]

    # condenser and evaporator pressure
    [p_evap, h_ev_liq, h_ev_vap] = sat_prop(refrigerant, t_evap)
    [p_cond, h_cond_liq, h_cond_vap] = sat_prop(refrigerant, t_cond)

    # compressor inlet temperature
    t_comp_in = t_evap + params[2]
    if t_comp_in < thermo_lim[0]:
        t_comp_in = thermo_lim[0]
    elif t_cond > thermo_lim[1]:
        t_comp_in = thermo_lim[1]

    h_comp_in = enth_sh(refrigerant, t_comp_in, p_evap)

    if hp_specification[1] == 'reciprocating':
        p_suc = p_evap - params[7]
        p_dis = p_cond + params[7]

        if p_suc < thermo_lim[2]:
            p_suc = thermo_lim[2]
        elif p_suc > thermo_lim[3]:
            p_suc = thermo_lim[3]

        if p_dis < thermo_lim[2]:
            p_dis = thermo_lim[2]
        elif p_dis > thermo_lim[3]:
            p_dis = thermo_lim[3]

        [v_suc, gamma_comp] = compressor_inlet_state_reci(
            h_comp_in, p_suc, refrigerant)
        [m_ref, compr_p_is] = reciprocating_compressor(
            params[3], params[6], p_suc, p_dis, gamma_comp, v_suc)
    elif hp_specification[1] == 'scroll':
        [v_suc, gamma_comp] = compressor_inlet_state(t_comp_in, p_evap,
                                                     refrigerant)
        [m_ref, compr_p_is] = scroll_compressor(
            params[3], params[6], p_evap,p_cond, gamma_comp, v_suc, params[7])
    elif hp_specification[1] == 'rotary':
        p_dis = p_cond + params[6]

        if p_dis < thermo_lim[2]:
            p_dis = thermo_lim[2]
        elif p_dis > thermo_lim[3]:
            p_dis = thermo_lim[3]

        [v_suc, gamma_comp] = compressor_inlet_state(t_comp_in, p_evap,
                                                     refrigerant)
        [m_ref, compr_p_is] = rotary_compressor(params[3], p_evap, p_dis,
                                                gamma_comp, v_suc)

    compr_p = (compr_p_is / params[4]) + params[5]

    # assume no subcooling isenthalpic expansion
    q_source = (m_ref * (h_ev_vap - h_cond_liq))

    q_load_calc = (m_ref * (h_ev_vap - h_cond_liq) + compr_p_is)

    return[m_ref, q_source, q_load_calc, compr_p]
