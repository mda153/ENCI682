# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:32:55 2022

@author: mda153
"""

import numpy as np
import o3seespy as o3
import matplotlib.pyplot as plt


def add_unconf_concrete_to_plot(ax, label='Unconf. concrete', c='b'):
    osi = o3.OpenSeesInstance(ndm=1, ndf=1, state=3)

    f_c28 = 17 * 1E3  # [kPa] 28 day compressive strength
    e_mod_conc = (3320 * np.sqrt(f_c28) + 6900) * 1E3  # [kPa] NZS3101 - 2006
    # f_conc_tens = (0.38 * np.sqrt(f_c28))  # [kPa] NZS3101 CL 5.2.4, direct tensile strength concrete
    f_conc_tens = 1.4 * 1E3  # [kPa] Defined as this for this section in Henderson Thesis
    eps_conc_tens = f_conc_tens / e_mod_conc  # Park and Paulay 1975 assume linear
    eps_unconf_conc_max = 0.004  # Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)

    eps_conc_crush = 0.003  # Concrte strain at crushing - Koopaee (2015) from NZS 3101?? taken from Fig 2.1 fc=20MPa curve

    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-f_c28, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush,
                                                  ec=e_mod_conc, fct=f_conc_tens,
                                                  et=eps_conc_tens)  # unconfined concrete properties
    # conc_unconf = o3.uniaxial_material.Concrete01(osi, fpc=-f_c28 * 0.8, epsc0=-eps_unconf_conc_max, fpcu=0.0,
    #                                               eps_u=-eps_conc_crush)
    peak_disps = np.array([-0.0001, -0.00005, -0.001])
    peak_disps = np.array([-0.005])
    rec_d = 0.00001
    disps = []
    init = 0
    for dd in range(len(peak_disps)):
        sgn = np.sign(peak_disps[dd] - init)
        disps += list(np.arange(init, peak_disps[dd], sgn * rec_d))
        init = peak_disps[dd]
    disp, react = o3.tools.run_uniaxial_disp_driver(osi, conc_unconf, disps)
    ax.plot(disp, react, label=label, c=c)



def add_conf_concrete_to_plot(ax, label='Conf. concrete', c='r'):
    osi = o3.OpenSeesInstance(ndm=1, ndf=1, state=3)

    f_c28 = 17 * 1E3  # [kPa] 28 day compressive strength
    e_mod_conc = (3320 * np.sqrt(f_c28) + 6900) * 1E3  # [kPa] NZS3101 - 2006
    # f_conc_tens = (0.38 * np.sqrt(f_c28))  # [kPa] NZS3101 CL 5.2.4, direct tensile strength concrete
    f_conc_tens = 1.4 * 1E3  # [kPa] Defined as this for this section in Henderson Thesis
    eps_conc_tens = f_conc_tens / e_mod_conc  # Park and Paulay 1975 assume linear
    eps_conf_conc_max = 0.005  # Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)
    eps_conc_crush_conf = 0.02

    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-f_c28, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf,
                                                ec=e_mod_conc, fct=f_conc_tens,
                                                et=eps_conc_tens)  # Confined concrete paramters
    peak_disps = np.array([-0.0001, -0.00005, -0.001])
    peak_disps = np.array([-0.005])
    rec_d = 0.00001
    disps = []
    init = 0
    for dd in range(len(peak_disps)):
        sgn = np.sign(peak_disps[dd] - init)
        disps += list(np.arange(init, peak_disps[dd], sgn * rec_d))
        init = peak_disps[dd]
    disp, react = o3.tools.run_uniaxial_disp_driver(osi, conc_conf, disps)
    ax.plot(disp, react, label=label, c=c)


def add_rebar_to_plot(ax, label="Rebar"):
    osi = o3.OpenSeesInstance(ndm=1, ndf=1, state=3)
    fy = 300 * 1E3  # [kPa] steel yield strength
    Es = 200 * 1E6  # [kPa] youngs modulus - steel
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.001)  # Reinforcing steel
    peak_disps = np.array([-0.0001, -0.00005, -0.001])
    peak_disps = np.array([-0.005])
    rec_d = 0.00001
    disps = []
    init = 0
    for dd in range(len(peak_disps)):
        sgn = np.sign(peak_disps[dd] - init)
        disps += list(np.arange(init, peak_disps[dd], sgn * rec_d))
        init = peak_disps[dd]
    disp, react = o3.tools.run_uniaxial_disp_driver(osi, rebar, disps)
    ax.plot(disp, react, label=label)


def add_mander_model_conf_concrete_to_plot(ax):
    # TODO: add the mander model equations
    pass

def add_mander_model_unconf_concrete_to_plot(ax, label="Mander"):
    # TODO: add the mander model equations
    
    # Mander unconfined concrete model
    
    fcu = []
    dels = 0.0001 #Delta strain for default material models (0.0001)
    espall = 0.0064 #Maximum uncon. conc. strain (usually 0.0064)
     
    fpc = 17 #Conc compressive strength (MPa)
    eco = 0.002 #Unconfined concrete strain (ususally 0.002 for normal weight or 0.004 for lightweight)
    Ec = 5000*np.sqrt(fpc) #Conc modulus of elasticity (MPa) in Mander script
    
    ec = ([0, dels, espall])
    Escu = fpc / eco
    ru = Ec / (Ec-Escu)
    xu = ec / (eco)
    
    for i in [1, len(ec)]:
        if ec(i)<2*eco:
            fcu(i) = fpc*xu(i)*ru/(ru-1 + xu(i)^ru)
    
        if (ec(i)>=2*eco and ec(i)<=espall):
            fcu(i) = fpc*(2*ru/(ru-1+2^ru))*(1-(ec(i) - 2*eco)/(espall - 2*eco))

        
        if ec(i) > espall:
            fcu(i) = 0
            
    ax.plot(ec, fcu, label=label)



def create():
    bf, ax = plt.subplots(nrows=2, sharex='col')
    add_unconf_concrete_to_plot(ax[0])
    add_conf_concrete_to_plot(ax[0])
    add_mander_model_conf_concrete_to_plot(ax[0])
    add_mander_model_unconf_concrete_to_plot(ax[0])

    add_rebar_to_plot(ax[1])
    ax[0].set_ylabel('Stress [Pa]')
    ax[1].set_ylabel('Stress [Pa]')
    ax[-1].set_xlabel('Strain')
    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    create()

