# -*- coding: utf-8 -*-
"""
Created on Wed May  4 16:32:55 2022

@author: mda153
"""

import numpy as np
import o3seespy as o3
import matplotlib.pyplot as plt


def add_unconf_concrete_to_plot(ax, fpc, espall, eps_unconf_conc_max, eps_conc_crush_unconf, Ec, label='O3SeesPy Unconfined Concrete (Concrete 04)'):
    osi = o3.OpenSeesInstance(ndm=1, ndf=1, state=3)
    
    from bwplot import cbox

    #f_c28 = 17 * 1E3  # [kPa] 28 day compressive strength
    #e_mod_conc = (3320 * np.sqrt(f_c28) + 6900) * 1E3  # [kPa] NZS3101 - 2006
    # f_conc_tens = (0.38 * np.sqrt(f_c28))  # [kPa] NZS3101 CL 5.2.4, direct tensile strength concrete
    #f_conc_tens = 1.4 * 1E3  # [kPa] Defined as this for this section in Henderson Thesis
    #eps_conc_tens = f_conc_tens / e_mod_conc  # Park and Paulay 1975 assume linear
    #eps_unconf_conc_max = 0.004  # Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)

    #eps_conc_crush = 0.003  # Concrte strain at crushing - Koopaee (2015) from NZS 3101?? taken from Fig 2.1 fc=20MPa curve
    
    ft = 0.8*np.sqrt(fpc)

    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush_unconf,
                                                  ec=Ec, fct=ft, et=0)  # unconfined concrete properties
    # conc_unconf = o3.uniaxial_material.Concrete01(osi, fpc=-f_c28 * 0.8, epsc0=-eps_unconf_conc_max, fpcu=0.0,
    #                                               eps_u=-eps_conc_crush)
    peak_disps = np.array([-0.0001, -0.00005, -0.001])
    peak_disps = np.array([-0.01])
    rec_d = 0.00001
    disps = []
    init = 0
    for dd in range(len(peak_disps)):
        sgn = np.sign(peak_disps[dd] - init)
        disps += list(np.arange(init, peak_disps[dd], sgn * rec_d))
        init = peak_disps[dd]
    disp, react = o3.tools.run_uniaxial_disp_driver(osi, conc_unconf, disps)
    ax.plot(disp, react, label=label, c=cbox(0), ls="--")



def add_conf_concrete_to_plot(ax, fpc, eps_conf_conc_max, eps_conc_crush_conf, Ec, label='O3SeesPy Confined Concrete (Concrete 04)'):
    osi = o3.OpenSeesInstance(ndm=1, ndf=1, state=3)
    
    from bwplot import cbox

    #f_c28 = 17 * 1E3  # [kPa] 28 day compressive strength
    #e_mod_conc = (3320 * np.sqrt(f_c28) + 6900) * 1E3  # [kPa] NZS3101 - 2006
    # f_conc_tens = (0.38 * np.sqrt(f_c28))  # [kPa] NZS3101 CL 5.2.4, direct tensile strength concrete
    #f_conc_tens = 1.4 * 1E3  # [kPa] Defined as this for this section in Henderson Thesis
    #eps_conc_tens = f_conc_tens / e_mod_conc  # Park and Paulay 1975 assume linear
    #eps_conf_conc_max = 0.005  # Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)
    #eps_conc_crush_conf = 0.02
    
    ft = 0.8*np.sqrt(fpc)

    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf,
                                                  ec=Ec, fct=ft, et=0)  # Confined concrete paramters
    peak_disps = np.array([-0.0001, -0.00005, -0.001])
    peak_disps = np.array([-0.01])
    rec_d = 0.00001
    disps = []
    init = 0
    for dd in range(len(peak_disps)):
        sgn = np.sign(peak_disps[dd] - init)
        disps += list(np.arange(init, peak_disps[dd], sgn * rec_d))
        init = peak_disps[dd]
    disp, react = o3.tools.run_uniaxial_disp_driver(osi, conc_conf, disps)
    # plt.legend(loc='upper left')
    ax.plot(disp, react, label=label, c=cbox(0))


def add_rebar_to_plot(ax, fy, Es, label="O3SeesPy Steel (Steel 01)"):
    osi = o3.OpenSeesInstance(ndm=1, ndf=1, state=3)
    from bwplot import cbox
    #fy = 300 * 1E3  # [kPa] steel yield strength
    #Es = 200 * 1E6  # [kPa] youngs modulus - steel
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.0015)  # Reinforcing steel
    peak_disps = np.array([-0.0001, -0.00005, -0.001])
    peak_disps = np.array([-0.1])
    rec_d = 0.00001
    disps = []
    init = 0
    for dd in range(len(peak_disps)):
        sgn = np.sign(peak_disps[dd] - init)
        disps += list(np.arange(init, peak_disps[dd], sgn * rec_d))
        init = peak_disps[dd]
    disp, react = o3.tools.run_uniaxial_disp_driver(osi, rebar, disps)
    ax.plot(disp, react, label=label, c=cbox(0))


def add_mander_model_conf_concrete_to_plot(ax, fpc, dels, eco, esm, s, Dh, clb, Ast, fy, Ec, wi, b, d, ncx, ncy, label="CUMBIA Confined Concrete (Mander et al. (1988))" ):
    # Mander Confined concrete model - Rectangular cross section
    from bwplot import cbox
    sp = s - Dh
    Ash = 0.25*np.pi*(Dh**2)
    
    bc = b - 2*clb + Dh
    dc = d - 2*clb + Dh
    Asx = ncx*Ash
    Asy = ncy*Ash
    Ac = bc*dc
    rocc = Ast/Ac
    rox = Asx/(s*dc)
    roy = Asy/(s*bc)
    ros = rox + roy
    n = 0
    for i in range(0, len(wi)):
        n += wi[i]**2 #takes prior value and adds next value
    ke = (((1 - n)/(6*bc*dc))*(1 - sp/(2*bc))*(1 - sp/(2*dc)))/ (1 - rocc)
    ro = 0.5*ros
    fpl = ke*ro*fy
    
    fpcc = (-1.254 + 2.254*np.sqrt(1 + 7.94*fpl/fpc) - 2*fpl/fpc)*fpc
    ecc = eco*(1 + 5*(fpcc/fpc-1))
    Esec = fpcc/ecc
    r = Ec/(Ec-Esec)
    ecu = 1.5*(0.004 + 1.4*ros*fy*esm/fpcc)
    
    ec = np.arange(0, ecu, dels)
    fc = []
    for j in range(0, len(ec)):
        x = (1/ecc)*ec[j]
        fc.append(fpcc*x*r/(r-1+x**r))
    
    ax.plot(ec, fc, label=label, c=cbox(1))
    

def add_mander_model_unconf_concrete_to_plot(ax, fpc, Ec, eco, espall, dels, label="CUMBIA Unconfined Concrete (Mander et al. (1988))"): #reads values from create()
    # Mander unconfined concrete model
    from bwplot import cbox
    fcu = []
    
    ec = np.arange(0, espall, dels)
    Escu = fpc / eco
    ru = Ec / (Ec-Escu)
    fcu_list = []

    
    for i in range(0, len(ec)):
        xu = ec[i] / (eco)
        if ec[i] < 2*eco:
            fcu = fpc*xu*ru/(ru-1 + xu**(ru))
        if ec[i] >= 2*eco and ec[i] <= espall:
                fcu = fpc*(2*ru/(ru-1+2**ru))*(1-(ec[i] - 2*eco)/(espall - 2*eco))
        if ec[i] > espall:
            fcu = 0
        fcu_list.append(fcu)

            
    ax.plot(ec, fcu_list, label=label, c=cbox(1), ls="--")
    
def raynor(ax, Es, fy, fsu, esh, esu, dels, C1, Ey, label="Raynor"):
    """

    :param Es:
    :param fy:
    :param fsu:
    :param esh:
    :param esu:
    :param dels: delta strain for default material models
    :param C1: defines strain hardening curve in the Raynor model [2-6]
    :param Ey: # slope of the yield plateau (MPa)
    :return:
    """
    es = np.linspace(0, esu, int(esu / dels + 1))
    fs = es * 0
    ey = fy / Es;
    fsh = fy + (esh - ey) * Ey
    
    es_list = []
    fs_list = []

    for i in range(len(es)):
        if es[i] < ey:
            fs[i] = Es * es[i]

        if es[i] >= ey and es[i] <= esh:
            fs[i] = fy + (es[i] - ey) * Ey

        if es[i] > esh:
            fs[i] = fsu - (fsu - fsh) * (((esu - es[i]) / (esu - esh)) ** C1)

    return es, fs
    es_list.append(es)
    fs_list.append(fs)
    
    


    print("es/fs:", es, fs)
    ax.plot(es_list, fs_list, label=label)


def king(ax, Es, fy, fsu, esh, esu, dels, label="CUMBIA Steel (King (1986))"):
    from bwplot import cbox
    r = esu - esh
    m = ((fsu / fy) * ((30 * r + 1) ** 2) - 60 * r - 1) / (15 * (r ** 2))
    es = np.linspace(0, esu, int(esu / dels + 1))
    fs = es * 0
    ey = fy / Es

    for i in range(len(es)):
        if es[i] < ey:
            fs[i] = Es * es[i]
        if es[i] >= ey and es[i] <= esh:
            fs[i] = fy
        if es[i] > esh:
            fs[i] = ((m * (es[i] - esh) + 2) / (60 * (es[i] - esh) + 2) + (es[i] - esh) * (60 - m) / (
                        2 * ((30 * r + 1) ** 2))) * fy

    #return es, fs
    ax.plot(es, fs, label=label, c=cbox(1))



def create():
    
    # from bwplot import cbox
    
    #Mander inputs - also used in opensees inputs
    fpc = 17.2 #Conc compressive strength (MPa)
    eps_unconf_conc_max = 0.002 # Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)
    eps_conf_conc_max = 0.002
    eps_conc_crush_unconf = 0.004 # Concrte strain at crushing unconf.
    eps_conc_crush_conf = 0.03 #conf conc strain at crushing conf.
    Ec = 5000*np.sqrt(fpc) #Conc modulus of elasticity (MPa) in Mander script for auto calc
    Ast = 452.3893421 #Total area of longtitudinal steel
    Dh = 6 #diameter of transverse reinforcement (mm)
    clb = 44 #cover to longtitudinal bars (mm)
    s = 600 #spacing of transverse steel (mm)
    fy = 300 #MPa
    Es = 200000 #Steel modulus
    fyh = 300
    eco = 0.002 #Unconfined concrete strain (ususally 0.002 for normal weight or 0.004 for lightweight)
    esm = 0.15 #max transverse steel strain (usually ~0.1 - 0.15)
    espall = 0.0064 #Maximum uncon. conc. strain (usually 0.0064)
    section = 1
    D = 1
    d = 600 #section depth (mm)
    b = 150 #section width (mm)
    ncx = 2 #number of legs, transverse steel x_dir (confinement)
    ncy = 2 #number of transverse steel legs, y-dir (shear)
    wi = [0, 0, 0, 0] #Vector with clear distances between - use zero for automatic calc using mander model
    dels = 0.0001 #Delta strain for default material models (0.0001)
    fsu = 450 #MPa
    esh = 0.008  # long steel strain for strain hardening (usually 0.008)*
    esu = 0.10  # long. steel maximum strain (usually ~0.10-0.15)*
    Ey = 350  # slope of the yield plateau (MPa)
    C1 = 3.5  # defines strain hardening curve in the Raynor model [2-6]
    
    
    
    
    bf, ax = plt.subplots(nrows=2, sharex='col')
    #Add concrete models to plot
    add_unconf_concrete_to_plot(ax[0], fpc, espall, eps_unconf_conc_max, eps_conc_crush_unconf, Ec)
    add_conf_concrete_to_plot(ax[0], fpc, eps_conf_conc_max, eps_conc_crush_conf, Ec)
    add_mander_model_conf_concrete_to_plot(ax[0], fpc, dels, eco, esm, s, Dh, clb, Ast, fy, Ec, wi, b, d, ncx, ncy)
    add_mander_model_unconf_concrete_to_plot(ax[0], fpc, Ec, eco, espall, dels) #This passes the variables to the relevant function
    
    #Add steel models to plot
    add_rebar_to_plot(ax[1], fy, Es)
    king(ax[1], Es, fy, fsu, esh, esu, dels)
    raynor(ax[1], Es, fy, fsu, esh, esu, dels, C1, Ey)
    
    # ax[0].plot(disp, force, c=cbox(0), label="CUMBIA")
    # ax[0].plot(disp1, force1, c=cbox(1), label=‘O3)
    
    
    ax[0].set_ylabel('Stress [MPa]')
    ax[1].set_ylabel('Stress [MPa]')
    ax[-1].set_xlabel('Strain')
    ax[0].set_ylim(0, 20)
    ax[1].set_ylim(0, 450)
    ax[-1].set_xlim(-0.001, 0.02)
    # plt.legend(bbox_to_anchor=(1.02, 0.1), loc='upper left', borderaxespad=0)
    ax[0].legend(bbox_to_anchor=(0, 1.02, 0.5, 0.5), loc='best', borderaxespad=0)
    ax[1].legend(loc='best', borderaxespad=0)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    create() #running the create function first, then runs the other ones



