# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 15:27:35 2022

@author: mda153
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import o3seespy as o3
import numpy as np
import o3seespy.extensions


def get_moment_curvature(axial_load=0, max_curve=0.15, num_incr=100):
    osi = o3.OpenSeesInstance(ndm=2, ndf=3, state=3)
    
    fpc = 17.2e6 #Conc compressive strength (Pa)
    eps_unconf_conc_max = 0.002 # Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)
    eps_conf_conc_max = 0.002
    eps_conc_crush_unconf = 0.004 # Concrte strain at crushing unconf.
    eps_conc_crush_conf = 0.03 #conf conc strain at crushing conf.
    Ec = 5000*(fpc**0.5)*1e3 #Conc modulus of elasticity (Pa) in Mander script for auto calc
    print("EC:",Ec)
    d_bar = 0.012 # longtitudinal bar diameter [m]
    total_num_bars = 4
    Ast = np.pi*(d_bar/2)**2*total_num_bars #Total area of longtitudinal steel (m^2)
    Dh = 0.006 #diameter of transverse reinforcement (m)
    clb = 0.044 #cover to longtitudinal bars (m)
    s = 0.3 #spacing of transverse steel (m)
    fy = 300e6 #Pa
    Es = 200e9 #Steel modulus (Pa)
    fyh = 300e6
    eco = 0.002 #Unconfined concrete strain (ususally 0.002 for normal weight or 0.004 for lightweight)
    esm = 0.15 #max transverse steel strain (usually ~0.1 - 0.15)
    espall = 0.0064 #Maximum uncon. conc. strain (usually 0.0064)
    # section = 1
    # D = 1
    d = 0.6 #section depth (mm)
    b = 0.15 #section width (mm)
    ncx = 2 #number of legs, transverse steel x_dir (confinement)
    ncy = 2 #number of transverse steel legs, y-dir (shear)
    wi = [0, 0, 0, 0] #Vector with clear distances between - use zero for automatic calc using mander model
    dels = 0.0001 #Delta strain for default material models (0.0001)
    

    
    ft = 0.8*np.sqrt(fpc)
    # ft = 1.4e3
    
    Z = (b*d**2)/6
    
    mom_crack = Z*ft
    print("mom_crack:", mom_crack)
    
    n_bars_tens = 2
    bars_total = 4
    bar_area = Ast/bars_total
    print("Bar area:", bar_area)
    j = 0.85 #Assumed
    d_eff = d - clb
    
    
    mom_yield = (n_bars_tens*bar_area*fy*j*d_eff)/1000
    print("mom_yeild:", mom_yield)
    
    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf, ec=Ec, fct=ft, et=0) #confined concrete properties
    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush_unconf, ec=Ec, fct=ft, et=0) #unconfined concrte properties
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.0015) #Steel reinforcing properties

    #Defining concrete beam parameters

    h = d #Section height
    b = b #Section width
    cover = clb #Cover depth
    gj = 1.0E10 # Is this G * J? 
    nf_core_y = 40
    nf_core_z = 40
    nf_cover_y = 50
    nf_cover_z = 50
    n_bars = 2 #Number of bars in a layer/row
    bars_total = 4
    bar_area = Ast/bars_total #Area of steel (m^2)
    print(bar_area)


    edge_y = h / 2.0
    edge_z = b / 2.0
    core_y = edge_y - cover
    core_z = edge_z - cover

    sect = o3.section.Fiber(osi, gj=gj)
    # define the core patch
    o3.patch.Quad(osi, conc_conf, nf_core_z, nf_core_y,  # core, counter-clockwise (diagonals at corners)
                  crds_i=[-core_y, core_z],
                  crds_j=[-core_y, -core_z],
                  crds_k=[core_y, -core_z],
                  crds_l=[core_y, core_z])

    o3.patch.Quad(osi, conc_unconf, 1, nf_cover_y,  # right cover, counter-clockwise (diagonals at corners)
                  crds_i=[-edge_y, edge_z],
                  crds_j=[-core_y, core_z],
                  crds_k=[core_y, core_z],
                  crds_l=[edge_y, edge_z])
    o3.patch.Quad(osi, conc_unconf, 1, nf_cover_y,  # left cover
                  crds_i=[-core_y, -core_z],
                  crds_j=[-edge_y, -edge_z],
                  crds_k=[edge_y, -edge_z],
                  crds_l=[core_y, -core_z])
    o3.patch.Quad(osi, conc_unconf, nf_cover_z, 1,  # bottom cover
                  crds_i=[-edge_y, edge_z],
                  crds_j=[-edge_y, -edge_z],
                  crds_k=[-core_y, -core_z],
                  crds_l=[-core_y, core_z])
    o3.patch.Quad(osi, conc_unconf, nf_cover_z, 1,  # top cover
                  crds_i=[core_y, core_z],
                  crds_j=[core_y, -core_z],
                  crds_k=[edge_y, -edge_z],
                  crds_l=[edge_y, edge_z])

    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[-core_y, core_z], end=[-core_y, -core_z])
    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[core_y, core_z], end=[core_y, -core_z]) #Assigns rebar layers

    spacing_y = 2 * core_y / (n_bars - 1)
    remaining_bars = n_bars - 1 #question here
    o3.layer.Straight(osi, rebar, remaining_bars, bar_area,
                      start=[core_y - spacing_y, core_z],
                      end=[-core_y + spacing_y, core_z])
    o3.layer.Straight(osi, rebar, remaining_bars, bar_area,
                      start=[core_y - spacing_y, -core_z],
                      end=[-core_y + spacing_y, -core_z])

    n1 = o3.node.Node(osi, 0.0, 0.0)
    n2 = o3.node.Node(osi, 0.0, 1)
    o3.Fix3DOF(osi, n1, 1, 1, 1)
    o3.Fix3DOF(osi, n2, 0, 0, 0)
    transf = o3.geom_transf.Linear2D(osi, [])
    integ = o3.beam_integration.Lobatto(osi, sect, big_n=5)
    ele = o3.element.DispBeamColumn(osi, [n1, n2], transf=transf, integration=integ)

    # ele = o3.element.ZeroLengthSection(osi, [n1, n2], sect)

    ndr = o3.recorder.NodeToArrayCache(osi, n2, dofs=[3], res_type='disp')
    nmr = o3.recorder.NodeToArrayCache(osi, n1, dofs=[3], res_type='reaction')
    ecr = o3.recorder.ElementToArrayCache(osi, ele, arg_vals=['deformation'])

    ts = o3.time_series.Constant(osi)
    o3.pattern.Plain(osi, ts)
    o3.Load(osi, n2, load_values=[axial_load, 0.0, 0.0])

    o3.system.BandGeneral(osi)
    o3.numberer.Plain(osi)
    o3.constraints.Transformation(osi)
    o3.test.NormUnbalance(osi, tol=1.0e-9, max_iter=10)
    o3.algorithm.Newton(osi)
    o3.integrator.LoadControl(osi, incr=0.0)
    o3.analysis.Static(osi)
    o3.analyze(osi, 1)

    #
    ts = o3.time_series.Linear(osi)
    o3.pattern.Plain(osi, ts)

    d_cur = max_curve / num_incr
    #
    variable_disp_cont = 0
    if variable_disp_cont:
        o3.Load(osi, n2, load_values=[0.0, 0.0, 1.0])
        o3.integrator.DisplacementControl(osi, n2, o3.cc.DOF2D_ROTZ, d_cur, 1, d_cur, d_cur)
    else:
        o3.SP(osi, n2, o3.cc.DOF2D_ROTZ, dof_values=[d_cur])
        o3.integrator.LoadControl(osi, incr=1)
    o3.analyze(osi, num_incr)
    o3.wipe(osi)
    rot = ndr.collect()
    # ele_curvature = ecr.collect()[:, 0]  # Not sure if this is actually curvature.

    moment = -nmr.collect()
    return moment, rot, mom_crack, mom_yield


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    mom, curve ,mom_crack, mom_yield = get_moment_curvature(num_incr=4000, max_curve=0.25)

    # m_cr = mom_crack
    # m_y = mom_yield
    # m_ult = mom[-1]/1000
    
    # print("Moments:", m_cr, m_y, m_ult)
    
    # plt.axhline(m_cr, c="green", label="Cracking Moment (Hand Calc)")
    # plt.axhline(m_y, c="orange", label="Yield Moment (Hand Calc)")
    # plt.axhline(m_ult, c="red", label="Ultimate Moment")
    
    plt.plot(curve, mom/1000)
    print(curve, mom)
    
    curve_list = curve.tolist()
    mom_list = mom.tolist()
    # print(curve_list, mom_list)
    
    np.savetxt('m_curv_opensees_vTEST_1930_59_Hend.csv', (mom_list, curve_list), fmt='%.4e', delimiter=',')# X is an array

    fpc = 17.2e6  # Conc compressive strength (Pa)
    Ec = 5000 * (fpc ** 0.5) * 1E3
    d = 0.6  # section depth (m)
    b = 0.15  # section width (m)
    ei_el = Ec * b * d ** 3 / 12
    moms = np.array([0, 60])
    curvs = moms / ei_el
    plt.plot(curvs, moms, c='k')
    print(curvs, moms)

    plt.grid()
    axes = plt.gca()
    axes.set_xlabel("Curvature [1/m]")
    axes.set_ylabel("Bending Moment [kNm]")
    axes.legend()
    
    

    plt.show()
  