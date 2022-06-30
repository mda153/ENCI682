# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 15:03:54 2022

@author: mda153
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 12:51:03 2022

@author: mda153
"""

import o3seespy as o3
import numpy as np
import matplotlib.pyplot as plt
import o3seespy.extensions
from loc_o3fibre import get_rc_fibre_section

def get_moment_curvature(axial_load, max_curve, num_incr):
    osi = o3.OpenSeesInstance(ndm=2, ndf=3, state=3)
    
    #Mander inputs - also used in opensees inputs - From View Material Behavior Script
    fpc = 30 #Conc compressive strength (MPa)
    eps_unconf_conc_max = 0.002 # Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)
    eps_conf_conc_max = 0.002
    eps_conc_crush_unconf = 0.004 # Concrte strain at crushing unconf.
    eps_conc_crush_conf = 0.03 #conf conc strain at crushing conf.
    Ec = 5000*np.sqrt(fpc) #Conc modulus of elasticity (MPa) in Mander script for auto calc
    Ast = 804.248 #Total area of longtitudinal steel
    Dh = 12 #diameter of transverse reinforcement (mm)
    clb = 44 #cover to longtitudinal bars (mm)
    s = 450 #spacing of transverse steel (mm)
    fy = 300 #MPa
    Es = 205000 #Steel modulus [MPa]
    fyh = 1
    eco = 0.002 #Unconfined concrete strain (ususally 0.002 for normal weight or 0.004 for lightweight)
    esm = 0.15 #max transverse steel strain (usually ~0.1 - 0.15)
    espall = 0.0064 #Maximum uncon. conc. strain (usually 0.0064)
    section = 1
    D = 1
    d = 600 #section depth (mm)
    b = 300 #section width (mm)
    ncx = 2 #number of legs, transverse steel x_dir (confinement)
    ncy = 2 #number of transverse steel legs, y-dir (shear)
    wi = [0, 0, 0, 0] #Vector with clear distances between - use zero for automatic calc using mander model
    dels = 0.0001 #Delta strain for default material models (0.0001)
    
    
    
    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf,
                                                  ec=Ec, fct=0, et=0)  # Confined concrete paramters
    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-espall, epscu=-eps_conc_crush_unconf,
                                                  ec=Ec, fct=0, et=0)  # unconfined concrete properties
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.01)  # Reinforcing steel

  
    
    #Defining concrete beam parameters

    h = d #Section height [mm]
    b = b #Section width [mm]
    cover = 40 #[mm] Concrete cover

    gj = 1.0E10 # linear elastic torsional stiffness
   
    nf_core_y = 16 #number of fibers in Y-Y dir - for the concrete core - (10 - 20 is a good number here)
    nf_core_z = 16 #number of fibers in Z-Z dir
    nf_cover_y = 20 #number of fibers Y-Y dir, including cover concrete
    nf_cover_z = 20 #number of fibers Z-Z dir, including cover concrete
    n_bars = 2
    bar_area = 201.0619298 #[mm^2] area of one 16mm diameter reinforcing bar
    
#Defining beam section
    edge_y = h / 2.0
    edge_z = b / 2.0
    core_y = edge_y - cover
    core_z = edge_z - cover



    sect = o3.section.Fiber(osi, gj=gj)
    # define the core patch
    o3.patch.Quad(osi, conc_conf, nf_core_z, nf_core_y,  # core, counter-clockwise
                  crds_i=[-core_y, core_z],
                  crds_j=[-core_y, -core_z],
                  crds_k=[core_y, -core_z],
                  crds_l=[core_y, core_z])

    o3.patch.Quad(osi, conc_unconf, 1, nf_cover_y,  # right cover, counter-clockwise (diagonals at corners of cover)
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

    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[-core_y, -core_z], end=[-core_y, core_z]) #defining bottom layer of reinforcement
    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[core_y, -core_z], end=[core_y, core_z]) #top layer of reinforcement

    spacing_y = 2 * core_y / (n_bars - 1)
    remaining_bars = n_bars - 1
    o3.layer.Straight(osi, rebar, remaining_bars, bar_area,
                      start=[-core_y, -core_z],
                      end=[-core_y, core_z]) #bottom layer of reinforcement
    o3.layer.Straight(osi, rebar, remaining_bars, bar_area,
                      start=[core_y, -core_z],
                      end=[core_y, core_z]) #top layer of reinforcement
    
    
#layer command generates a row of fibers along a geometric-arc, in a straight line
  #layer command generates a row of fibers along a geometric-arc, in a straight line
   

    n1 = o3.node.Node(osi, 0.0, 0.0)
    n2 = o3.node.Node(osi, 0.0, 0.0)
    o3.Fix3DOF(osi, n1, 1, 1, 1)
    o3.Fix3DOF(osi, n2, 0, 1, 0)
    ele = o3.element.ZeroLengthSection(osi, [n1, n2], sect)

    nd = o3.recorder.NodeToArrayCache(osi, n2, dofs=[3], res_type='disp')
    nm = o3.recorder.NodeToArrayCache(osi, n1, dofs=[3], res_type='reaction')
    


    ts = o3.time_series.Constant(osi)
    o3.pattern.Plain(osi, ts)
    o3.Load(osi, n2, load_values=[axial_load, 0.0, 0.0])

    o3.system.BandGeneral(osi)
    o3.numberer.Plain(osi)
    o3.constraints.Plain(osi)
    o3.test.NormUnbalance(osi, tol=1.0e-9, max_iter=10)
    o3.algorithm.Newton(osi)
    o3.integrator.LoadControl(osi, incr=0.0)
    o3.analysis.Static(osi)
    o3.analyze(osi, 1)

    #
    ts = o3.time_series.Linear(osi)
    o3.pattern.Plain(osi, ts)
    o3.Load(osi, n2, load_values=[0.0, 0.0, 1.0])

    d_cur = max_curve / num_incr

    o3.integrator.DisplacementControl(osi, n2, o3.cc.DOF2D_ROTZ, d_cur, 1, d_cur, d_cur)
    o3.analyze(osi, num_incr)
    o3.wipe(osi)
    curvature = nd.collect()
    moment = -nm.collect()
    print(moment)
    print(curvature)
    return moment, curvature, d, b, Ec, sect ,conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z



if __name__ == '__main__':

    mom, curve, d, b, Ec,sect,conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z = get_moment_curvature(axial_load=0,max_curve=0.002,num_incr=500)
    plt.plot(curve*1000, mom/1000/1000)
    plt.grid()
    axes = plt.axes()
    axes.set_xlabel("Curvature [1/m]")
    axes.set_ylabel("Moment [kNm]")
    
    

    plt.show()


def run(lf=10, dx=0.5, xd=(3,7), yd=(-0.2,-0.2), ksoil=2.5e3, udl=0.03e3, axial_load=0, max_curve=0.003, num_incr=500):
    """
    Run an analysis imposing a uniform load, on a foundation with soil springs

    Then impose a displacement on the ends of some of the soil springs

    :param lf:
        Foundation length [m]
    :param dx:
        Target spacing of springs [m]
    :param s_depth:
        Foundation section depth [m]
    :param ksoil:
        Subgrade stiffness of soil [N/m3]
    :param udl:
        Uniform distributed load on foundation [Pa]
    :param xd:
        (x0, x1) positions of displaced section of soil
    :param yd:
        (y0, y1) displacements of section of soil at x0 and x1 (note should be -ve)
    :return:
    """
    moment, curvature, d, b, Ec,sect ,conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z = get_moment_curvature(axial_load, max_curve, num_incr)
    print(moment)
    print(curvature)
    print(d)
    print(b)
    print(Ec)
    osi = o3.OpenSeesInstance(ndm=2, ndf=3)  # 2D with 3 DOF
    s_width = b # section width [m]
    s_depth = d #section depth [m]
    
    nnodes = max(int(lf / dx) + 1, 2)  # number of nodes
    nx = np.linspace(0, lf, nnodes)  # x-position of nodes
    x_trib = np.ones_like(nx)  # Trib area
    x_trib[1:] += np.diff(nx)
    x_trib[:-1] += np.diff(nx)

    ks = ksoil * x_trib * s_width  # N/m (TODO: define this better, see Henderson thesis for details)
    py = 1.0e3 * x_trib * s_width  # N Capacity (TODO: define this better, see Henderson thesis for details)

    transf = o3.geom_transf.Linear2D(osi, [])  # Using simple transform - no P-delta effects

    fd_nds = []  # foundation nodes
    sl_nds = []  # soil nodes
    sl_eles = []
    fd_eles = []
    "Soil Nodes"
    for i in range(len(nx)):
        fd_nds.append(o3.node.Node(osi, nx[i], 0.0))
        sl_nds.append(o3.node.Node(osi, nx[i], 0.0))
        o3.Mass(osi, fd_nds[-1], 1.0, 1.0, 1.0)
        mat = o3.uniaxial_material.ElasticPP(osi, ks[i], 1*py[i] / ks[i], 0.001 * py[i] / ks[i]) #remove tension capacity
        #mat = o3.uniaxial_material.Elastic(osi, ks[i])
        sl_eles.append(o3.element.ZeroLength(osi, [fd_nds[-1], sl_nds[-1]], mats=[mat], dirs=[o3.cc.Y]))
        "Foundation nodes"
        if i != 0:
            e_mod = Ec  # Pa
            area = s_width * s_depth
            i_z = s_depth ** 3 * s_width / 12
            bot_sect = o3.section.Elastic2D(osi, e_mod, area, i_z)
            # bot_sect = get_moment_curvature(osi, sect ,conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z)
            integ = o3.beam_integration.Lobatto(osi, bot_sect, big_n=5) 
            # fd_eles.append(o3.element.ElasticBeamColumn2D(osi, [fd_nds[i], fd_nds[i-1]], area=area, e_mod=e_mod,
            #                                               iz=i_z, transf=transf))
            fd_eles.append(o3.element.DispBeamColumn(osi, [fd_nds[i], fd_nds[i-1]], transf=transf, integration=integ))
        



    o3.Fix3DOFMulti(osi, sl_nds, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)  # Fix all the soil nodes
    o3.Fix3DOF(osi, fd_nds[0], o3.cc.FIXED, o3.cc.FREE, o3.cc.FREE)

    # Define load
    ploads = udl * s_width * x_trib
    ts0 = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts0)
    for i in range(nnodes):
        o3.Load(osi, fd_nds[i], [0, -ploads[i], 0])

    # Run initial static analysis
    o3.constraints.Transformation(osi)
    o3.test_check.NormDispIncr(osi, tol=1.0e-6, max_iter=35, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.numberer.RCM(osi)
    o3.system.SparseGeneral(osi)
    n_steps_gravity = 10
    d_gravity = 1. / n_steps_gravity
    o3.integrator.LoadControl(osi, d_gravity, num_iter=n_steps_gravity)
    o3.analysis.Static(osi)
    o3.analyze(osi, num_inc=n_steps_gravity)
    o3.load_constant(osi, time=0.0)
    ndisps = []
    for i in range(nnodes):
        ndisps.append(o3.get_node_disp(osi, fd_nds[i], dof=o3.cc.Y))
    print('y_disps: ', [f'{dy:.3}' for dy in ndisps])
    assert np.isclose(min(ndisps), max(ndisps)), (min(ndisps), max(ndisps), -udl / ksoil)
    assert np.isclose(min(ndisps), -udl / ksoil), (min(ndisps), max(ndisps), -udl / ksoil)

    ts0 = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts0)
    xd0 = xd[0]  # start position of displaced section
    xd1 = xd[1]  # end position of displaced section
    yd0 = yd[0]  # target displacement at start
    yd1 = yd[1]  # target displacement at end
    ind0 = np.argmin(abs(nx - xd0))
    ind1 = np.argmin(abs(nx - xd1))
    xd0n = nx[ind0]  # actual position of node where displacement will be imposed
    xd1n = nx[ind1]
    for i in range(ind0, ind1+1):
        xi = nx[i]
        ydi = yd0 + (yd1 - yd0) / (xd1n - xd0n) * (xi - xd0n)  # linear interpolate
        o3.SP(osi, sl_nds[i], o3.cc.Y, [ydi])  # impose displacement
    ymin = min([yd0, yd1])
    max_ind = [ind0, ind1][np.argmin([yd0, yd1])]  # use the node that has the largest displacement
    o3.integrator.DisplacementControl(osi, fd_nds[max_ind], dof=o3.cc.Y, incr=ymin / 100, num_iter=10)
    ndisps = [[]]
    for j in range(nnodes):
        ndisps[0].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
    for i in range(100):
        o3.analyze(osi, 1)
        ndisps.append([])
        for j in range(nnodes):
            ndisps[-1].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
        print('y_disps: ', [f'{dy:.3}' for dy in ndisps[-1]])
        
    print(nx)
    plt.plot(nx, ndisps[0], label='Initial Foundation')
    plt.plot(nx, ndisps[49], label='Foundation @ 50%')
    plt.plot(nx, ndisps[99], label='Foundation @ 100%')
    plt.plot([xd0n, xd1n], [yd0, yd1], c='k', label='Imposed')
    plt.xlabel('Foundation Length (m)')
    plt.ylabel('Vertical Displacement (m)')
    plt.grid()

    plt.legend()
    plt.show()
    




if __name__ == '__main__':
    run(lf=10, dx=0.5)