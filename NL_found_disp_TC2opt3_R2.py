# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 12:26:03 2022

@author: mda153
"""

"""
Created on Mon Mar  7 12:51:03 2022

@author: mda153
"""

import o3seespy as o3
import numpy as np
import matplotlib.pyplot as plt
import sfsimodels as sm


def get_rc_fibre_section(osi, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z):
    """

    y - vertical
    z - horizontal

    :param osi:
    :param ele_nodes:
    :param sm_sect:
    :param conc_conf_mat:
    :param conc_unconf_mat:
    :param rebar_mat:
    :param nf_core_y:
    :param nf_core_z:
    :param nf_cover_y:
    :param nf_cover_z:
    :return:
    """
    import numpy as np
    fc = 30.0e6 #Pa
    fy = 300.0e6  #Pa
    Es = 200.0e9  # Steel modulus [Pa]
    


    col_sect = sm.sections.RCDetailedSection(depth=0.6, width=0.3)
    col_sect.layer_depths = [0.04, 0.56]  # 40mm cover
    col_sect.bar_diams = [[0.016, 0.016, 0.016], [0.016, 0.016, 0.016]]  # 16mm bars
    col_sect.bar_centres = [[0.04, 0.56], [0.04, 0.56]]
    rc_mat = sm.materials.ReinforcedConcreteMaterial(fc=fc, fy=fy,
                                                     e_mod_steel=Es, poissons_ratio=0.18)
    
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.01)  # Reinforcing steel
    area_steel = 0.025 ** 2 / 4 * np.pi * 8
    a_sect = col_sect.depth * col_sect.width
    rho_steel = area_steel / a_sect
    col_sect.rc_mat = rc_mat

    sm_sect = col_sect

    edge_y = sm_sect.depth / 2.0
    edge_z = sm_sect.width / 2.0
    core_y = edge_y - sm_sect.cover_h
    core_z = edge_z - sm_sect.cover_w

    if hasattr(sm_sect, 'gj'):
        gj = sm_sect.gj
    else:
        gj = None
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
    _pi = 3.14159265
    indiv = 0
    if indiv:
        for k in range(len(sm_sect.layer_depths)):
            n_bars = len(sm_sect.bar_centres[k])
            y_pos = sm_sect.layer_depths[k] - sm_sect.depth / 2
            for i in range(n_bars):
                area = sm_sect.bar_diams[k][i] ** 2 / 4 * _pi
                o3.section.gen_fibre_section(osi, y_pos, z_loc=sm_sect.bar_centres[k][i], area=area, mat=rebar)
    else:
        import numpy as np
        for k in range(len(sm_sect.layer_depths)):
            n_bars = len(sm_sect.bar_centres[k])
            bar_area = np.mean(sm_sect.bar_diams[k]) ** 2 / 4 * _pi
            y_pos = sm_sect.layer_depths[k] - sm_sect.depth / 2
            o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[y_pos, core_z], end=[y_pos, -core_z])
    return sect


def run(lf, dx, xd, yd, ksoil, udl, axial_load, max_curve, num_incr, stype="rc"):
    print('running: ', stype)
    d = 0.600  # section depth (m)
    b = 0.300  # section width (m)

    # Mander inputs - also used in opensees inputs - From View Material Behavior Script
    fpc = 30.0e6  # Conc compressive strength (Pa)
    # eps_unconf_conc_max = 0.002  # Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)
    eps_conf_conc_max = 0.002
    eps_conc_crush_unconf = 0.004  # Concrte strain at crushing unconf.
    eps_conc_crush_conf = 0.03  # conf conc strain at crushing conf.
    Ec = 5000 * np.sqrt(fpc)  # Conc modulus of elasticity (Pa) in Mander script for auto calc
    # Ast = 804.248  # Total area of longtitudinal steel
    # Dh = 12  # diameter of transverse reinforcement (mm)
    # clb = 44  # cover to longtitudinal bars (mm)
    # s = 450  # spacing of transverse steel (mm)
    fy = 300.0e6  #Pa
    Es = 200.0e9  # Steel modulus [Pa]
    # fyh = 1
    # eco = 0.002  # Unconfined concrete strain (ususally 0.002 for normal weight or 0.004 for lightweight)
    # esm = 0.15  # max transverse steel strain (usually ~0.1 - 0.15)
    espall = 0.0064  # Maximum uncon. conc. strain (usually 0.0064)
    # section = 1
    # D = 1

    nf_core_y = 16  # number of fibers in Y-Y dir - for the concrete core - (10 - 20 is a good number here)
    nf_core_z = 16  # number of fibers in Z-Z dir
    nf_cover_y = 20  # number of fibers Y-Y dir, including cover concrete
    nf_cover_z = 20  # number of fibers Z-Z dir, including cover concrete
    # n_bars = 2
    # bar_area = 201.0619298  # [mm^2] area of one 16mm diameter reinforcing bar

    osi = o3.OpenSeesInstance(ndm=2, ndf=3)  # 2D with 3 DOF
    s_width = b # section width [m]
    s_depth = d #section depth [m]

    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf,
                                                ec=Ec, fct=0, et=0)  # Confined concrete paramters
    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-espall, epscu=-eps_conc_crush_unconf,
                                                  ec=Ec, fct=0, et=0)  # unconfined concrete properties
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.01)  # Reinforcing steel


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
    # spring_mats = []
    "Soil Nodes"
    for i in range(len(nx)):
        fd_nds.append(o3.node.Node(osi, nx[i], 0.0))
        sl_nds.append(o3.node.Node(osi, nx[i], 0.0))
        o3.Mass(osi, fd_nds[-1], 1.0, 1.0, 1.0)

        dettach = 1
        mat_base = o3.uniaxial_material.ElasticPP(osi, 1*ks[i], 1*py[i] / ks[i], -1 * py[i] / ks[i]) #low tension stiffness
        if dettach:
            mat_obj2 = o3.uniaxial_material.Elastic(osi, 1000 * ks[i], eneg=0.001 * ks[i])
            # mat_obj2 = o3.uniaxial_material.Elastic(osi, 0.001 * ks[i], eneg=1000 * ks[i])
            mat = o3.uniaxial_material.Series(osi, [mat_base, mat_obj2])
        else:
            mat = o3.uniaxial_material.ElasticPP(osi, ks[i], 1*py[i] / ks[i], -1 * py[i] / ks[i]) #remove tension capacity
        #mat = o3.uniaxial_material.Elastic(osi, ks[i])
        sl_eles.append(o3.element.ZeroLength(osi, [fd_nds[-1], sl_nds[-1]], mats=[mat], dirs=[o3.cc.Y]))
        "Foundation nodes"
        if i != 0:
            e_mod = Ec # Pa
            area = s_width * s_depth
            i_z = s_depth ** 3 * s_width / 12

            if stype == "elastic":
                bot_sect = o3.section.Elastic2D(osi, e_mod, area, i_z)
                integ = o3.beam_integration.Lobatto(osi, bot_sect, big_n=5)
                # fd_eles.append(o3.element.ElasticBeamColumn2D(osi, [fd_nds[i], fd_nds[i-1]], area=area, e_mod=e_mod, iz=i_z, transf=transf))

            elif stype == "rc":
                bot_sect = get_rc_fibre_section(osi, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z)
                integ = o3.beam_integration.Lobatto(osi, bot_sect, big_n=5)
            else:
                raise ValueError('set stype to elastic or rc')
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
    assert np.isclose(min(ndisps), -udl / ksoil, rtol=0.01), (min(ndisps), max(ndisps), -udl / ksoil)

    ts0 = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts0)
    yd = np.array(yd) -udl / ksoil
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
    max_steps = 1000
    fd_incs = ymin / max_steps
    o3.integrator.DisplacementControl(osi, fd_nds[max_ind], dof=o3.cc.Y, incr=fd_incs, num_iter=100) 
    
    ndr = o3.recorder.NodesToArrayCache(osi, fd_nds, dofs=[o3.cc.Y], res_type='disp') #Node displacement recorder
    efr = o3.recorder.ElementsToArrayCache(osi, fd_eles, arg_vals=['section', 1, 'force'])  # Not working for NL but there are some other inputs needed
    print("hi found elsse", len(fd_eles))
    # print(efr.parameters)
    ndisps = [[]]
    mom = [[]]
    print("nndoes")
    print(nnodes)
    for j in range(nnodes):
        ndisps[0].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
    for j in range(nnodes - 1):
        mom[0].append(o3.get_ele_response(osi, fd_eles[j], 'force')[2])
    print('hi', ndisps)
    for i in range(max_steps + 20):
        fail = o3.analyze(osi, 1)
        if fail:
            raise ValueError()
        ndisps.append([])
        mom.append([])
        for j in range(nnodes):
            ndisps[-1].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
            
        for j in range(nnodes-1):
            # print('force: ', o3.get_ele_response(osi, fd_eles[j], 'force'))
            mom[-1].append(o3.get_ele_response(osi, fd_eles[j], 'force')[2])
            # mom.append(o3.get_ele_response(osi, fd_eles[4], 'force')[2])
            
        # print('y_disps: ', [f'{dy:.3}' for dy in ndisps[-1]])
        y_disp_max = o3.get_node_disp(osi, sl_nds[max_ind], dof=o3.cc.Y)
        if -y_disp_max > -ymin:
            
            print('break: ', y_disp_max, ymin)
            break

    for j in range(nnodes - 1):
        print('force: ', o3.get_ele_response(osi, fd_eles[j], 'force'))
    o3.wipe(osi)  # if you are using a recorder you need to have this line
    node_disps = ndr.collect()
    forces = efr.collect()
    return nx, node_disps, forces, xd0n,xd1n, yd0, yd1, mom


def create(): #creates plot
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
        Uniform distributed load on foundation [Pa] #want to do this in N/m 
    :param xd:
        (x0, x1) positions of displaced section of soil
    :param yd:
        (y0, y1) displacements of section of soil at x0 and x1 (note should be -ve)
    :return:
    """
    bf, ax = plt.subplots(nrows=2)
    #grabs relevant values from run function
    stypes = ["elastic", "rc"]
    cols = ['r', 'b']
    # stypes = [stypes[0]]
    for ss, stype in enumerate(stypes):
        nx_elastic, ndisps_elastic, forces_el, xd0n_elastic, xd1n_elastic, yd0_elastic, yd1_elastic,mom_el = run(lf=10, dx=0.5, xd=(3,7), yd=(-0.1,-0.1), ksoil=25000e3, udl=8.59e3, axial_load=0, max_curve=0.003, num_incr=500, stype=stype)
        # nx_rc, ndisps_rc, xd0n_rc, forces_rc, xd1n_rc, yd0_rc, yd1_rc, mom_rc = run(lf=10, dx=0.5, xd=(3,7), yd=(-0.1,-0.1), ksoil=2.5e3, udl=0.03e3, axial_load=0, max_curve=0.003, num_incr=500, stype="rc")

        # print(forces_el)
        # print(forces_rc)

        # print(mom_el)
        print('hi - len mom el',len(mom_el))
        print('mom_rc')
        # print(mom_rc)
        # print('hi len mom rc ', len(mom_rc))
        print(len(nx_elastic))

        # length_found = 10 #foundation lenght in m

        "elastic model plotting"
        ax[0].plot(nx_elastic, ndisps_elastic[0], label=f'Initial {stype}', c=cols[ss], ls='-.')
        # ax[0].plot(nx_elastic, ndisps_elastic[int(len(ndisps_elastic) / 2)], label='Linear @ 50%', c="b", ls = "--") #Plotting for 50% but this is no longer 50%!! Depends on array length!!
        ax[0].plot(nx_elastic, ndisps_elastic[-1], label=f'Final {stype}', c=cols[ss])

        # ax[0].plot([xd0n_elastic, xd1n_elastic], [yd0_elastic, yd1_elastic], c='r', label='Imposed')

        # print(forces_el)

        # mom = forces_el[:, 2::6]
        # print(mom)

        el_x = (nx_elastic[1:] + nx_elastic[:-1])/2
        print(el_x)
        print('hi len exl', len(el_x), el_x)
        print("he len mom_el", len(mom_el), mom_el)
        ax[1].plot(el_x, mom_el[-1], label=f'{stype}', c=cols[ss])
        # ax[1].plot(el_x, mom_el[20:40], label='Linear', c="r") #plots

        "dispbeam column plotting"

        # ax[0].plot(nx_rc, ndisps_rc[0], label='Initial Non-Linear', c="k", ls="-.") #non linear
        # ax[0].plot(nx_rc, ndisps_rc[int(len(ndisps_elastic) / 2)], label='Non-Linear', c="k", ls="--")
        # ax[0].plot(nx_rc, ndisps_rc[-1], label='Non-Linear Final', c="k")


        # mom = forces_rc[:, 2::6]
        # print(len(mom))

        # print(mom)
        el_x = (nx_elastic[1:] + nx_elastic[:-1]) / 2
    # print(el_x)
    # print("hi len mom rc", len(mom_rc), mom_rc)
    #
    # ax[1].plot(el_x, mom_rc[10:30], label='Non-Linear',c='k')
    # ax[1].plot(el_x, mom_rc[10:30], label='Non-Linear', c="k")
    # ax[0].plot([xd0n_rc, xd1n_rc], [yd0_rc, yd1_rc], c='r', label='Imposed')
    ax[0].set_xlabel('Foundation Length (m)')
    ax[0].set_ylabel('Settlement (m)')
    ax[1].set_ylabel("Bending Moment (Nm)")
    ax[1].set_xlabel("Foundation Length (m)")
    # plt.grid()
    
    ax[1].legend(bbox_to_anchor =(1,1))
    ax[0].legend(bbox_to_anchor =(1,1))

    plt.show()





if __name__ == '__main__':
    
    # "get_moment_curvature function"
    # mom, curve, d, b, Ec, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z = get_moment_curvature(axial_load=0,max_curve=0.002,num_incr=500)
    # plt.plot(curve*1000, mom/1000/1000)
    # plt.grid()
    # axes = plt.axes()
    # axes.set_xlabel("Curvature [1/m]")
    # axes.set_ylabel("Moment [kNm]")
    
    # plt.show()
    
    
    "Impose foundation displacement function"
    create()
    # run(lf=10, dx=0.5)
