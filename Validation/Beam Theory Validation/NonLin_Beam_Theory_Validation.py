# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 09:34:56 2022

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
    # import numpy as np
    fc = 17.0e6 #[Pa] concrete compressive strength
    fy = 300.0e6  #[Pa] steel yield strength
    Es = 200.0e9  # [Pa] Steel modulus of elasticity 
    depth = 0.6 #[m] section depth
    width = 0.15 #[m] section width
    
    


    col_sect = sm.sections.RCDetailedSection(depth=depth, width=width)
    col_sect.layer_depths = [0.044, 0.556]  # 44mm cover (m units)
    col_sect.bar_diams = [[0.012, 0.012], [0.012, 0.012]]  # 12mm bars - two 12mm bars directly beside each other at top of section - model this as 1 x 24mm bar
    col_sect.bar_centres = [[0.044, 0.106], [0.044, 0.106]]
    rc_mat = sm.materials.ReinforcedConcreteMaterial(fc=fc, fy=fy,
                                                     e_mod_steel=Es, poissons_ratio=0.18)
    
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.01)  # Reinforcing steel

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
        # import numpy as np
        for k in range(len(sm_sect.layer_depths)):
            n_bars = len(sm_sect.bar_centres[k])
            bar_area = np.mean(sm_sect.bar_diams[k]) ** 2 / 4 * _pi
            y_pos = sm_sect.layer_depths[k] - sm_sect.depth / 2
            o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[y_pos, core_z], end=[y_pos, -core_z])
    return sect


def run(lf, dx, xd, yd, ksoil, udl, num_incr, stype="Non-Linear"):
    print('running: ', stype)
    # spring_spacing = 0.5 # [m] spring spacing

    imposed_displacement = 0.1 #[m]
    d = 0.6  # section depth [m]
    b = 0.15  # section width [m]

    # Mander inputs - also used in opensees inputs - From View Material Behavior Script
    fpc = 17.0e6  # Conc compressive strength (Pa)
    ft = 1.40e6 #[Pa] concrete tensile strength
    eps_conf_conc_max = 0.002 
    eps_conc_crush_unconf = 0.004  # Concrte strain at crushing unconf.
    eps_conc_crush_conf = 0.03  # conf conc strain at crushing conf.
    Ec = 5000 * np.sqrt(fpc)*1e3  # Concrete modulus of elasticity (Pa) - eqn in Mander script for auto calc

    fy = 300.0e6  #Pa
    Es = 200.0e9  # Steel modulus [Pa]

    espall = 0.0064  # Maximum uncon. conc. strain (usually 0.0064)
 

    nf_core_y = 8  # number of fibers in Y-Y dir - for the concrete core - (10 - 20 is a good number here)
    nf_core_z = 8  # number of fibers in Z-Z dir
    nf_cover_y = 10  # number of fibers Y-Y dir, including cover concrete
    nf_cover_z = 10  # number of fibers Z-Z dir, including cover concrete
    # n_bars = 2
    # bar_area = 201.0619298  # [mm^2] area of one 16mm diameter reinforcing bar

    osi = o3.OpenSeesInstance(ndm=2, ndf=3)  # 2D with 3 DOF
    s_width = b # section width [m]
    s_depth = d #section depth [m]

    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf,
                                                ec=Ec, fct=ft, et=0)  # Confined concrete paramters
    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-espall, epscu=-eps_conc_crush_unconf,
                                                  ec=Ec, fct=ft, et=0)  # unconfined concrete properties
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.01)  # Reinforcing steel


    nnodes = max(int(lf / dx) + 1, 2)  # number of nodes
    nx = np.linspace(0, lf, nnodes)  # x-position of nodes
    x_trib = np.ones_like(nx)  # Trib area
    x_trib[1:] += np.diff(nx)
    x_trib[:-1] += np.diff(nx)
    print(x_trib)
    ks = ksoil * s_width * x_trib #soil spring stiffness N/m
    
    print('hi ks', ks)
    
    #Could not find this in henderson thesis - does not expicitly mention it
    #Surely the force resistance should be related to the stiffness of each spring then multiply this by the displacement it is subject to

    # ks = ksoil * x_trib *s_width # N/m (TODO: define this better, see Henderson thesis for details) spring stiffness
    #Spring stiffness should be defined as ks = ksoil * section width * spring spacing
    
    py = ks * imposed_displacement*1e0  # N Capacity (TODO: define this better, see Henderson thesis for details) spring force resistance
    # thesis does not define this... I could not find it
    print("hi py", py)
    

    transf = o3.geom_transf.Linear2D(osi, [])  # Using simple transform - no P-delta effects

    fd_nds = []  # foundation nodes
    sl_nds = []  # soil nodes
    sl_eles = [] #soil elements
    fd_eles = [] #foundation beam elements

    "Soil Nodes"
    for i in range(len(nx)):
        fd_nds.append(o3.node.Node(osi, nx[i], 0.0))
        sl_nds.append(o3.node.Node(osi, nx[i], 0.0))
        o3.Mass(osi, fd_nds[-1], 1.0, 1.0, 1.0)

        dettach = 1
        mat_base = o3.uniaxial_material.ElasticPP(osi, 1*ks[i], 1*py[i] / ks[i], -1 * py[i] / ks[i]) #low tension stiffness
        if dettach:
            mat_obj2 = o3.uniaxial_material.Elastic(osi, 10000 * ks[i], eneg=0.00001 * ks[i])
            # mat_obj2 = o3.uniaxial_material.ENT(osi, 1000 * ks[i])
            # mat_obj2 = o3.uniaxial_material.Elastic(osi, 0.001 * ks[i], eneg=1000 * ks[i])
            mat = o3.uniaxial_material.Series(osi, [mat_base, mat_obj2])
        else:
            mat = o3.uniaxial_material.ElasticPP(osi, ks[i], 1*py[i] / ks[i], -1 * py[i] / ks[i]) #remove tension capacity
        #mat = o3.uniaxial_material.Elastic(osi, ks[i])
        sl_eles.append(o3.element.ZeroLength(osi, [fd_nds[-1], sl_nds[-1]], mats=[mat], dirs=[o3.cc.Y]))
        "Foundation nodes"
        if i != 0:
            fc = 17.0e6 #[Pa]
            e_mod = 5000*np.sqrt(fc)*1e3 # Pa
            
            area = s_width * s_depth
            i_z = s_depth ** 3 * s_width / 12

            if stype == "Elastic":
                bot_sect = o3.section.Elastic2D(osi, e_mod, area, i_z)
                integ = o3.beam_integration.Lobatto(osi, bot_sect, big_n=5)
                # fd_eles.append(o3.element.ElasticBeamColumn2D(osi, [fd_nds[i], fd_nds[i-1]], area=area, e_mod=e_mod, iz=i_z, transf=transf))

            elif stype == "Non-Linear":
                bot_sect = get_rc_fibre_section(osi, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z)
                integ = o3.beam_integration.Lobatto(osi, bot_sect, big_n=5)
            else:
                raise ValueError('set stype to Elastic or Non-Linear')
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
    mom = []
    shear = []
    for i in range(nnodes):
        ndisps.append(o3.get_node_disp(osi, fd_nds[i], dof=o3.cc.Y))
    # print('y_disps: ', [f'{dy:.3}' for dy in ndisps])
    assert np.isclose(min(ndisps), max(ndisps)), (min(ndisps), max(ndisps), -udl / ksoil)
    assert np.isclose(min(ndisps), -udl / ksoil, rtol=0.001), (min(ndisps), max(ndisps), -udl / ksoil)

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
    ydiff = ymin - -udl / ksoil
    max_ind = [ind0, ind1][np.argmin([yd0, yd1])]  # use the node that has the largest displacement

    # o3.integrator.DisplacementControl(osi, fd_nds[max_ind], dof=o3.cc.Y, incr=fd_incs, num_iter=10)
    o3.integrator.LoadControl(osi, incr=1 / num_incr) #better load control funciton
    
    ndr = o3.recorder.NodesToArrayCache(osi, fd_nds, dofs=[o3.cc.Y], res_type='disp') # Node displacement recorder
    efr = o3.recorder.ElementsToArrayCache(osi, fd_eles, arg_vals=['section', 1, 'force']) #element force recorder - records moments
    fr = o3.recorder.ElementsToArrayCache(osi, fd_eles, arg_vals=['section', 1, 'force']) #element force recorder - records moments
    soil_ele_forces_recorder = o3.recorder.ElementsToArrayCache(osi, sl_eles, arg_vals=['localForce']) #records the soil spring force

    ndisps = [[]]
    mom = [[]]
    shear = [[]]

    for j in range(nnodes):
        ndisps[0].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
        
    for j in range(nnodes - 1):
        mom[0].append(o3.get_ele_response(osi, fd_eles[j], 'force')[5])
        
    for j in range(nnodes-1):
        shear[0].append(o3.get_ele_response(osi, fd_eles[j], 'force')[4])

    for i in range(num_incr):
        fail = o3.analyze(osi, 1)
        if fail:
            raise ValueError()
        ndisps.append([])
        mom.append([])
        shear.append([])
        
        for j in range(nnodes):
            ndisps[-1].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
            
        for j in range(nnodes-1):
            # print('force: ', o3.get_ele_response(osi, fd_eles[j], 'force'))
            mom[-1].append(o3.get_ele_response(osi, fd_eles[j], 'force')[5])
            
        for j in range(nnodes-1):
            shear[-1].append(o3.get_ele_response(osi, fd_eles[j], 'force')[4])

        # print('y_disps: ', [f'{dy:.3}' for dy in ndisps[-1]])
        y_disp_max = o3.get_node_disp(osi, sl_nds[max_ind], dof=o3.cc.Y)
        if -y_disp_max > -ymin:
            
            # print('break: ', y_disp_max, ymin)
            break

    # for j in range(nnodes - 1):
        # print('force: ', o3.get_ele_response(osi, fd_eles[j], 'force'))
    o3.wipe(osi)  # if you are using a recorder you need to have this line
    node_disps = ndr.collect()
    forces = efr.collect()
    spring_forces = soil_ele_forces_recorder.collect()
    # print('ndisps', ndisps)
    # print("moment", mom)
    # print("shear", shear)
    # print("SPRINGZZ", spring_forces[-1], len(spring_forces))
    # print("max moment:", min(mom[-1]))
    # print("max shear:", min(shear[-1]))
    # print("min moment:", max(mom[-1]))
    # print("conc_modulus:", Ec)
    
    xd = xd
    L = xd[1] - xd[0]
    xd_list = list(xd)
    xd_list[0] = 0
    xd_list[1] = L
    print("LENGTH:", L, xd_list)
    # L = 6 #[m]
    w = udl #[kN/m]
    I = i_z
    E = Ec
    
    x_incrs = np.arange(xd_list[0]+2, xd_list[1]-2+dx, dx, dtype=None)
    print(x_incrs)
    
    M_pin = []
    V_pin = []
    disp_pin = []
    
    for x in x_incrs:
        M = -(w*x)/2*(L-x)
        M_pin.append(M)
        
        # print(M_pinned)
        
    for x in x_incrs:
        V = w*(L/2-x)
        V_pin.append(V)
        
    for x in x_incrs:
        disp = -((w*x)/(24*E*I))*(L**3-2*L*x**2+x**3)
        disp_pin.append(disp)
        
    print("max:", max(V_pin), min(M_pin), min(disp_pin))
    
    M_fix = []
    V_fix = []
    disp_fix = []
    
    for x in x_incrs:
        M = -(w/12)*((6*L*x)-(L**2)-(6*x**2))
        M_fix.append(M)
        
        # print(M_pinned)
        
    for x in x_incrs:
        V = w*(L/2-x)
        V_fix.append(V)
        
    for x in x_incrs:
        disp = -((w*x**2)/(24*E*I))*(L-x)**2
        disp_fix.append(disp)
    
    
    return nx, node_disps, forces, xd0n,xd1n, yd0, yd1, mom, shear, ndisps, spring_forces, M_fix, V_fix, disp_fix, M_pin, V_pin, disp_pin, x_incrs


def create(): #creates plot
    """
    Run an analysis imposing a uniform load, on a foundation with soil springs

    Then impose a displacement on the ends of some of the soil springs

    :param lf:
        Foundation length [m]
    :param dx:
        Target spacing of springs [m]
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
    bf, ax = plt.subplots(nrows=3, squeeze=True, sharex=(True))
    #grabs relevant values from run function
    stypes = ["Elastic"]
    cols = ['r', 'b']
       
    m_crack = -12000
    
    m_yield = -32000
    
    m_ult = -48000 #{Nm}

    support_disp = -0.1 #[m]
    
    xd = (2, 8)
    # stypes = [stypes[0]]
    for ss, stype in enumerate(stypes):
        nx, ndisps, forces, xd0n, xd1n, yd0, yd1,mom, shear, ndisps, spring_forces, M_fix, V_fix, disp_fix, M_pin, V_pin, disp_pin, x_incrs = run(lf=10, dx=0.5, xd=xd, yd=(support_disp,support_disp), ksoil=8.333e6, udl=8600, num_incr=500, stype=stype)
        # nx_rc, ndisps_rc, xd0n_rc, forces_rc, xd1n_rc, yd0_rc, yd1_rc, mom_rc, sh_rc, ndisps = run(lf=10, dx=0.5, xd=(3,7), yd=(-0.1,-0.1), ksoil=2.5e3, udl=0.03e3, axial_load=0, max_curve=0.003, num_incr=500, stype="rc")
        #ksoil=25000kN/m=> 25000e3N/m and 
        # udl = 0.45kPa => 484Pa - this is based on house weight results from henderson - scaling down the whole house weight to the proportion of the weight acting on foundation beam

        "elastic model plotting"
        # print("SPRINGGG", spring_forces)
        ax[2].plot(nx, ndisps[0], label=f'Initial {stype}', c=cols[ss], ls='-.')
        # ax[0].plot(nx, ndisps[int(len(ndisps) / 2)], label='Linear @ 50%', c="b", ls = "--") #Plotting for 50% but this is no longer 50%!! Depends on array length!!
        
        # ndisp_fin = ndisps[-1]
        # ndisp_fin[0] = 0
        # ndisp_fin[-1] = 0
        
        moment = mom[-1]
        moment[-1] = moment[0]
        
        incr = 0.5
        moment.insert(21, 0)
        
        shear_force = shear[-1]
        # shear_force.append(-shear_force[0])
        
        
        print("moment at first point:", moment[-1], len(mom[-1]))
        print("shear at first point:", shear[-1], len(shear[-1]))
        
        print("length nodes_x:", len(nx),"length node_disps:", len(ndisps[-1]),"length bending_moment:", len(mom[-1]), "length shear:",len(shear[-1]))
        
        
        
        
        ax[0].plot(nx, ndisps[-1], label=f'{stype}', c=cols[ss])
        
        # ax[2].plot(nx, py)

        # ax[0].plot([xd0n, xd1n], [yd0, yd1], c='r', label='Imposed')


        el_x = (nx[1:] + nx[:-1])/2
        print("el_x:", el_x)
        
    
        ax[1].plot(nx, moment, label=f'{stype}', c=cols[ss])
        # ax[1].plot(nx, moment, label="Elastic", c=cols[ss]) #plots
        
        
        ax[1].plot(x_incrs, M_pin, label='Pin-Pin')
        ax[1].plot(x_incrs, M_fix, label='Fix-Fix')
        
        ax[0].plot(x_incrs, disp_pin, label='Pin-Pin')
        ax[0].plot(x_incrs, disp_fix, label='Fix-Fix')
        
        ax[2].plot(x_incrs, V_pin, label='Pin')
        ax[2].plot(x_incrs, V_fix, label='Fix')

        
        print("max_moment:", min(mom[-1]))
        print("max_disp:", min(ndisps[-1]))
        print("max shear:", min(shear[-1]))
        
        fmt = "x"
        # ax[2].plot(el_x, py, label=f'Initial {stype}', c=cols[ss], ls='-.')

        ax[2].plot(nx, spring_forces[-1],fmt ,label="Node location", c=cols[ss])
        # ax[2].plot(nx, spring_forces[-1], fmt)
        
        # ax[2].plot(nx, spring)
        ax[2].plot(el_x, shear[-1], label=f'{stype}', c=cols[ss])
        
        # ax[2].plot(el_x, shear_force, label=f'{stype}', c=cols[ss], drawstyle="steps-mid")
        # ax[2].plot(el_x, shear_force, label=f'{stype}', alpha=1, c = "grey")
        
        # plt.step(x, y + 1, where='mid', label='mid')
        # ax[1].plot(nx, mom[-1], label=f'{stype}', c=cols[ss])
        # ax[2].plot(nx, shear_force, label=f'{stype}', c=cols[ss])
        # ax[2].plot(nx, shear[-1], label=f'{stype}', c=cols[ss])


    ax[0].axvspan(*xd, color=(0.3, 0.3, 0.3, 0.5), zorder=0)
    ax[1].axvspan(*xd, color=(0.3, 0.3, 0.3, 0.5), zorder=0)
    ax[2].axvspan(*xd, color=(0.1, 0.1, 0.1, 0.2), zorder=0)
    # ax[0].plot([xd0n_elastic, xd1n_elastic], [yd0_elastic, yd1_elastic], c='r', label='Imposed') #plots imposed displacement
    # ax[1].axhline(m_crack, ls='-.', label="Cracking Moment", c="g")
    # ax[1].axhline(m_yield, ls="-.", label='Yield Moment', c="orange")
    # ax[1].axhline(m_ult, ls="-.", label="Ultimate Moment", c="k")
    
    # ax[0].set_xlabel('Foundation Length (m)')
    ax[0].set_ylabel('Displacement (m)')
    ax[1].set_ylabel("Moment (Nm)")
    ax[2].set_xlabel("Foundation Beam Length (m)")
    ax[2].set_ylabel("Spring force")
    # ax[2].set_ylabel("Force (N)")
    # plt.grid()
    
    ax[1].legend(bbox_to_anchor =(1,1))
    ax[0].legend(bbox_to_anchor =(1,1))
    ax[2].legend(bbox_to_anchor =(1,1))
    
    # ax[2].legend(bbox_to_anchor=(1,1))
    
    ax[0].grid()
    ax[1].grid()
    ax[2].grid()
    # ax[2].grid()
    

    plt.show()





if __name__ == '__main__':

    "Impose foundation displacement function"
    create()