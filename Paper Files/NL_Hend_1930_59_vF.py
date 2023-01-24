# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:10:23 2022

@author: mda153
"""
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
    fc = 17.2e6 #[Pa] concrete compressive strength
    fy = 300.0e6  #[Pa] steel yield strength
    Es = 200.0e9  # [Pa] Steel modulus of elasticity 
    depth = 0.6  #[m] section depth
    width = 0.15 #[m] section width
    
    


    col_sect = sm.sections.RCDetailedSection(depth=depth, width=width)
    col_sect.layer_depths = [0.044, 0.556]  # 44mm cover (m units)
    col_sect.bar_diams = [[0.012, 0.012], [0.012, 0.012]]  # 12mm bars - two 12mm bars directly beside each other at top of section
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


    imposed_displacement = 0.1 #[m]
    d = 0.6  # section depth [m]
    b = 0.15  # section width [m]

    # Mander inputs - also used in opensees inputs - From View Material Behavior Script
    fpc = 17.2e6  # Conc compressive strength (Pa)

    Z = b*d**2/6
    ft = 0.56*np.sqrt(fpc)*1e3 #CUMBIA tensile strenght of conc Montejo et al. 2007 script for auto calc
    m_crack = Z*ft

    Ec = 5000 * np.sqrt(fpc)*1e3  # Concrete modulus of elasticity (Pa) - eqn in Mander script for auto calc
    
    eps_unconf_conc_max = 0.0064 # Ultimate strain for unconfined concrete
    eps_conf_conc_max = 0.015 #Ultimate strain for conf. concrete
    eps_conc_crush_unconf = 0.004 # Concrte strain at crushing unconf.
    eps_conc_crush_conf = 0.004 #conf conc strain at crushing conf.

    fy = 300.0e6  #Pa
    Es = 200.0e9  # Steel modulus [Pa]
    
    nf_core_y = 16  # number of fibers in Y-Y dir - for the concrete core - (10 - 20 is a good number here)
    nf_core_z = 16  # number of fibers in Z-Z dir
    nf_cover_y = 20  # number of fibers Y-Y dir, including cover concrete
    nf_cover_z = 20  # number of fibers Z-Z dir, including cover concrete


    osi = o3.OpenSeesInstance(ndm=2, ndf=3)  # 2D with 3 DOF
    s_width = b # section width [m]
    s_depth = d #section depth [m]
    
    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf, ec=Ec, fct=ft ,et=0) #confined concrete properties
 
    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush_unconf, ec=Ec, fct=ft, et=0) #unconfined concrte properties
 
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.01) #Steel reinforcing properties

    nnodes = max(int(lf / dx) + 1, 2)  # number of nodes
    nx = np.linspace(0, lf, nnodes)  # x-position of nodes
    x_trib = np.zeros_like(nx)  # Trib area
    x_trib[1:] += np.diff(nx) / 2
    x_trib[:-1] += np.diff(nx) / 2

    ks = ksoil * s_width *x_trib  #soil spring stiffness N/m
    
    py = ks * imposed_displacement*1e3  # N Capacity (TODO: define this better, see Henderson thesis for details) spring force resistance


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
            mat = o3.uniaxial_material.Series(osi, [mat_base, mat_obj2])
        else:
            mat = o3.uniaxial_material.ElasticPP(osi, ks[i], 1*py[i] / ks[i], -1 * py[i] / ks[i]) #remove tension capacity
        sl_eles.append(o3.element.ZeroLength(osi, [fd_nds[-1], sl_nds[-1]], mats=[mat], dirs=[o3.cc.Y]))
        "Foundation nodes"
        if i != 0:
            fc = 17.2e6 #[Pa]
            e_mod = 5000*np.sqrt(fc)*1e3 # Pa
            
            area = s_width * s_depth
            i_z = s_depth ** 3 * s_width / 12

            if stype == "Elastic":
                bot_sect = o3.section.Elastic2D(osi, e_mod, area, i_z)
                integ = o3.beam_integration.Lobatto(osi, bot_sect, big_n=5)

            elif stype == "Non-Linear":
                bot_sect = get_rc_fibre_section(osi, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z)
                integ = o3.beam_integration.Lobatto(osi, bot_sect, big_n=5)
            else:
                raise ValueError('set stype to Elastic or Non-Linear')
            fd_eles.append(o3.element.DispBeamColumn(osi, [fd_nds[i-1], fd_nds[i]], transf=transf, integration=integ))

    o3.Fix3DOFMulti(osi, sl_nds, o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)  # Fix all the soil nodes
    o3.Fix3DOF(osi, fd_nds[0], o3.cc.FIXED, o3.cc.FREE, o3.cc.FREE)

    # Define load
    ploads = udl * x_trib
    ts0 = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts0)
    for i in range(nnodes):
        o3.Load(osi, fd_nds[i], [0, -ploads[i], 0])

    # Run initial static analysis
 
    o3.constraints.Transformation(osi)
    # o3.test_check.NormDispIncr(osi, tol=1.0e-6, max_iter=35, p_flag=0)
    o3.test.NormUnbalance(osi, tol=1.0e-4, max_iter=35, p_flag=0)
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
    assert np.isclose(min(ndisps), max(ndisps)), (min(ndisps), max(ndisps), -udl / (ksoil * b))
    assert np.isclose(min(ndisps), -udl / (ksoil * b), rtol=0.001), (min(ndisps), max(ndisps), -udl / (ksoil * b))

    total_force = sum([o3.get_ele_response(osi, ele, 'localForce')[0] for ele in sl_eles])
    # print('total_force under initial load: ', total_force)
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
        mom.append([-o3.get_ele_response(osi, fd_eles[0], 'force')[2]])
        shear.append([])

        for j in range(nnodes):
            ndisps[-1].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
            
        for j in range(nnodes-1):
            forces = o3.get_ele_response(osi, fd_eles[j], 'force')

            mom[-1].append(-forces[5])
            shear[-1].append(forces[4])


        y_disp_max = o3.get_node_disp(osi, sl_nds[max_ind], dof=o3.cc.Y)
        if -y_disp_max > -ymin:
            
     
            break
    total_force = sum([o3.get_ele_response(osi, ele, 'localForce')[0] for ele in sl_eles])

    o3.wipe(osi)  # if you are using a recorder you need to have this line
    node_disps = ndr.collect()
    forces = efr.collect()
    spring_forces = soil_ele_forces_recorder.collect()
    
    return nx, node_disps, forces, xd0n,xd1n, yd0, yd1, mom, shear, ndisps, spring_forces, m_crack


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
        Uniform distributed load along foundation beam [N/m]
    :param xd:
        (x0, x1) positions of displaced section of soil
    :param yd:
        (y0, y1) displacements of section of soil at x0 and x1 (note should be -ve)
    :return:
    """
    bf, ax = plt.subplots(nrows=3, squeeze=True, sharex='col', figsize=(6, 8))
    #grabs relevant values from run function
    stypes = ["Elastic", "Non-Linear"]
    cols = ['r', 'b']

    support_disp = -0.1 #[m]
    udl = 6900 # kN/m
    lf = 15  # length of foundation
    dx = 0.1 #spring spacing
    xd = (0, 15)
    i = 20 #loss of support [m]
    los = dx * i * 2  # in increments of 2 since spring in centre
    xd = (lf / 2 - los / 2, lf / 2 + los / 2) #This here controls where the LOS starts and finshes [m] along foundation
    # stypes = [stypes[0]]
    for ss, stype in enumerate(stypes):
        nx, node_disps, forces, xd0n, xd1n, yd0, yd1,mom, shear, ndisps, spring_forces, m_crack = run(lf=lf, dx=dx, xd=xd, yd=(support_disp,support_disp), ksoil=2.5*1E7, udl=udl, num_incr=1000, stype=stype)


        "elastic model plotting"
        ax[2].plot(nx, ndisps[0], label=f'Initial {stype}', c=cols[ss], ls='-.')

        
        moment = mom[-1]
        moment[-1] = moment[0]
        


        ind = np.argmin(abs(nx - xd[0]))
        diff_disp = (node_disps[-1] - node_disps[0])
        ax[0].plot(nx, diff_disp - diff_disp[ind], label=f'{stype}', c=cols[ss])

        max_sag_disp = abs(min(diff_disp - diff_disp[ind]))
        
        
        LOS_length = 4

        el_x = (nx[1:] + nx[:-1])/2

        ax[1].plot(nx, moment, label=f'{stype}', c=cols[ss])
 
        
        max_disps = np.array(ndisps[-1])
        
        # slope = abs(max_disps[7] - max_disps[6])/(nx[7] - nx[6])
        # print("The Maximum Slope Is:", slope)
        
        # print("max_vert_disp:", min(max_disps), i/2)
        
        LOS_length = 4
        
        drift = (abs(min(max_disps)))/(LOS_length/2)*100

        print("Maximum saging displacement:", max_sag_disp)
        print("max positive moment:", min(mom[-1]),"max_negative_moment", max(mom[-1]))
        print("max and min shear:", min(shear[-1]), max(shear[-1]))
        print("Maximum Drift Demand Is:", drift)

        ax[2].plot(el_x, shear[-1], label=f'{stype}', c=cols[ss])


    ax[0].axvspan(*xd, color=(0.3, 0.3, 0.3, 0.5), zorder=0)
    ax[1].axvspan(*xd, color=(0.3, 0.3, 0.3, 0.5), zorder=0)
    ax[2].axvspan(*xd, color=(0.1, 0.1, 0.1, 0.2), zorder=0)
    # ax[0].plot([xd0n_elastic, xd1n_elastic], [yd0_elastic, yd1_elastic], c='r', label='Imposed') #plots imposed displacement
    fc = 17.2e6
    ft = 0.56*(fc**0.5)*1e3
    b = 0.15
    h = 0.6
    
    Z = (b*h**2)/6
    m_crack = ft*Z
    # print("m_crack:", m_crack*1000)
    

    h = 0.6
    d_bar = 0.012
    n_bar_total = 2 #in tension
    cover = 0.044
    As = np.pi*(d_bar/2)**2*n_bar_total
    fy = 300e6
    j = 0.85
    d_eff = cover - h
    
    m_y = As*fy*j*d_eff
    # print('m_y estimate:', m_y)
    
    m_ult = 55000
    
    ax[1].axhline(m_crack, ls='-.', label="Cracking Moment", c="g")
    ax[1].axhline(-m_crack, ls='-.', label="Cracking Moment", c="g")
    ax[1].axhline(m_y, ls="-.", label='Yield Moment', c="orange")
    ax[1].axhline(-m_y, ls="-.", label='Yield Moment', c="orange")
    ax[1].axhline(-m_ult, ls="-.", label="Ultimate Moment", c="k")
    
    # ax[0].set_xlabel('Foundation Length (m)')
    ax[0].set_ylabel('Displacement (m)')
    ax[1].set_ylabel("Moment (Nm)")
    ax[2].set_xlabel("Foundation Beam Length (m)")
    ax[2].set_ylabel("Shear Force (N)")
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
