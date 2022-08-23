# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 12:55:24 2022

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


def run(lf=10, dx=0.5, s_depth=0.6, xd=(4,6), yd=(-0.25,-0.25), ksoil=3.846e6, udl=8600):
    """
    Run an analysis imposing a uniform load and a target displacement, on a foundation with soil springs

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
        Uniform distributed load on foundation [Pa] - based on house load - Henderson then multiplied by area ratio of founation to total house area
    :param xd:
        (x0, x1) positions of displaced section of soil
    :param yd:
        (y0, y1) displacements of section of soil at x0 and x1 (note should be -ve)
    :return:
    """

    osi = o3.OpenSeesInstance(ndm=2, ndf=3)  # 2D with 3 DOF
    s_width = 0.13  # [m] section width
    imposed_displacement = 0.1

    nnodes = max(int(lf / dx) + 1, 2)  # number of nodes
    nx = np.linspace(0, lf, nnodes)  # x-position of nodes
    x_trib = np.ones_like(nx)  # Trib area
    x_trib[1:] += np.diff(nx)
    x_trib[:-1] += np.diff(nx)
    
    # print("Hi xtrib", x_trib)
    
    # x_trib = int(0.5, 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
    
    nnodes = max(int(lf / dx) + 1, 2)  # number of nodes
    nx = np.linspace(0, lf, nnodes)  # x-position of nodes
    spring_spacing = np.ones_like(nx)  # Trib area
    spring_spacing[1:] += np.diff(nx)
    spring_spacing[:-1] += np.diff(nx)
    # spring_spacing = (spring_spacing)/4
    
    # print("spring spacing", spring_spacing)

    ks = ksoil * x_trib * s_width # N/m (TODO: define this better, see Henderson thesis for details)
    

    py = ks * imposed_displacement  # N Capacity (TODO: define this better, see Henderson thesis for details)
    
    #Could not find this in henderson thesis - does not expicitly mention it
    #Surely the force resistance should be related to the stiffness of each spring then multiply this by the displacement it is subject to
    
    print("ks", ks, "py", py)

    transf = o3.geom_transf.Linear2D(osi, [])  # Using simple transform - no P-delta effects

    fd_nds = []  # foundation nodes
    sl_nds = []  # soil nodes
    sl_eles = [] # soil elemeents
    fd_eles = [] # foundation elements
    for i in range(len(nx)):
        fd_nds.append(o3.node.Node(osi, nx[i], 0.0))
        sl_nds.append(o3.node.Node(osi, nx[i], 0.0))
        o3.Mass(osi, fd_nds[-1], 1.0, 1.0, 1.0)
        dettach = 1
        mat_base = o3.uniaxial_material.ElasticPP(osi, 1 * ks[i], 1 * py[i] / ks[i],
                                                  -1 * py[i] / ks[i])  # low tension stiffness
        if dettach:
            mat_obj2 = o3.uniaxial_material.Elastic(osi, 1000 * ks[i], eneg=0.001 * ks[i])
            # mat_obj2 = o3.uniaxial_material.Elastic(osi, 0.001 * ks[i], eneg=1000 * ks[i])
            mat = o3.uniaxial_material.Series(osi, [mat_base, mat_obj2])
        else:
            mat = o3.uniaxial_material.ElasticPP(osi, ks[i], 1 * py[i] / ks[i], -1 * py[i] / ks[i])
        # mat = o3.uniaxial_material.ElasticPP(osi, ks[i], 1*py[i] / ks[i], -0.001 * py[i] / ks[i]) #remove tension capacity
        #mat = o3.uniaxial_material.Elastic(osi, ks[i])
        sl_eles.append(o3.element.ZeroLength(osi, [fd_nds[-1], sl_nds[-1]], mats=[mat], dirs=[o3.cc.Y]))
        if i != 0:
            fc = 17.2e6 #[Pa] concrete compressive strength
            e_mod = 5000*np.sqrt(fc)*1e3
            print('e_mod', e_mod)# [Pa] Elastic modulus of concrete
            # e_mod = 26e6
            area = s_width * s_depth
            i_z = s_depth ** 3 * s_width / 12
            fd_eles.append(o3.element.ElasticBeamColumn2D(osi, [fd_nds[i], fd_nds[i-1]], area=area, e_mod=e_mod, iz=i_z, transf=transf))

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
    assert np.isclose(min(ndisps), -udl / ksoil, rtol=0.001), (min(ndisps), max(ndisps), -udl / ksoil)

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
    o3.integrator.DisplacementControl(osi, fd_nds[max_ind], dof=o3.cc.Y, incr=ymin / 100, num_iter=100)
    print("ymin", ymin)
    
    ndr = o3.recorder.NodesToArrayCache(osi, fd_nds, dofs=[o3.cc.Y], res_type='disp') # Node displacement recorder
    efr = o3.recorder.ElementsToArrayCache(osi, fd_eles, arg_vals=['section', 1, 'force']) #element force recorder -using this to record moment
    
    mom = [[]]
    
    ndisps = [[]]
    for j in range(nnodes):
        ndisps[0].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
        
    for j in range(nnodes - 1):
        mom[0].append(o3.get_ele_response(osi, fd_eles[j], 'force')[5])
        
    for i in range(100):
        fail = o3.analyze(osi, 1)
        if fail:
            raise ValueError()
        ndisps.append([])
        mom.append([])
        for j in range(nnodes):
            ndisps[-1].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))
        for j in range(nnodes-1):
            # print('force: ', o3.get_ele_response(osi, fd_eles[j], 'force'))
            mom[-1].append(o3.get_ele_response(osi, fd_eles[j], 'force')[5])
            
        print('mom', mom)
        print('y_disps: ', [f'{dy:.3}' for dy in ndisps[-1]])
        y_disp_max = o3.get_node_disp(osi, sl_nds[max_ind], dof=o3.cc.Y)
        if -y_disp_max > -ymin:
            break
        

    o3.wipe(osi)  # if you are using a recorder you need to have this line
    node_disps = ndr.collect()
    forces = efr.collect()
    
    m_henderson = -17806
    
    # print('forces',mom, len(forces))
    
    bf, ax = plt.subplots(nrows=2, squeeze=True)
    
    # el_x = (nx_elastic[1:] + nx_elastic[:-1])/2
    
    el_x = (nx[1:] + nx[:-1])/2 #averages nodes to plot moment - elements in x
    
    # len_found = ([0:10:0.1])
    

    ax[0].plot(nx, ndisps[0], label='Initial Foundation')
    # ax[0].plot(nx, ndisps[int(len(ndisps) / 2)], label='Foundation @ 50%')
    ax[0].plot(nx, ndisps[-1], label='Foundation @ 100%', c='k')
    # ax[0].plot([xd0n, xd1n], [yd0, yd1], c='k', label='Imposed') #plots imposed displacement
    ax[1].plot(el_x, mom[-1], label="Demand")
    
    ax[1].axhline(m_henderson , ls="-.", label="Henderson SAP2000 Elastic Model", c="k")
    

    
        
    ax[0].set_xlabel('Foundation Length (m)')
    ax[0].set_ylabel('Vertical Disp. (m)')
    ax[1].set_ylabel("Bending Moment (Nm)")
    ax[1].set_xlabel("Foundation Length (m)")
    
    ax[1].legend(bbox_to_anchor =(1,0))
    ax[0].legend(bbox_to_anchor =(1,1))
    
    ax[0].grid()
    ax[1].grid()
        
    plt.legend()
    plt.show()
    


    
 
    



if __name__ == '__main__':
    run(10, 0.5, s_depth=0.6)