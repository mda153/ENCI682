# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import o3seespy as o3
import numpy as np
import o3seespy.extensions


def get_moment_curvature(axial_load=0, max_curve=0.001, num_incr=500):
    osi = o3.OpenSeesInstance(ndm=2, ndf=3, state=3)


    f_c28 = 2.569  #28 day compressive strength [ksi]
    e_mod = 57000.0 * np.sqrt(f_c28 * 1000.0) / 1e3 #Is used here in the Concrete 04 model => Concrete 01 does not account for this (concrete modulus)
    f_conc_tens = 0.36 * np.sqrt(f_c28)  # NZS3101 CL 5.2.6, concrte tensile strength
    eps_conc_tens = f_conc_tens / e_mod  # Park and Paulay 1975 assume linear
    eps_unconf_conc_max = 0.004  #Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)
    eps_conf_conc_max = 0.06 #Ultimate strain for confined concrete Preistly et al. 2007, lower bound (0.06 - 0.09)
    eps_conc_crush = 0.03  #Concrte strain at crushing - Koopaee (2015) from NZS 31010?? taken from Fig 2.1 fc=20MPa curve
    
    
    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-f_c28, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush, ec=e_mod, fct=f_conc_tens, et=eps_conc_tens) #confined concrete properties
    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-f_c28, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush, ec=e_mod, fct=f_conc_tens, et=eps_conc_tens) #unconfined concrte properties
    rebar = o3.uniaxial_material.Steel01(osi, fy=43.511, e0=29732.736, b=0.02) #Steel reinforcing properties
    
    #Defining concrete beam parameters

    h = 23.622 #Section height
    b = 5.906 #Section width
    cover = 1.732 #Cover depth
    gj = 1.0E10 # Is this G * J?
    nf_core_y = 8
    nf_core_z = 8
    nf_cover_y = 10
    nf_cover_z = 10
    n_bars = 2 #Number of bars
    bar_area = 0.351 #Area of steel (inches^2)

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
    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[core_y, core_z], end=[core_y, -core_z])

    spacing_y = 2 * core_y / (n_bars - 1)
    remaining_bars = n_bars - 1
    o3.layer.Straight(osi, rebar, remaining_bars, bar_area,
                      start=[core_y - spacing_y, core_z],
                      end=[-core_y + spacing_y, core_z])
    o3.layer.Straight(osi, rebar, remaining_bars, bar_area,
                      start=[core_y - spacing_y, -core_z],
                      end=[-core_y + spacing_y, -core_z])

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
    return moment, curvature


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    mom, curve = get_moment_curvature()
    plt.plot(curve, mom)
    axes = plt.axes()
    axes.set_xlabel("Curvature [1/in^2]")
    axes.set_ylabel("Bending Moment [kip*in]")
    
    

    plt.show()
