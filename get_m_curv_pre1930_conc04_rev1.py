# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 10:25:48 2022

@author: mda153
"""

import o3seespy as o3
import numpy as np
import o3seespy.extensions


def get_moment_curvature(axial_load=0, max_curve=0.001, num_incr=500):
    osi = o3.OpenSeesInstance(ndm=2, ndf=3, state=3)

    f_c28 = 2.466  #28 day compressive strength [ksi]
    e_mod = 57000.0 * np.sqrt(f_c28 * 1000.0) / 1e3 #Is used here in the Concrete 04 model => Concrete 01 does not account for this (concrete modulus)
    f_conc_tens = 0.36 * np.sqrt(f_c28)  # NZS3101 CL 5.2.6, concrte tensile strength
    eps_conc_tens = f_conc_tens / e_mod  # Park and Paulay 1975 assume linear
    eps_unconf_conc_max = 0.004  #Ultimate strain for unconfined concrete Priestly et al. 2007, lower bound (0.004 - 0.005)

    eps_conc_crush = 0.03  #Concrte strain at crushing - Koopaee (2015) from NZS 3101?? taken from Fig 2.1 fc=20MPa curve
    
    

    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-f_c28, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush, ec=e_mod, fct=f_conc_tens, et=eps_conc_tens) #unconfined concrte properties
 
    #Defining concrete beam parameters

    h = 23.622 #Section height
    b = 5.906 #Section width

    gj = 1.0E10 # linear elastic torsional stiffness
   
    nf_y = 2 #number of fibers in Y-Y dir
    nf_z = 2 #number of fibers in Z-Z dir
 
    
    
#Defining beam section
    edge_y = h / 2.0
    edge_z = b / 2.0



    sect = o3.section.Fiber(osi, gj=gj) #defining a single fibre
    # define the core patch
    #Patch command generates a number of fibers over a geometric cross-section
    #patch.Quad(material tag, number of fibers in Z-Z dir, number of fibers in Y-Y dir)

    o3.patch.Rect(osi, conc_unconf, nf_z, nf_y,  # defining square rectangular cross section
                  crds_i=[-edge_y, -edge_z],
                  crds_j=[edge_y, edge_z])
    
    
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
    return moment, curvature


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    mom, curve = get_moment_curvature()
    plt.plot(curve, mom)
    axes = plt.axes()
    axes.set_xlabel("Curvature [1/in^2]")
    axes.set_ylabel("Bending Moment [kip*in]")
    
    

    plt.show()