
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import o3seespy as o3
import numpy as np
import o3seespy.extensions


def get_moment_curvature(axial_load=0, max_curve=0.2, num_incr=1000):
    osi = o3.OpenSeesInstance(ndm=2, ndf=3, state=3)
    
    fpc = 17.2e6 #Conc compressive strength (Pa)
    eps_unconf_conc_max = 0.0064 # Ultimate strain for unconfined concrete
    eps_conf_conc_max = 0.015 #Ultimate strain for conf. concrete
    eps_conc_crush_unconf = 0.004 # Concrte strain at crushing unconf.
    eps_conc_crush_conf = 0.004 #conf conc strain at crushing conf.
    Ec = 5000*(fpc**0.5)*1e3 #Conc modulus of elasticity (Pa) in Montejo et al. 2007 script for auto calc
    print("Estimated concrete modulus:",Ec)
    d_bar = 0.012 # longtitudinal bar diameter [m]
    total_num_bars = 4
    Ast = (np.pi*(d_bar/2)**2)*total_num_bars
    # print(Ast)#Total area of longtitudinal steel (m^2)
    cover = 0.044 #cover to longtitudinal bars (m)
    fy = 300e6 #Pa
    Es = 200e9 #Steel modulus (Pa)

    h = 0.6 #section height (mm)
    b = 0.15 #section width (mm)
  
    ft = 0.56*np.sqrt(fpc)*1e3 #CUMBIA tensile strenght of conc Montejo et al. 2007 script for auto calc
    
    Z = (b*h**2)/6
    
    mom_crack = Z*ft
    print("estimated cracking moment:", mom_crack/1000)

    print("estimated concrete tensile strength:", ft)
      
    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf, ec=Ec, fct=ft ,et=0) #confined concrete properties
 
    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush_unconf, ec=Ec, fct=ft, et=0) #unconfined concrte properties
 
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.01) #Steel reinforcing properties
    
    gj = 1.0E10 # Is this G * J? 
    nf_core_y = 16
    nf_core_z = 16
    nf_cover_y = 20
    nf_cover_z = 20
    
    n_bars = 2 #Number of bars in a layer/row
    bars_total = 4
    bar_area = Ast/bars_total #Area of steel (m^2)
    
    j = 0.85 #Assumed
    d_eff = h - cover
    mom_yield = (n_bars*bar_area*fy*j*d_eff)/1000
    print("estimated_yeild_moment:", mom_yield)
        
    print("Bar_area (single bar):", bar_area, "Number of bars per layer", n_bars)

    edge_y = h / 2.0
    edge_z = b / 2.0
    
    core_y = edge_y - cover
    core_z = edge_z - cover


    sect = o3.section.Fiber(osi, gj=gj)
    # define the core fibre element patch
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
    
    #Assigns rebar layers
    
    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[-core_y+d_bar/2, core_z], end=[-core_y+d_bar/2, -core_z])
    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[core_y-d_bar/2, core_z], end=[core_y-d_bar/2, -core_z]) #Assigns rebar layers
    
    #define element and nodes

    n1 = o3.node.Node(osi, 0.0, 0.0)
    n2 = o3.node.Node(osi, 0.0, 1.0)

    o3.Fix3DOF(osi, n1, 1, 1, 1)
    o3.Fix3DOF(osi, n2, 0, 0, 0)
    transf = o3.geom_transf.Linear2D(osi, [])
    integ = o3.beam_integration.Lobatto(osi, sect, big_n=5)
    ele = o3.element.DispBeamColumn(osi, [n1, n2], transf=transf, integration=integ)
    
#Recording analysis results at nodes
    ndr = o3.recorder.NodeToArrayCache(osi, n2, dofs=[3], res_type='disp')
    nmr = o3.recorder.NodeToArrayCache(osi, n1, dofs=[3], res_type='reaction')
    ecr_comp = o3.recorder.ElementToArrayCache(osi, ele, arg_vals=['section', 'fiber', h/2, 0, 'stressStrain']) #Recording stress-Strain behaviour compression fiber
    ecr_tens = o3.recorder.ElementToArrayCache(osi, ele, arg_vals=['section', 'fiber', -h/2, 0, 'stressStrain']) #Recording stress-Strain behaviour tension fiber
    ecr_steel = o3.recorder.ElementToArrayCache(osi, ele, arg_vals=['section', 'fiber', -h/2, 0, 2,'stressStrain']) #Recording stress-Strain behaviour steel fiber
    
    

#Define constant axial load
    ts = o3.time_series.Constant(osi)
    o3.pattern.Plain(osi, ts)
    o3.Load(osi, n2, load_values=[axial_load, 0.0, 0.0])
    
    
#Define analysis parameters
    o3.system.BandGeneral(osi)
    o3.numberer.Plain(osi)
    o3.constraints.Transformation(osi)
    o3.test.NormUnbalance(osi, tol=1.0e-9, max_iter=10)
    o3.algorithm.Newton(osi)
    o3.integrator.LoadControl(osi, incr=0.0)
    o3.analysis.Static(osi)
    o3.analyze(osi, 1)

#Define reference moment ...
    ts = o3.time_series.Linear(osi)
    o3.pattern.Plain(osi, ts)
#Computing curvaure increments
    d_cur = max_curve / num_incr
    #
    variable_disp_cont = 0
    if variable_disp_cont:
        o3.Load(osi, n2, load_values=[0, 0, 1])
        o3.integrator.DisplacementControl(osi, n2, o3.cc.DOF2D_ROTZ, d_cur, 1, d_cur, d_cur)
    else:
        o3.SP(osi, n2, o3.cc.DOF2D_ROTZ, dof_values=[d_cur])
        o3.integrator.LoadControl(osi, incr=1)
    o3.analyze(osi, num_incr)
    o3.wipe(osi)
    rot = ndr.collect()
    moment = -nmr.collect()

    #Estimate yield curvature
    
    rho_long = Ast/(b*d_eff) #Longtitudinal reinforcement ratio
    n = Es/Ec #Modular ratio
    k = np.sqrt((rho_long**2*n**2) + (2*rho_long*n)) - (rho_long*n) #Ratio of depth to N.A., Pujol et. al.
    esp_yeild = fy/Es #Yeild strain
    curv_yeild = esp_yeild/((1-k)*d_eff) #Estimated yeild curvature
    print("Estimated yeild curvature:", curv_yeild)

    return moment, rot, mom_crack, mom_yield


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    mom, curve ,mom_crack, mom_yield = get_moment_curvature(num_incr=100, max_curve=0.2)

    plt.plot(curve, mom/1000)

    curve_list = curve.tolist()
    mom_list = mom.tolist()

    
    np.savetxt('m_curv_opensees_vTEST_1930_59_Hend.csv', (mom_list, curve_list), fmt='%.4e', delimiter=',')# X is an array

    # fpc = 17.2e6  # Conc compressive strength (Pa)
    # Ec = 5000 * (fpc ** 0.5) * 1E3
    # d = 0.6  # section depth (m)
    # b = 0.15  # section width (m)
    # ei_el = Ec * b * d ** 3 / 12/1000
    # moms = np.array([0, 60])
    # curvs = moms / ei_el
    # plt.plot(curvs, moms, c='k')
    # print(curvs, moms)

    plt.grid()
    
    axes = plt.gca()
    axes.set_xlabel("Curvature [1/m]")
    axes.set_ylabel("Bending Moment [kNm]")
    # axes.set_xlim([0, 0])
    # axes.set_xticks([0.0001])
    axes.legend()
    
    

    plt.show()