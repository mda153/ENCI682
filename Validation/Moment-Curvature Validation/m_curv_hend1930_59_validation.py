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
    Ec = 5000*(fpc**0.5)*1e3 #Conc modulus of elasticity (Pa) in Mander script for auto calc
    print("EC:",Ec)
    d_bar = 0.012 # longtitudinal bar diameter [m]
    total_num_bars = 4
    Ast = np.pi*(d_bar/2)**2*total_num_bars
    # print(Ast)#Total area of longtitudinal steel (m^2)
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
    
    fpcu = eps_conc_crush_conf*Ec
    # fpcu = 50*fpc
    print("CUMBIA Conc Crush Strength:", fpcu)

    

    
    # ft = 0.8*np.sqrt(fpc)*1e3
    ft = 0.56*np.sqrt(fpc)*1e3 #CUMBIA tensile strenght of conc
    # ft = 1.4e3
    
    Z = (b*d**2)/6
    
    mom_crack = Z*ft
    print("mom_crack:", mom_crack/1000)
    
    n_bars_tens = 2
    bars_total = 4
    bar_area = Ast/bars_total
    print("Bar area:", bar_area)
    j = 0.85 #Assumed
    d_eff = d - clb
    

    print("ft:", ft)
    
    
    mom_yield = (n_bars_tens*bar_area*fy*j*d_eff)/1000
    print("mom_yeild:", mom_yield)
    
    conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf, ec=Ec, fct=ft ,et=0) #confined concrete properties
    # conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-fpc, epsc0=-eps_conf_conc_max, fpcu=-fpcu, eps_u=-eps_conc_crush_conf)
    conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush_unconf, ec=Ec, fct=ft, et=0) #unconfined concrte properties
    # conc_unconf = o3.uniaxial_material.Concrete01(osi, fpc=-fpc, epsc0=-eps_unconf_conc_max, fpcu=-0, eps_u=-eps_conc_crush_unconf)
    rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.00001) #Steel reinforcing properties
    # rebar = o3.uniaxial_material.Steel02(osi, fy=fy, e0=Es, b=0.0001, params=[3.5, 0.005, 0.005]) #Steel reinforcing properties
    
    # conc_conf = rebar
    # conc_unconf = rebar
    # rebar = conc_conf
    #Defining concrete beam parameters

    h = d #Section height
    b = b #Section width
    cover = clb #Cover depth
    gj = 1.0E10 # Is this G * J? 
    nf_core_y = 16
    nf_core_z = 16
    nf_cover_y = 20
    nf_cover_z = 20
    n_bars = 2 #Number of bars in a layer/row
    bars_total = 4
    bar_area = Ast/bars_total #Area of steel (m^2)
    print(bar_area)


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

    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[-core_y, core_z], end=[-core_y, -core_z])
    o3.layer.Straight(osi, rebar, n_bars, bar_area, start=[core_y, core_z], end=[core_y, -core_z]) #Assigns rebar layers

    # spacing_y = 2 * core_y / (n_bars - 1)
    spacing_y = 2 * core_y / (n_bars - 1)
    # print(spacing_y)
    spacing_y = 0.062 #[m] spacing between reinforcing bars (horizontal)
    print(spacing_y)
    remaining_bars = n_bars - 1 #question here
    o3.layer.Straight(osi, rebar, remaining_bars, bar_area,
                      start=[core_y - spacing_y, core_z],
                      end=[-core_y + spacing_y, core_z])
    o3.layer.Straight(osi, rebar, remaining_bars, bar_area,
                      start=[core_y - spacing_y, -core_z],
                      end=[-core_y + spacing_y, -core_z])
    
    #define element and nodes

    n1 = o3.node.Node(osi, 0.0, 0.0)
    n2 = o3.node.Node(osi, 0.0, 1.0)
    # n3 = o3.node.Node(osi, 0.0, 1)
    o3.Fix3DOF(osi, n1, 1, 1, 1)
    o3.Fix3DOF(osi, n2, 0, 0, 0)
    transf = o3.geom_transf.Linear2D(osi, [])
    integ = o3.beam_integration.Lobatto(osi, sect, big_n=5)
    ele = o3.element.DispBeamColumn(osi, [n1, n2], transf=transf, integration=integ)
    # ele = [1, 2, 3, 5]
    print(ele)
    # ele = o3.element.DispBeamColumn(osi, [n1, n3], transf=transf, integration=integ)
    # ele = []

    # ele = o3.element.ZeroLengthSection(osi, [n1, n2], sect)
    
#Recording analysis results at nodes
    ndr = o3.recorder.NodeToArrayCache(osi, n2, dofs=[3], res_type='disp')
    nmr = o3.recorder.NodeToArrayCache(osi, n1, dofs=[3], res_type='reaction')
    ecr_comp = o3.recorder.ElementToArrayCache(osi, ele, arg_vals=['section', 'fiber', h/2, 0, 'stressStrain']) #Recording stress-Strain behaviour compression fiber
    ecr_tens = o3.recorder.ElementToArrayCache(osi, ele, arg_vals=['section', 'fiber', -h/2, 0, 'stressStrain']) #Recording stress-Strain behaviour tension fiber
    ecr_steel = o3.recorder.ElementToArrayCache(osi, ele, arg_vals=['section', 'fiber', -h/2, 0, 2,'stressStrain']) #Recording stress-Strain behaviour steel fiber
    
    # ecr = o3.recorder.ElementToArrayCache(osi, ele)
    
    #Consider changing ele to the element number 
    # ecr_file = o3.recorder.ElementToArrayCache(osi, ele, arg_vals=['section', 1, 'force']) #Recording stress-Strain behaviour 
    # ecr = o3.recorder.ElementToArrayCache(osi, ele)
    # recorder('Element', '-file', "results/pushover/1brace.out",'-time','-ele', 11422,11423,700001,11425,11426,"section", 2,'fiber',str(0.05) , str(0), 'stressStrain')
    
    # comp_fibre = o3.recorder.ElementsToArrayCache(osi, ele, arg_vals=['section', 1, 'stressStrain'])
    
    # ndr = o3.recorder.NodesToArrayCache(osi, fd_nds, dofs=[o3.cc.Y], res_type='disp') # Node displacement recorder
    # efr = o3.recorder.ElementsToArrayCache(osi, ele, arg_vals=['section', 1, 'stressStrain']) #element force recorder - records moments
    # fr = o3.recorder.ElementsToArrayCache(osi, fd_eles, arg_vals=['section', 1, 'force']) #element force recorder - records moments
    # soil_ele_forces_recorder = o3.recorder.ElementsToArrayCache(osi, sl_eles, arg_vals=['localForce']) #records the soil spring force
    

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
    # ele_curvature = ecr.collect()[:, 0]  # Not sure if this is actually curvature.
    ele_SS_comp = ecr_comp.collect()
    ele_SS_tens = ecr_tens.collect()
    ele_SS_steel = ecr_steel.collect()
    # ecr_file = ecr_file.collect()

    moment = -nmr.collect()
    # comp_stress_strain = comp_fibre.collect()
    print("element_SS:", ele_SS_comp,ele_SS_tens,ele_SS_steel)
    #Estimate yield curvature
    
    rho_long = Ast/(b*d_eff) #Longtitudinal reinforcement ratio
    # print(rho_long)
    n = Es/Ec #Modular ratio
    k = np.sqrt((rho_long**2*n**2) + (2*rho_long*n)) - (rho_long*n) #Ratio of depth to N.A., Pujol et. al.
    # print(k)
    
    
    esp_yeild = fy/Es #Yeild strain
    curv_yeild = esp_yeild/((1-k)*d_eff) #Estimated yeild curvature
    print("Estimated yeild curvature:", curv_yeild)
    
    
    return moment, rot, mom_crack, mom_yield, ele_SS_comp, ele_SS_tens, ele_SS_steel


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    mom, curve ,mom_crack, mom_yield, ele_SS_comp, ele_SS_tens, ele_SS_steel = get_moment_curvature(num_incr=1000, max_curve=0.2)

    # m_cr = mom_crack
    # m_y = mom_yield
    # m_ult = mom[-1]/1000
    
    # print("Moments:", m_cr, m_y, m_ult)
    
    # plt.axhline(m_cr, c="green", label="Cracking Moment (Hand Calc)")
    # plt.axhline(m_y, c="orange", label="Yield Moment (Hand Calc)")
    # plt.axhline(m_ult, c="red", label="Ultimate Moment")

    plt.plot(curve, mom/1000)
    # plt.plot(ele_SS_comp, c="r")
    # plt.imshow(ele_SS_comp)
    # plt.plot(ele_SS_tens, c="b")
    # plt.plot(ele_SS_steel, c="k")
    
    # plt.plot(ele_SS)
    # plt.plot(comp_stress_strain)
    # print(curve, mom)
    
    curve_list = curve.tolist()
    mom_list = mom.tolist()
    # print(curve_list, mom_list)
    
    np.savetxt('m_curv_opensees_vTEST_1930_59_Hend.csv', (mom_list, curve_list), fmt='%.4e', delimiter=',')# X is an array

    fpc = 17.2e6  # Conc compressive strength (Pa)
    Ec = 5000 * (fpc ** 0.5) * 1E3
    d = 0.6  # section depth (m)
    b = 0.15  # section width (m)
    ei_el = Ec * b * d ** 3 / 12/1000
    moms = np.array([0, 60])
    curvs = moms / ei_el
    # plt.plot(curvs, moms, c='k')
    print(curvs, moms)

    plt.grid()
    
    axes = plt.gca()
    axes.set_xlabel("Curvature [1/m]")
    axes.set_ylabel("Bending Moment [kNm]")
    # axes.set_xlim([0, 0])
    # axes.set_xticks([0.0001])
    axes.legend()
    
    

    plt.show()
