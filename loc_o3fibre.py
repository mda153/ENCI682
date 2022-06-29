import numpy as np
import o3seespy as o3


def get_rc_fibre_section(osi, sm_sect, conc_conf_mat, conc_unconf_mat, rebar_mat, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z):
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
    if osi is None:
        osi = o3.OpenSeesInstance(ndm=2, ndf=3, state=3)

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
    o3.patch.Quad(osi, conc_conf_mat, nf_core_z, nf_core_y,  # core, counter-clockwise (diagonals at corners)
                  crds_i=[-core_y, core_z],
                  crds_j=[-core_y, -core_z],
                  crds_k=[core_y, -core_z],
                  crds_l=[core_y, core_z])

    o3.patch.Quad(osi, conc_unconf_mat, 1, nf_cover_y,  # right cover, counter-clockwise (diagonals at corners)
                  crds_i=[-edge_y, edge_z],
                  crds_j=[-core_y, core_z],
                  crds_k=[core_y, core_z],
                  crds_l=[edge_y, edge_z])
    o3.patch.Quad(osi, conc_unconf_mat, 1, nf_cover_y,  # left cover
                  crds_i=[-core_y, -core_z],
                  crds_j=[-edge_y, -edge_z],
                  crds_k=[edge_y, -edge_z],
                  crds_l=[core_y, -core_z])
    o3.patch.Quad(osi, conc_unconf_mat, nf_cover_z, 1,  # bottom cover
                  crds_i=[-edge_y, edge_z],
                  crds_j=[-edge_y, -edge_z],
                  crds_k=[-core_y, -core_z],
                  crds_l=[-core_y, core_z])
    o3.patch.Quad(osi, conc_unconf_mat, nf_cover_z, 1,  # top cover
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
                area = sm_sect.bar_diams[k][i] ** 2 / 4 *_pi
                o3.section.gen_fibre_section(osi, y_pos, z_loc=sm_sect.bar_centres[k][i], area=area, mat=rebar_mat)
    else:
        import numpy as np
        for k in range(len(sm_sect.layer_depths)):
            n_bars = len(sm_sect.bar_centres[k])
            bar_area = np.mean(sm_sect.bar_diams[k]) ** 2 / 4 * _pi
            y_pos = sm_sect.layer_depths[k] - sm_sect.depth / 2
            o3.layer.Straight(osi, rebar_mat, n_bars, bar_area, start=[y_pos, core_z], end=[y_pos, -core_z])
    return sect


def get_mom_curve(sm_sect, max_curve, axial_load, num_incr=100):
    osi = o3.OpenSeesInstance(ndm=2, ndf=3, state=3)
    n1 = o3.node.Node(osi, 0.0, 0.0)
    n2 = o3.node.Node(osi, 0.0, 0.0)
    o3.Fix3DOF(osi, n1, 1, 1, 1)
    o3.Fix3DOF(osi, n2, 0, 1, 0)
    # TODO: update fibres to Mander mdoel Concrete07 or ConfinedConcrete
    conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-sm_sect.rc_mat.fc, epsc0=-0.005, fpcu=-sm_sect.rc_mat.fc*0.85, eps_u=-0.02)
    conc_unconf = o3.uniaxial_material.Concrete01(osi, fpc=-sm_sect.rc_mat.fc * 0.8, epsc0=-0.004, fpcu=0.0, eps_u=-0.006)
    rebar = o3.uniaxial_material.Steel01(osi, fy=sm_sect.rc_mat.fy, e0=sm_sect.rc_mat.e_mod_steel, b=0.02)
    sect = get_rc_fibre_section(osi, sm_sect, conc_conf, conc_unconf, rebar, 5, 5, 5, 5)
    ele = o3.element.ZeroLengthSection(osi, [n1, n2], sect)

    nd = o3.recorder.NodeToArrayCache(osi, n2, dofs=[o3.cc.DOF2D_ROTZ], res_type='disp')
    nm = o3.recorder.NodeToArrayCache(osi, n1, dofs=[o3.cc.DOF2D_ROTZ], res_type='reaction')

    ts = o3.time_series.Constant(osi)
    o3.pattern.Plain(osi, ts)
    o3.Load(osi, n2, load_values=[-axial_load, 0.0, 0.0])

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
    o3.extensions.to_py_file(osi)
    o3.wipe(osi)
    curvature = nd.collect()
    moment = -nm.collect()
    return moment, curvature


def get_mom_curve_w_vert_disp(sm_sect, max_curve, axial_load, num_incr=100):
    osi = o3.OpenSeesInstance(ndm=2, ndf=3, state=3)
    n1 = o3.node.Node(osi, 0.0, 0.0)
    n2 = o3.node.Node(osi, 0.0, 0.0)
    o3.Fix3DOF(osi, n1, 1, 1, 1)
    o3.Fix3DOF(osi, n2, 0, 1, 0)
    # TODO: update fibres to Mander mdoel Concrete07 or ConfinedConcrete
    conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-sm_sect.rc_mat.fc, epsc0=-0.005, fpcu=-sm_sect.rc_mat.fc*0.85, eps_u=-0.02)
    conc_unconf = o3.uniaxial_material.Concrete01(osi, fpc=-sm_sect.rc_mat.fc * 0.8, epsc0=-0.004, fpcu=0.0, eps_u=-0.006)
    rebar = o3.uniaxial_material.Steel01(osi, fy=sm_sect.rc_mat.fy, e0=sm_sect.rc_mat.e_mod_steel, b=0.02)
    sect = get_rc_fibre_section(osi, sm_sect, conc_conf, conc_unconf, rebar, 5, 5, 5, 5)
    ele = o3.element.ZeroLengthSection(osi, [n1, n2], sect)

    nd = o3.recorder.NodeToArrayCache(osi, n2, dofs=[o3.cc.DOF2D_ROTZ], res_type='disp')
    nd_vert = o3.recorder.NodeToArrayCache(osi, n2, dofs=[o3.cc.X], res_type='disp')
    nm = o3.recorder.NodeToArrayCache(osi, n1, dofs=[o3.cc.DOF2D_ROTZ], res_type='reaction')

    ts = o3.time_series.Constant(osi)
    o3.pattern.Plain(osi, ts)
    o3.Load(osi, n2, load_values=[-axial_load, 0.0, 0.0])

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
    o3.extensions.to_py_file(osi)
    o3.wipe(osi)
    curvature = nd.collect()
    vert_disp = nd_vert.collect()
    moment = -nm.collect()
    return moment, curvature, vert_disp


def run_example():

    import sfsimodels as sm
    col_sect = sm.sections.RCDetailedSection(depth=0.45, width=0.45)
    col_sect.layer_depths = [0.06, 0.2, 0.34]
    col_sect.bar_diams = [[0.025, 0.025, 0.025], [0.025, 0.025], [0.025, 0.025, 0.025]]
    col_sect.bar_centres = [[0.06, 0.2, 0.34], [0.06, 0.34], [0.06, 0.2, 0.34]]
    rc_mat = sm.materials.ReinforcedConcreteMaterial(fc=30e6, fy=300e6,
                                                     e_mod_steel=200e9, poissons_ratio=0.18)
    area_steel = 0.025 ** 2 / 4 * np.pi * 8
    a_sect = col_sect.depth * col_sect.width
    rho_steel = area_steel / a_sect
    print('rho_steel: ', rho_steel)
    col_sect.rc_mat = rc_mat
    import matplotlib.pyplot as plt
    nloads = [-100e3, 0, 200e3, 400e3, 700e3, 2400e3]
    bf, ax = plt.subplots(nrows=2)
    for nload in nloads:
        print('p rat: ', nload / a_sect / 1e6)
        mom, curve, vdisp = get_mom_curve_w_vert_disp(col_sect, axial_load=nload, max_curve=0.04)
        ax[0].plot(curve, mom, label=f'N={nload/1e3:.0f}kN ({nload/(col_sect.rc_mat.fc * a_sect):.2f})')
        ax[1].plot(curve, vdisp, label=f'N={nload / 1e3:.0f}kN ({nload / (col_sect.rc_mat.fc * a_sect):.2f})')

    ax[0].legend()
    plt.show()






if __name__ == '__main__':
    run_example()