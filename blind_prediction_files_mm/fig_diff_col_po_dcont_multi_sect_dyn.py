import eqdes.concrete
import numpy as np
import o3seespy as o3
from loc_o3fibre import get_rc_fibre_section
from bwplot import cbox
import sfsimodels as sm
import pandas as pd
import all_paths as ap


def pushover(hinge, disps):
    osi = o3.OpenSeesInstance(ndm=2, state=3)

    clear_height = 2.667  # m  Line 3
    sdepth = 0.356  # m (section depth)
    enlarged_end_sections = 0.508  # m
    bar_area = 819e-6  # m2  (Para 2)
    trans_bar_area = 71.0e-6  # m2 (Para 2)
    trans_spacing_ends = 0.076  # m
    ends_length = 0.305  # m
    trans_spacing_mid = 0.326  # m
    f_c28 = 17.71e6  # Pa
    fy_steel_long = 478.0e6  # Pa
    fy_steel_trans = 379.0e6  # Pa
    fu_steel_trans = 568.0e6  # Pa
    fu_steel_long = 680.0e6  # Pa
    e_mod_steel = 200.0e9  # Pa (From FAQs E=29000ksi)
    r_strain_hardening = 5  # FAQ: hardening strain of 4 to 6 times yield strain at the strain-hardening
    conc_unit_weight = 24.0e3
    top_weight = 1.9812 * enlarged_end_sections ** 2 * conc_unit_weight

    # TODO: see Eq 2-10 strength increase due to foundation connection (CL 2.6.5.6)
    a_gross = sdepth ** 2
    phi_o_fy = 1.25  # CL 2.6.5.6
    nload_base = a_gross * clear_height * conc_unit_weight + top_weight
    f_os_conn = phi_o_fy + 2 * (nload_base / f_c28 / a_gross - 0.1) ** 2
    # Overstrength from base connection is 1.6%  - ignore effect
    phi_y = 2 * fy_steel_long / (e_mod_steel * sdepth)  # Eq C2-2/2 in CL 2.6.1.3.4 # is h the full section depth -pg 25

    d_bar = np.sqrt(4 * bar_area / np.pi)
    d_bar_trans = np.sqrt(4 * trans_bar_area / np.pi)
    INCH_TO_METRE = 0.0254
    bar_centres = 1.5 * INCH_TO_METRE + d_bar_trans
    a_cv = sdepth * (sdepth - bar_centres)
    v_max = min([0.2 * f_c28, 8e6])
    print('v_max: ', v_max)
    v_base_max = a_cv * v_max
    print('v_base_max: ', v_base_max)
    fc = f_c28

    # define a reinforced concrete section using sfsimodels package
    spacings = [trans_spacing_ends, trans_spacing_mid, trans_spacing_mid, trans_spacing_ends]
    sects = []
    for i in range(4):
        col_sect = sm.sections.RCDetailedSection()
        col_sect.depth = sdepth  # m Line 2
        col_sect.width = sdepth  # m Line 2
        col_sect.layer_depths = [bar_centres, sdepth / 2, sdepth - bar_centres]
        col_sect.bar_diams = [[d_bar, d_bar, d_bar], [d_bar, d_bar], [d_bar, d_bar, d_bar]]
        col_sect.bar_centres = [[bar_centres, sdepth / 2, sdepth - bar_centres],
                                [bar_centres, sdepth - bar_centres],
                                [bar_centres, sdepth / 2, sdepth - bar_centres]]
        rc_mat = sm.materials.ReinforcedConcreteMaterial(fc=f_c28, fy=fy_steel_long,
                                                         e_mod_steel=e_mod_steel, poissons_ratio=0.18)
        col_sect.db_trans = d_bar_trans
        col_sect.fy_trans = fy_steel_trans
        col_sect.nb_trans_x = 2
        col_sect.nb_trans_y = 2
        col_sect.spacing_trans = spacings[i]
        rc_mat.unit_weight = 24.0e3
        col_sect.rc_mat = rc_mat
        f_crack = 0.4
        fc_conf, eps_conf = eqdes.concrete.calc_confined_via_mander_1988(col_sect, fc_unconf=f_c28)
        col_sect.fc_conf = fc_conf
        col_sect.eps_conf = eps_conf
        sects.append(col_sect)
    e_mod_conc = sm.materials.calc_e_mod_conc_via_park_and_paulay_1975(f_c28)
    f_conc_tens = 0.36 * np.sqrt(fc)  # NZS3101 CL 5.2.6
    eps_conc_tens = f_conc_tens / e_mod_conc  # Park and Paulay 1975 assume linear
    eps_unconf_conc_max = 0.002  # From Pauley Fig 2.1
    eps_conc_crush = 0.03  # Koopaee (2015) from NZS 31010?? taken from Fig 2.1 fc=20MPa curve
    eps_ult = 0.15  # From https://cdn.ymaws.com/concretenz.org.nz/resource/resmgr/docs/conf/2013/s5_p4.pdf
    r = 5
    epy_y = fy_steel_long / e_mod_steel
    eps_sh = r * epy_y
    lsr = trans_spacing_ends / d_bar
    height = clear_height  # m
    t_1 = 0.7  # s
    xi = 0.03  # critical damping ratio

    # Establish nodes
    ypos = [0, ends_length, clear_height / 2, clear_height - ends_length, clear_height]
    dh = np.diff(ypos)
    node_masses = np.zeros_like(ypos)
    node_masses[:-1] = dh / 2 * a_gross * conc_unit_weight / 9.8
    node_masses[1:] = dh / 2 * a_gross * conc_unit_weight / 9.8
    node_masses[-1] += top_weight / 9.8
    nodes = []
    for i in range(len(ypos)):
        nodes.append(o3.node.Node(osi, 0, ypos[i]))

    # Fix bottom node
    o3.Fix3DOF(osi, nodes[0], o3.cc.FIXED, o3.cc.FIXED, o3.cc.FIXED)
    o3.Fix3DOF(osi, nodes[-1], o3.cc.FREE, o3.cc.FREE, o3.cc.FIXED)

    # nodal mass (weight / g):
    o3.Mass(osi, nodes[-1], top_weight / 9.8, 0., 0.)

    # Define material
    transf = o3.geom_transf.Linear2D(osi, [])  # can change for P-delta effects
    area = 1.0
    nf_core_y = 10
    nf_core_z = 10
    nf_cover_y = 10
    nf_cover_z = 10

    e_mod = sm.materials.calc_e_mod_conc_via_park_and_paulay_1975(f_c28)
    v_eles = []
    for i in range(len(spacings)):
        col_sect = sects[i]
        iz = col_sect.width * col_sect.depth ** 3 / 12 * f_crack
        col_sect = sects[i]
        if hinge == 'elastic':

            bot_sect = o3.section.Elastic2D(osi, e_mod, area, iz)
        elif hinge == 'concentrated-nl':
            m_cap = 3 * bar_area * fy_steel_long * sdepth * 0.7
            # m_cap = 150.0e3  # Nm
            e0 = e_mod * iz
            print('e0: ', e0)
            print('m_cap: ', m_cap)
            mat_flex_lhs = o3.uniaxial_material.Steel01(osi, m_cap, e0=e_mod * iz, b=0.005)
            mat_axial = o3.uniaxial_material.Elastic(osi, e_mod * area)
            mat_shear = o3.uniaxial_material.Elastic(osi, e_mod * area)
            bot_sect = o3.section.Aggregator(osi, mats=[[mat_axial, o3.cc.P], [mat_shear, o3.cc.V_Y], [mat_flex_lhs, o3.cc.M_Z]])
        elif hinge == 'rc-fibre01-og':

            conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-col_sect.rc_mat.fc, epsc0=-0.005, fpcu=-col_sect.rc_mat.fc * 0.85, eps_u=-0.02)
            conc_unconf = o3.uniaxial_material.Concrete01(osi, fpc=-col_sect.rc_mat.fc * 0.8, epsc0=-0.004, fpcu=0.0, eps_u=-0.006)
            rebar = o3.uniaxial_material.Steel01(osi, fy=col_sect.rc_mat.fy, e0=col_sect.rc_mat.e_mod_steel, b=0.02)

            bot_sect = get_rc_fibre_section(osi, col_sect, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z)
        elif hinge == 'rc-fibre01':

            conc_conf = o3.uniaxial_material.Concrete01(osi, fpc=-fc_conf, epsc0=-eps_conf, fpcu=-fc_conf * 0.85, eps_u=-0.02)
            conc_unconf = o3.uniaxial_material.Concrete01(osi, fpc=-col_sect.rc_mat.fc, epsc0=-0.004, fpcu=0.0, eps_u=-0.006)
            rebar = o3.uniaxial_material.Steel01(osi, fy=col_sect.rc_mat.fy, e0=col_sect.rc_mat.e_mod_steel, b=0.001)

            bot_sect = get_rc_fibre_section(osi, col_sect, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z)
        elif hinge == 'rc-fibre04':

            conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fc_conf, epsc=-eps_conf, epscu=-eps_conc_crush, ec=e_mod_conc, fct=f_conc_tens, et=eps_conc_tens)
            conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-f_c28, epsc=-eps_unconf_conc_max, epscu=-eps_conc_crush, ec=e_mod_conc, fct=f_conc_tens, et=eps_conc_tens)
            rebar = o3.uniaxial_material.Steel01(osi, fy=col_sect.rc_mat.fy, e0=col_sect.rc_mat.e_mod_steel, b=0.01)
            # rebar = o3.uniaxial_material.ReinforcingSteelDMBuck(osi, fy_steel_long, fu_steel_long, e_mod_steel, e_mod_sh=0.05*e_mod_steel, eps_sh=eps_sh, eps_ult=eps_ult, lsr=lsr)
            # rebar = o3.uniaxial_material.ReinforcingSteel(osi, fy_steel_long, fu_steel_long, e_mod_steel,
            #                                                     e_mod_sh=0.05 * e_mod_steel, eps_sh=eps_sh * 4,
            #                                                     eps_ult=eps_ult)

            bot_sect = get_rc_fibre_section(osi, col_sect, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z, nf_cover_y, nf_cover_z)
        else:
            raise ValueError(f'hinge type is not valid: {hinge}, must be "elastic", "concentrated-nl" or "rc-fibre04"')

        # ph_len = 0.3  # m
        # bot_sect = o3.section.Elastic2D(osi, e_mod, area, iz)
        # centre_sect = o3.section.Elastic2D(osi, e_mod, area, iz)
        # integ = o3.beam_integration.HingeMidpoint(osi, bot_sect, ph_len, top_sect, ph_len, centre_sect)
        integ = o3.beam_integration.Lobatto(osi, bot_sect, big_n=5)
        # vert_ele = o3.element.NonlinearBeamColumn(osi, nodes[i: i+2], 6, bot_sect, transf, max_iter=10)
        vert_ele = o3.element.DispBeamColumn(osi, nodes[i: i + 2], transf, integ)
        v_eles.append(vert_ele)

    # set damping based on first eigen mode
    angular_freq = o3.get_eigen(osi, solver='fullGenLapack', n=1)[0] ** 0.5
    response_period = 2 * np.pi / angular_freq
    print('response_period: ', response_period)
    beta_k = 2 * xi / angular_freq
    o3.rayleigh.Rayleigh(osi, alpha_m=0.0, beta_k=beta_k, beta_k_init=0.0, beta_k_comm=0.0)
    ts0 = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts0)
    for i in range(1, len(nodes)):
        o3.Load(osi, nodes[i], [0, -node_masses[i] * 9.8, 0])

    o3.constraints.Transformation(osi)
    o3.test_check.NormDispIncr(osi, tol=1.0e-6, max_iter=35, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.numberer.RCM(osi)
    o3.system.SparseGeneral(osi)
    n_steps_gravity = 10
    d_gravity = 1. / n_steps_gravity
    o3.integrator.LoadControl(osi, d_gravity, num_iter=10)
    # o3.rayleigh.Rayleigh(osi, a0, a1, 0.0, 0.0)
    o3.analysis.Static(osi)
    o3.analyze(osi, num_inc=n_steps_gravity)
    o3.load_constant(osi, time=0.0)
    o3.wipe_analysis(osi)
    print('init_disp [top]: ', o3.get_node_disp(osi, nodes[-1], o3.cc.DOF2D_Y))

    # Start cyclic loading

    o3.constraints.Transformation(osi)
    o3.test_check.NormDispIncr(osi, tol=1.0e-5, max_iter=35, p_flag=0)
    o3.algorithm.Newton(osi)
    o3.numberer.RCM(osi)
    o3.system.FullGeneral(osi)
    o3.integrator.Newmark(osi, gamma=0.5, beta=0.25)
    xi_init = 0.02
    if xi_init:
        f2 = 10.0  # Hz
        period = 1.0
        omega_1 = 2 * np.pi / period
        # beta_k = 2 * xi_init / omega_1
        omega_2 = 2 * np.pi * f2
        a0 = 2 * xi_init * omega_1 * omega_2 / (omega_1 + omega_2)
        a1 = 2 * xi_init / (omega_1 + omega_2)
        o3.rayleigh.Rayleigh(osi, alpha_m=a0, beta_k=0.0, beta_k_init=0.0, beta_k_comm=a1)

    # o3.rayleigh.Rayleigh(osi, a0, a1, 0.0, 0.0)
    o3.analysis.Transient(osi)

    o3.record(osi)

    disps = list(disps)
    d_per_step = 0.0001
    d_per_dt = 0.0001
    dt = d_per_step / d_per_dt
    diff_disps = np.diff(disps, prepend=0)
    time_incs = np.abs(diff_disps) / d_per_dt
    approx_n_steps = time_incs / dt
    time_incs = np.where(approx_n_steps < 8, 8 * dt, time_incs)
    approx_n_steps = time_incs / dt
    assert min(approx_n_steps) >= 8, approx_n_steps
    curr_time = o3.get_time(osi)
    times = list(np.cumsum(time_incs) + curr_time)
    disps.append(disps[-1])
    times.append(1e10)

    disps = list(disps)
    n_steps_p2 = int((times[-2] - curr_time) / dt) + 10
    verbose = 1
    if verbose:
        print('n_steps: ', n_steps_p2)
    times.insert(0, curr_time)
    disps.insert(0, 0.0)
    top_node = nodes[-1]
    dof = o3.cc.DOF2D_X

    if top_node:
        init_disp = o3.get_node_disp(osi, top_node, dof)
        disps = list(np.array(disps) + init_disp)
        ts0 = o3.time_series.Path(osi, time=times, values=disps, factor=1)
        pat0 = o3.pattern.Plain(osi, ts0)
        if verbose:
            print('times: ', times)
            print('disps: ', disps)
        o3.SP(osi, top_node, dof=dof, dof_values=[1.0])

    rot_hinge = []
    mom_hinge = []
    applied_load = []
    col_top_xdisp = []
    col_top_ydisp = []
    sfs_f = [[], []]
    sfs_d = [[], []]
    ts0 = o3.time_series.Linear(osi, factor=1)
    o3.pattern.Plain(osi, ts0)
    target_d_inc = 0.0001
    o3.Load(osi, nodes[-1], [1.0, 0, 0])
    top_ss_node = nodes[-1]
    fail = 0
    for i in range(n_steps_p2):

        fail = o3.analyze(osi, 1, dt=1)
        if fail:
            print('Model failed')
            break
        o3.gen_reactions(osi)
        # print(o3.get_ele_response(osi, v_eles[0], 'force')[2], o3.get_ele_response(osi, v_eles[-1], 'force')[2])
        print(o3.get_ele_response(osi, v_eles[0], 'chordRotation')[1], o3.get_ele_response(osi, v_eles[-1], 'chordRotation')[1])
        # o3.get_ele_response(osi, vert_ele, 'sectionDisplacements')
        rot_hinge.append(o3.get_ele_response(osi, v_eles[0], 'chordRotation')[1])  # note first val in strain
        mom_hinge.append(o3.get_ele_response(osi, v_eles[0], 'force')[2])
        col_top_xdisp.append(o3.get_node_disp(osi, top_ss_node, o3.cc.DOF2D_X))
        col_top_ydisp.append(o3.get_node_disp(osi, top_ss_node, o3.cc.DOF2D_Y))
        applied_load.append(-o3.get_ele_response(osi, v_eles[0], 'force')[0])

    o3.extensions.to_py_file(osi)
    for i in range(30):
        print(1, i, o3.get_ele_response(osi, vert_ele, 'stressStrain', extra_args=['section', '1', 'fiber', f'{i}']))
    for i in range(30):
        print(2, i, o3.get_ele_response(osi, vert_ele, 'stressStrain', extra_args=['section', '2', 'fiber', f'{i}']))
    return np.array(rot_hinge), np.array(mom_hinge), np.array(col_top_xdisp), np.array(applied_load), np.array(col_top_ydisp)


def run():

    import matplotlib.pyplot as plt
    import engformat as ef
    clear_height = 2.667  # m  Line 3
    bf, ax = plt.subplots(nrows=3, ncols=2)
    ax = ax.flatten()
    hinge_types = ['elastic', 'concentrated-nl', 'rc-fibre04']
    hinge_types = ['rc-fibre04', 'rc-fibre01']
    hinge_types = ['concentrated-nl', 'rc-fibre01', 'rc-fibre04']
    hinge_types = ['rc-fibre04', 'concentrated-nl']
    hinge_types = ['rc-fibre04']
    # hinge_types = ['rc-fibre01']
    df = pd.read_excel(ap.ROOT_DIR + 'docs/loading_protocol.xlsx', skiprows=[1])
    print(df[:3])
    disps = np.array(df['disp_mm']) / 1e3
    # disps = disps[:20]
    # disps = [0.15]
    for i, htype in enumerate(hinge_types):
        rot_hinge, mom_hinge, col_top_xdisp, applied_load, col_top_ydisp = pushover(htype, disps)
        v_base = 290.0e3
        v_slope = (290 - 360) * 1e3 / (0.172 - 0.053)
        mdisps = np.maximum.accumulate(abs(col_top_xdisp))
        v_base = 360e3 + v_slope * (mdisps - 0.053)
        vals = np.where(abs(applied_load) > abs(v_base))
        if len(vals[0]) == 0:
            ind = None
        else:
            ind = vals[0][0]
        v_max = v_base[ind]
        ind75 = np.where(abs(applied_load) > 0.75 * abs(v_max))[0][0]

        ax[0].plot(rot_hinge[:ind], mom_hinge[:ind] / 1e3, c=cbox(i), label=htype)
        ax[1].plot(rot_hinge[:ind], col_top_ydisp[:ind] * 1e3, c=cbox(i))
        ax[2].plot(col_top_xdisp[:ind], applied_load[:ind] / 1e3, c=cbox(i))
        ax[3].plot(col_top_xdisp[:ind], c=cbox(i))
        ax[5].plot(mom_hinge[:ind] / 1e3, c=cbox(i))
        ax[4].plot(applied_load[:ind] / 1e3, c=cbox(i))
        print('v_max: ', v_max)
        print('col_top_xdisp[ind75]: ', col_top_xdisp[ind75])
        print('drift75: ', col_top_xdisp[ind75] / clear_height)
        print('max-drift: ', max(abs(col_top_xdisp[:ind])) / clear_height)
        print('duct: ', max(abs(col_top_xdisp[:ind])) / abs(col_top_xdisp[ind75] / 0.75))

        # full
        ax[0].plot(rot_hinge, mom_hinge / 1e3, c=cbox(i), ls=':')
        ax[1].plot(rot_hinge, col_top_ydisp * 1e3, c=cbox(i), ls=':')
        ax[2].plot(col_top_xdisp, applied_load / 1e3, c=cbox(i), ls=':')
        ax[2].plot(col_top_xdisp, v_base / 1e3, c='r')
        ax[2].plot(col_top_xdisp, -v_base / 1e3, c='r')
        ax[3].plot(col_top_xdisp, c=cbox(i), ls=':')
        ax[5].plot(mom_hinge / 1e3, c=cbox(i), ls=':')
        ax[4].plot(applied_load / 1e3, c=cbox(i), ls=':')
        ax[4].plot(v_base / 1e3, c='r')
        ax[4].plot(-v_base / 1e3, c='r')

    # ax[3].plot(disps, c='k')
    for axplot in ax:
        ef.xy(axplot)
    ax[0].legend()
    ax[0].set_ylabel('Moment [kNm]')
    ax[1].set_ylabel('Vertical disp. [mm]')
    ax[2].set_ylabel('Base shear [kN]')
    ax[4].set_ylabel('Base shear [kN]')
    ax[3].set_ylabel('Col top disp [m]')
    ax[4].set_ylabel('Moment [kNm]')
    ax[0].set_xlabel('Base chord rot. [rad]')
    ax[3].set_xlabel('Steps')
    ax[1].set_xlabel('Base chord rot. [rad]')
    ax[2].set_xlabel('Top of column xdisp [m]')
    bf.tight_layout()
    plt.show()


if __name__ == '__main__':



    run()
