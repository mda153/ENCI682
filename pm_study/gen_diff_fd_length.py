import os
import o3seespy as o3
import numpy as np
import matplotlib.pyplot as plt
import sfsimodels as sm
import loc_o3fibre
import loc_fd_analysis
from bwplot import cbox
import all_paths as ap


def create(save_results=0, save_fig=0, show=0):  # creates plot
    name = os.path.splitext(os.path.basename(os.path.realpath(__file__)))[0]
    name = name.replace('gen_', '')
    version = 'v1'
    out_folder = f'{ap.PM_STUDY_RESULTS}{name}/{version}/'
    try:
        os.makedirs(out_folder)
    except FileExistsError:
        pass

    bf, ax = plt.subplots(nrows=3, squeeze=True)
    # grabs relevant values from run function
    lfs = [10, 20, 30]

    m_crack = -10900
    m_yield = -32000
    m_ult = -48000

    support_disp = -0.5

    len_loss_sup = 2.0

    # stypes = [stypes[0]]
    for ss, lf in enumerate(lfs):
        prefix = f'LF{lf}'.replace('.', 'p')
        depth = 0.6  # [m] section depth
        width = 0.13  # [m] section width

        fpc = 17.2e6  # Conc compressive strength (Pa)
        ft = 1.40e6  # [Pa] concrete tensile strength
        eps_conf_conc_max = 0.002
        eps_conc_crush_unconf = 0.004  # Concrte strain at crushing unconf.
        eps_conc_crush_conf = 0.03  # conf conc strain at crushing conf.
        Ec = 5000 * np.sqrt(fpc) * 1e3  # Concrete modulus of elasticity (Pa) - eqn in Mander script for auto calc

        fy = 300.0e6  # Pa
        Es = 200.0e9  # Steel modulus [Pa]

        espall = 0.0064  # Maximum uncon. conc. strain (usually 0.0064)

        sm_sect = sm.sections.RCDetailedSection(depth=depth, width=width)
        sm_sect.layer_depths = [0.044, 0.056, 0.556]  # 44mm cover (m units)
        sm_sect.bar_diams = [[0.012], [0.012], [0.012, 0.012]]  # 12mm bars - two 12mm bars directly beside each other at top of section - model this as 1 x 24mm bar
        sm_sect.bar_centres = [[width / 2], [width / 2], [0.044, 0.086]]

        osi = o3.OpenSeesInstance(ndm=2, ndf=3)  # 2D with 3 DOF

        conc_conf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-eps_conf_conc_max, epscu=-eps_conc_crush_conf,
                                                    ec=Ec, fct=ft, et=0)  # Confined concrete paramters
        conc_unconf = o3.uniaxial_material.Concrete04(osi, fc=-fpc, epsc=-espall, epscu=-eps_conc_crush_unconf,
                                                      ec=Ec, fct=ft, et=0)  # unconfined concrete properties
        rebar = o3.uniaxial_material.Steel01(osi, fy=fy, e0=Es, b=0.01)  # Reinforcing steel

        nf_core_y = 8  # number of fibers in Y-Y dir - for the concrete core - (10 - 20 is a good number here)
        nf_core_z = 8  # number of fibers in Z-Z dir
        nf_cover_y = 10  # number of fibers Y-Y dir, including cover concrete
        nf_cover_z = 10  # number of fibers Z-Z dir, including cover concrete

        fd_sect = loc_o3fibre.get_rc_fibre_section(osi, sm_sect, conc_conf, conc_unconf, rebar, nf_core_y, nf_core_z,
                                                   nf_cover_y, nf_cover_z)
        xd = ((lf - len_loss_sup) / 2, (lf + len_loss_sup) / 2)
        yd = (support_disp, support_disp)
        oo = loc_fd_analysis.run(osi, lf=lf, dx=0.5, xd=xd, yd=yd, ksoil=3.847e6, udl=8600, num_incr=500, sm_sect=sm_sect, o3_sect=fd_sect)

        ax[0].plot(oo.nx, oo.ndisps[0], label=f'Initial {lf}m', c=cbox(ss), ls='-.')
        # ax[0].plot(nx, ndisps[int(len(ndisps) / 2)], label='Linear @ 50%', c="b", ls = "--") #Plotting for 50% but this is no longer 50%!! Depends on array length!!
        ax[0].plot(oo.nx, oo.ndisps[-1], label=f'Final {lf}m', c=cbox(ss))
        el_x = (oo.nx[1:] + oo.nx[:-1] ) / 2

        ax[1].plot(el_x, oo.mom[-1], label=f'{lf}m', c=cbox(ss))  # plots

        ax[2].plot(oo.nx, oo.spring_forces[-1], label=f'{lf}m', c=cbox(ss))

        ax[0].axvspan(*xd, color=(0.3, 0.3, 0.3, 0.5), zorder=0)
        if save_results:
            outs = ['nx', 'node_disps', 'forces', 'mom', 'shear', 'ndisps', 'spring_forces']
            for item in outs:
            # for item in oo.__dict__:
                print(item)
                np.savetxt(out_folder + f'{prefix}_{item}.txt', getattr(oo, item))

    # ax[0].plot([xd0n_elastic, xd1n_elastic], [yd0_elastic, yd1_elastic], c='r', label='Imposed') #plots imposed displacement
    ax[1].axhline(m_crack, ls='-.', label="Cracking Moment", c="g")
    ax[1].axhline(m_yield, ls="-.", label='Yield Moment', c="orange")
    ax[1].axhline(m_ult, ls="-.", label="Ultimate Moment", c="k")

    ax[0].set_xlabel('Foundation Length (m)')
    ax[0].set_ylabel('Vertical Disp. (m)')
    ax[1].set_ylabel("Bending Moment (Nm)")
    ax[1].set_xlabel("Foundation Length (m)")
    # ax[2].set_ylabel("Force (N)")
    # plt.grid()

    ax[1].legend(bbox_to_anchor=(1, 1))
    ax[0].legend(bbox_to_anchor=(1, 1))
    # ax[2].legend(bbox_to_anchor=(1,1))

    ax[0].grid()
    ax[1].grid()
    # ax[2].grid()

    if save_fig:

        ffp = out_folder + name + '.png'
        bf.savefig(ffp, dpi=150)
    if show:
        plt.show()


if __name__ == '__main__':

    create(save_results=1, save_fig=1, show=1)

