import o3seespy as o3
import numpy as np


def run(osi, lf, dx, xd, yd, ksoil, udl, num_incr, sm_sect, o3_sect):
    # spring_spacing = 0.5 # [m] spring spacing

    s_width = sm_sect.width # section width [m]
    s_depth = sm_sect.depth  # section depth [m]

    nnodes = max(int(lf / dx) + 1, 2)  # number of nodes
    nx = np.linspace(0, lf, nnodes)  # x-position of nodes
    x_trib = np.ones_like(nx)  # Trib area
    x_trib[1:] += np.diff(nx)
    x_trib[:-1] += np.diff(nx)

    ks = ksoil * s_width * x_trib  # soil spring stiffness N/m

    print('hi ks', ks)

    # Could not find this in henderson thesis - does not expicitly mention it
    # Surely the force resistance should be related to the stiffness of each spring then multiply this by the displacement it is subject to

    # ks = ksoil * x_trib *s_width # N/m (TODO: define this better, see Henderson thesis for details) spring stiffness
    # Spring stiffness should be defined as ks = ksoil * section width * spring spacing
    imposed_displacement = 0.1  # [m]
    py = ks * imposed_displacement  # N Capacity (TODO: define this better, see Henderson thesis for details) spring force resistance
    # thesis does not define this... I could not find it

    transf = o3.geom_transf.Linear2D(osi, [])  # Using simple transform - no P-delta effects

    fd_nds = []  # foundation nodes
    sl_nds = []  # soil nodes
    sl_eles = []  # soil elements
    fd_eles = []  # foundation beam elements

    "Soil Nodes"
    for i in range(len(nx)):
        fd_nds.append(o3.node.Node(osi, nx[i], 0.0))
        sl_nds.append(o3.node.Node(osi, nx[i], 0.0))
        o3.Mass(osi, fd_nds[-1], 1.0, 1.0, 1.0)

        dettach = 1
        mat_base = o3.uniaxial_material.ElasticPP(osi, 1* ks[i], 1 * py[i] / ks[i],
                                                  -1 * py[i] / ks[i])  # low tension stiffness
        if dettach:
            mat_obj2 = o3.uniaxial_material.Elastic(osi, 10000 * ks[i], eneg=0.000001 * ks[i])
            mat = o3.uniaxial_material.Series(osi, [mat_base, mat_obj2])
        else:
            mat = o3.uniaxial_material.ElasticPP(osi, ks[i], 1 * py[i] / ks[i],
                                                 -1 * py[i] / ks[i])  # remove tension capacity
        sl_eles.append(o3.element.ZeroLength(osi, [fd_nds[-1], sl_nds[-1]], mats=[mat], dirs=[o3.cc.Y]))
        "Foundation nodes"
        if i != 0:
            fc = 17.2e6  # [Pa]
            e_mod = 5000 * np.sqrt(fc) * 1e3  # Pa

            area = s_width * s_depth
            i_z = s_depth ** 3 * s_width / 12

            integ = o3.beam_integration.Lobatto(osi, o3_sect, big_n=5)
            fd_eles.append(o3.element.DispBeamColumn(osi, [fd_nds[i], fd_nds[i - 1]], transf=transf, integration=integ))

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

    yd = np.array(yd) - udl / ksoil
    xd0 = xd[0]  # start position of displaced section
    xd1 = xd[1]  # end position of displaced section
    yd0 = yd[0]  # target displacement at start
    yd1 = yd[1]  # target displacement at end
    ind0 = np.argmin(abs(nx - xd0))
    ind1 = np.argmin(abs(nx - xd1))
    xd0n = nx[ind0]  # actual position of node where displacement will be imposed
    xd1n = nx[ind1]
    for i in range(ind0, ind1 + 1):
        xi = nx[i]
        ydi = yd0 + (yd1 - yd0) / (xd1n - xd0n) * (xi - xd0n)  # linear interpolate
        o3.SP(osi, sl_nds[i], o3.cc.Y, [ydi])  # impose displacement

    ymin = min([yd0, yd1])
    ydiff = ymin - -udl / ksoil
    max_ind = [ind0, ind1][np.argmin([yd0, yd1])]  # use the node that has the largest displacement

    # o3.integrator.DisplacementControl(osi, fd_nds[max_ind], dof=o3.cc.Y, incr=fd_incs, num_iter=10)
    o3.integrator.LoadControl(osi, incr=1 / num_incr)

    ndr = o3.recorder.NodesToArrayCache(osi, fd_nds, dofs=[o3.cc.Y], res_type='disp')  # Node displacement recorder
    efr = o3.recorder.ElementsToArrayCache(osi, fd_eles,
                                           arg_vals=['section', 1, 'force'])  # element force recorder - records moments
    soil_ele_forces_recorder = o3.recorder.ElementsToArrayCache(osi, sl_eles, arg_vals=[
        'localForce'])  # records the soil spring force

    ndisps = [[]]
    mom = [[]]
    shear = [[]]

    for j in range(nnodes):
        ndisps[0].append(o3.get_node_disp(osi, fd_nds[j], dof=o3.cc.Y))

    for j in range(nnodes - 1):
        mom[0].append(o3.get_ele_response(osi, fd_eles[j], 'force')[5])

    for j in range(nnodes - 1):
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

        for j in range(nnodes - 1):
            # print('force: ', o3.get_ele_response(osi, fd_eles[j], 'force'))
            mom[-1].append(o3.get_ele_response(osi, fd_eles[j], 'force')[5])

        for j in range(nnodes - 1):
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
    print('ndisps', ndisps)
    print("moment", mom)
    print("shear", shear)
    print("SPRINGZZ", spring_forces[-1], len(spring_forces))

    return Output(nx, node_disps, forces, [xd0n, xd1n], [yd0, yd1], mom, shear, ndisps, spring_forces)


class Output(object):
    def __init__(self, nx, node_disps, forces, xdn, yd, mom, shear, ndisps, spring_forces):
        self.nx = nx
        self.node_disps = node_disps
        self.forces = forces
        self.xdn = xdn
        self.yd = yd
        self.mom = mom
        self.shear = shear
        self.ndisps = ndisps
        self.spring_forces = spring_forces
