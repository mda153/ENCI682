# =============================================================================
# Start of User Input Data
# =============================================================================

# input data:

name = 'Outputs'  # identifies actual work, the output file will be name.xls

# section properties:
D = 1067 - 22 * 4  # section diameter (mm)
clb = 50 - 22  # cover to longitudinal bars (mm)

# member properties

L = 1651  # member clear length (mm)
bending = 'single'  # single or double
ductilitymode = 'uniaxial'  # biaxial or uniaxial

# reinforcement details:

nbl = 30  # number of longitudinal bars
Dbl = 32  # long. bar diameter (mm)
Dh = 16  # diameter of transverse reinf. (mm)
stype = 'spirals'  # 'spirals' or 'hoops'*
s = 100  # spacing of transverse steel (mm)*

# aplieed loads:

P = 1325  # axial load kN (-) tension (+)compression

# material models (input the 'name' of the file with the stress-strain relationship
# to use the default models: Mander model for confined or unconfined  concrete type 'mc' or 'mu'.
# For lightweight confined concrete type 'mclw'
# King model for the steel 'ks', Raynor model for steel 'ra':

confined = 'mc'
unconfined = 'mu'
rebar = 'ks'

# material properties

fpc = 1.3 * 28  # concrete compressive strength (MPa)
Ec = 0  # concrete modulus of elasticity (MPa) or
# input 0 for automatic calculation using
# 5000(fpc)^0.5
eco = 0.002  # unconfined strain (usually 0.002 for normal weight or 0.004 for lightweight)*
esm = 0.10  # max transv. steel strain (usually ~0.10-0.15)*
espall = 0.0064  # max uncon. conc. strain (usually 0.0064)

fy = 1.1 * 414  # long steel yielding stress (MPa)
fyh = 1.1 * 414  # transverse steel yielding stress (MPa)
Es = 200000  # steel modulus of elasticity
fsu = 1.4 * fy  # long steel max stress (MPa)*
esh = 0.008  # long steel strain for strain hardening (usually 0.008)*
esu = 0.10  # long. steel maximum strain (usually ~0.10-0.15)*

Ey = 350  # slope of the yield plateau (MPa)
C1 = 3.5  # defines strain hardening curve in the Raynor model [2-6]

# this information is used only if the default material models are selected

# strain limits for yield surface (interaction diagram);

csid = 0.004  # concrete
ssid = 0.015  # steel

# Deformation Limit States:

ecser = 0.004;
esser = 0.015  # concrete (ecser) and steel (esser) serviceability strain
ecdam = 0.018;
esdam = 0.060  # concrete (ecser) and steel (esser) damage control strain
# (to use the 2/3 of the ultimate concrete strain just tipe 'twth'
# temperature information (in case of freezing conditions)
temp = 40  # temperature of the specimen in celsius
kLsp = 0.022  # constant to calculate Lsp = kLsp*fy*Dbl
# (usually 0.022 at ambient temp. or 0.011 at -40C)

# analysis control parameters:
itermax = 1000  # max number of iterations
ncl = 40  # # of concrete layers
tolerance = 0.001  # x fpc x Ag
dels = 0.0001  # delta strain for default material models

# =============================================================================
# End of User Input Data
# =============================================================================

