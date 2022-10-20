from cumbiapy import steel
import matplotlib.pyplot as plt


def create_raynor():
    fy = 1.1 * 414  # long steel yielding stress (MPa)
    fyh = 1.1 * 414  # transverse steel yielding stress (MPa)
    Es = 200000  # steel modulus of elasticity
    fsu = 1.4 * fy  # long steel max stress (MPa)*
    esh = 0.008  # long steel strain for strain hardening (usually 0.008)*
    esu = 0.10  # long. steel maximum strain (usually ~0.10-0.15)*

    Ey = 350  # slope of the yield plateau (MPa)
    C1 = 3.5  # defines strain hardening curve in the Raynor model [2-6]
    dels = 0.0001  # delta strain for default material models
    es, fs = steel.raynor(Es, fy, fsu, esh, esu, dels, C1, Ey)
    return es, fs


def create_king():
    fy = 1.1 * 414  # long steel yielding stress (MPa)
    fyh = 1.1 * 414  # transverse steel yielding stress (MPa)
    Es = 200000  # steel modulus of elasticity
    fsu = 1.4 * fy  # long steel max stress (MPa)*
    esh = 0.008  # long steel strain for strain hardening (usually 0.008)*
    esu = 0.10  # long. steel maximum strain (usually ~0.10-0.15)*

    Ey = 350  # slope of the yield plateau (MPa)
    C1 = 3.5  # defines strain hardening curve in the Raynor model [2-6]
    dels = 0.0001  # delta strain for default material models
    es, fs = steel.king(Es, fy, fsu, esh, esu, dels)
    return es, fs


def show_models():
    es, fs = create_raynor()
    plt.plot(es, fs, label='Raynor')
    es, fs = create_king()
    plt.plot(es, fs, label='King')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    show_models()



