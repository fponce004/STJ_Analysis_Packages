import numpy as np
from scipy.special import wofz

def funVoigt_p(p, x):
    return sum([funVoigt01(x, *np.abs(p[4*val:4*(val+1)])) for val in range(int(len(p)/4))])

def funVoigt01(x, a1, g1_1, g1_2, c1):
    return a1*np.real(wofz(((x - c1) + 1j*g1_2)/(g1_1*np.sqrt(2))))/(g1_1*np.sqrt(2*np.pi))

def funVoigt02(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt03(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt04(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt05(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3, a4,
               g4_1, g4_2, c4, a5, g5_1, g5_2, c5):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt06(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt07(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt08(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt09(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt10(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt11(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt12(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt13(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt14(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt15(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt16(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt17(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt18(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt19(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
               a19, g19_1, g19_2, c19):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
                           a19, g19_1, g19_2, c19]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt20(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
               a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
                           a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt21(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
               a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20, a21, g21_1, g21_2, c21):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
                           a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20, a21, g21_1, g21_2, c21]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt22(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
               a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20, a21, g21_1, g21_2, c21,
               a22, g22_1, g22_2, c22):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
                           a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20, a21, g21_1, g21_2, c21,
                           a22, g22_1, g22_2, c22]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt23(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
               a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20, a21, g21_1, g21_2, c21,
               a22, g22_1, g22_2, c22, a23, g23_1, g23_2, c23):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
                           a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20, a21, g21_1, g21_2, c21,
               a22, g22_1, g22_2, c22, a23, g23_1, g23_2, c23]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])

def funVoigt24(x, a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
               a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
               a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
               a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
               a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
               a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
               a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20, a21, g21_1, g21_2, c21,
               a22, g22_1, g22_2, c22, a23, g23_1, g23_2, c23, a24, g24_1, g24_2, c24):
    par = np.abs(np.array([a1, g1_1, g1_2, c1, a2, g2_1, g2_2, c2, a3, g3_1, g3_2, c3,
                           a4, g4_1, g4_2, c4, a5, g5_1, g5_2, c5, a6, g6_1, g6_2, c6,
                           a7, g7_1, g7_2, c7, a8, g8_1, g8_2, c8, a9, g9_1, g9_2, c9,
                           a10, g10_1, g10_2, c10, a11, g11_1, g11_2, c11, a12, g12_1, g12_2, c12,
                           a13, g13_1, g13_2, c13, a14, g14_1, g14_2, c14, a15, g15_1, g15_2, c15,
                           a16, g16_1, g16_2, c16, a17, g17_1, g17_2, c17, a18, g18_1, g18_2, c18,
                           a19, g19_1, g19_2, c19, a20, g20_1, g20_2, c20, a21, g21_1, g21_2, c21,
               a22, g22_1, g22_2, c22, a23, g23_1, g23_2, c23, a24, g24_1, g24_2, c24]))
    return sum([funVoigt01(x, *par[4*val:4*(val+1)]) for val in range(int(len(par)/4))])
