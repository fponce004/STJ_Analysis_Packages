import numpy as np

def fitPul00(t, to, A, td, tr):
    #to => offset
    #A  => Amplitude positive or negative
    #td => decay time
    #tr => rise time
    t = np.array(t)
    td = abs(td)
    tr = abs(tr)
    i = 0
    if to > t.max():
        i = len(t)
    else:
        i = np.where(to < t)[0][0]
    pre = np.zeros(i)
    post = (A*(1 - np.exp(-(t[i:]-to)/tr))*np.exp(-(t[i:]-to)/td))
    return np.insert(pre, pre.size, post)


def fitPul01(t, b, m, to_1, A_1, td_1, tr_1):
    return b + m*t + fitPul00(t, to_1, A_1, td_1, tr_1)

def fitPul02(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul03(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2, to_3, A_3, td_3, tr_3]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul04(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2, to_3, A_3, td_3, tr_3,
           to_4, A_4, td_4, tr_4]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul05(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
             to_5, A_5, td_5, tr_5):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
           to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
           to_5, A_5, td_5, tr_5]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul06(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
             to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
           to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
           to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul07(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
             to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
             to_7, A_7, td_7, tr_7):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
           to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
           to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
           to_7, A_7, td_7, tr_7]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul08(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
             to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
             to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
           to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
           to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
           to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul09(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
             to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
             to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8,
             to_9, A_9, td_9, tr_9):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
           to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
           to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
           to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8,
           to_9, A_9, td_9, tr_9]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul10(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
             to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
             to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8,
             to_9, A_9, td_9, tr_9, to_10, A_10, td_10, tr_10):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
           to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
           to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
           to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8,
           to_9, A_9, td_9, tr_9, to_10, A_10, td_10, tr_10]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul11(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
             to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
             to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8,
             to_9, A_9, td_9, tr_9, to_10, A_10, td_10, tr_10,
             to_11, A_11, td_11, tr_11):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
           to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
           to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
           to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8,
           to_9, A_9, td_9, tr_9, to_10, A_10, td_10, tr_10,
           to_11, A_11, td_11, tr_11]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])

def fitPul12(t, b, m, to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
             to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
             to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
             to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8,
             to_9, A_9, td_9, tr_9, to_10, A_10, td_10, tr_10,
             to_11, A_11, td_11, tr_11, to_12, A_12, td_12, tr_12):
    par = [to_1, A_1, td_1, tr_1, to_2, A_2, td_2, tr_2,
           to_3, A_3, td_3, tr_3, to_4, A_4, td_4, tr_4,
           to_5, A_5, td_5, tr_5, to_6, A_6, td_6, tr_6,
           to_7, A_7, td_7, tr_7, to_8, A_8, td_8, tr_8,
           to_9, A_9, td_9, tr_9, to_10, A_10, td_10, tr_10,
           to_11, A_11, td_11, tr_11, to_12, A_12, td_12, tr_12]
    return b + m*t + sum([fitPul00(t, *par[4*p: 4*(p+1)]) for p in range(int(len(par)/4))])
