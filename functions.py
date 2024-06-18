import numpy as np
import pandas as pd

def prime(x):
    return 1000 * np.log(x / 1000 + 1)


def unprime(x):
    return (np.exp(x / 1000) - 1) * 1000


def Dp17O(d17O, d18O):
    return (prime(d17O) - 0.528 * prime(d18O)) * 1000


def d17O(d18O, Dp17O):
    return unprime(Dp17O / 1000 + 0.528 * prime(d18O))


def calculate_theta(d18O_A, Dp17O_A, d18O_B, Dp17O_B):
    a18 = (d18O_B + 1000) / (d18O_A + 1000)
    a17 = (d17O(d18O_B, Dp17O_B) + 1000) / (d17O(d18O_A, Dp17O_A) + 1000)

    theta = round(np.log(a17) / np.log(a18), 4)

    return theta


def apply_theta(d18O_A, Dp17O_A, d18O_B=None, shift_d18O=None, theta=None):

    if d18O_B == None:
        d18O_B = d18O_A + shift_d18O
    d17O_A = d17O(d18O_A, Dp17O_A)

    a18 = (d18O_B + 1000) / (d18O_A + 1000)
    a17 = a18**theta

    d17O_B = a17 * (d17O_A + 1000) - 1000
    Dp17O_B = Dp17O(d17O_B, d18O_B)

    return Dp17O_B


def mix_d17O(d18O_A, d17O_A=None, D17O_A=None, d18O_B=None, d17O_B=None, D17O_B=None, step=100):
    ratio_B = np.arange(0, 1+1/step, 1/step)

    if d17O_A is None:
        d17O_A = unprime(D17O_A/1000 + 0.528 * prime(d18O_A))

    if d17O_B is None:
        d17O_B = unprime(D17O_B/1000 + 0.528 * prime(d18O_B))

    mix_d18O = ratio_B * float(d18O_B) + (1 - ratio_B) * float(d18O_A)
    mix_d17O = ratio_B * float(d17O_B) + (1 - ratio_B) * float(d17O_A)
    mix_D17O = Dp17O(mix_d17O, mix_d18O)
    xB = ratio_B * 100

    df = pd.DataFrame(
        {'mix_d17O': mix_d17O, 'mix_d18O': mix_d18O, 'mix_Dp17O': mix_D17O, 'xB': xB})
    return df


def B_from_a(a, A):
    return (A + 1000) / a - 1000


def A_from_a(a, B):
    return (B + 1000) * a - 1000


def epsilon(d18O_A, d18O_B):
    epsilon = ((d18O_A + 1000) / (d18O_B + 1000) - 1) * 1000
    return epsilon
