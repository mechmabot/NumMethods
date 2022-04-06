import numpy as np
from matplotlib import pyplot as plt
from  module import *
from sys import argv

def FUNCTION_K(x , t):
    return 1 / (((np.sin((x+t)/2)**2))+0.25*(np.cos((x+t)/2)**2))

def FUCNTION_RIGHT_PART(x):
    return 1 / (16 * np.pi) * (5 + 3 * np.cos(2*x))

def FUNCTION_TRUE(x):
    return 1 / (160*np.pi) * (25 + 27 * np.cos(2*x))

a = 0
b = 2 * np.pi
koef = -1/(4*np.pi)
Epsilon = 10**-12
Epsilon_1 = 10**-6
Epsilon_0 = 10**-6
m = 1
n = 3*m + 1
h = (b-a)/(n-1)
koef_lambda = 3/8*h*koef
p = 6   
cur_Epsilon = np.inf
plot = True
img_save_dir = "C:\\Users\\ukhin\\Desktop\\Study\\ЧМы\\task_2\\images\\"
test_name = "TEST 13"
msg = f"Right part is 1/(16*pi)*(5 + 3*cos(2*x))\nK is 1/(((sin((x+t)/2)**2))+0.25*(cos((x+t)/2)**2))\n [A, B] = [{a}, {b}]\nEpsilon = {Epsilon}\n"

run(msg, test_name, a, b, n, p, koef, Epsilon, Epsilon_1, Epsilon_0, m, koef_lambda, cur_Epsilon, plot, FUNCTION_K, FUCNTION_RIGHT_PART, FUNCTION_TRUE, img_save_dir)
