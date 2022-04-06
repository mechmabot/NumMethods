import numpy as np
from matplotlib import pyplot as plt
from  module import *
from sys import argv

def FUNCTION_K(x , t):
    return np.sin(x*t)

def FUCNTION_RIGHT_PART(x):
    return 1 + ((1/x) * (np.cos(x/2) - 1))

def FUNCTION_TRUE(x):
   return 1

a = 10**-15
b = 0.5
koef = 1
Epsilon = 10**-20
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
test_name = "TEST 12"
msg = f"Right part is sin(x*t)\nK is 1+((1/x) * (cos(x/2) - 1))\n [A, B] = [{a}, {b}]\nEpsilon = {Epsilon}\n"

run(msg, test_name, a, b, n, p, koef, Epsilon, Epsilon_1, Epsilon_0, m, koef_lambda, cur_Epsilon, plot, FUNCTION_K, FUCNTION_RIGHT_PART, FUNCTION_TRUE, img_save_dir)
