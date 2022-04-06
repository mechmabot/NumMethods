import numpy as np
from matplotlib import pyplot as plt
#import sys


def Koefs_formulas(n):
    Coefs = np.zeros(n)
    Coefs[0] = 1
    Coefs[-1] = 1
    for i in range(4, n-2, 3):
        Coefs[i-1] = 2
    for i in range(2, n-1, 3):
        Coefs[i-1] = 3
        Coefs[i] = 3
    return Coefs

def Matrix_to_solve(a, b, n, FUNCTION_K): # n = 3m + 1
    matrix = []
    h = (b-a)/(n-1)
    X = [a + i*h for i in range(n)]
    S = X.copy()
    Coefs = Koefs_formulas(n)
    for i in range(n):
        matrix.append([
            Coefs[j] * FUNCTION_K(X[i], S[j]) for j in range(n)
        ])
    return np.array(matrix)

def Right_part_to_solve(a, b, n, FUCNTION_RIGHT_PART):

    h = (b-a)/(n-1)

    X = [a + i*h for i in range(n)]

    F = [FUCNTION_RIGHT_PART(x) for x in X]

    return np.array(F)



def Solve_system(M, F, n, k_lmbd):
    E = np.eye(n)
    #U = np.dot(np.linalg.inv(E - k_lmbd*M),F)
    U = np.linalg.solve(E - k_lmbd*M, F)
    return U


def Calculate_X(a,b,U_n,x,FUCNTION_RIGHT_PART, FUNCTION_K, koef) :
    n = U_n.shape[0]
    h = (b-a)/(n-1)
    S = np.array([a + i*h for i in range(n)])

    A = Koefs_formulas(n)

    f = FUCNTION_RIGHT_PART(x)

    t1 = np.array([FUNCTION_K(x,S[j]) for j in range(len(S))])
    t2 = t1 * A * U_n * 3 / 8 * h

    return t2.sum()*koef + f


# Вычисление разности значений функций в точке 'x' для двух разных U
def Difference_U(a,b,U,U_prev, x,  FUCNTION_RIGHT_PART,FUNCTION_K, koef):
    Ux = Calculate_X(a,b,U,x,  FUCNTION_RIGHT_PART=FUCNTION_RIGHT_PART,FUNCTION_K=FUNCTION_K, koef=koef)
    Ux_prev = Calculate_X(a,b,U_prev,x,  FUCNTION_RIGHT_PART=FUCNTION_RIGHT_PART,FUNCTION_K=FUNCTION_K, koef=koef)
    return Ux - Ux_prev

# Интеграл Симпсона для разностной функции двух восстановленных
# функций по U
def Integrate_difference(t,U,U_prev,a_ist, b_ist,  FUCNTION_RIGHT_PART,FUNCTION_K, koef):
    h1 = (t[1]-t[0])/2
    S = h1/3 * (Difference_U(a_ist, b_ist,U,U_prev,t[0],  FUCNTION_RIGHT_PART,FUNCTION_K, koef) + 4 * Difference_U(a_ist, b_ist,U,U_prev,(t[0]+t[1])/2,  FUCNTION_RIGHT_PART,FUNCTION_K, koef) + Difference_U(a_ist, b_ist,U,U_prev,t[1],  FUCNTION_RIGHT_PART,FUNCTION_K, koef))
    return abs(S)

# Вычисление разности значений функций в точке 'x' для восстановленой
# по U и реальной функцией
def Error(a,b,U, x, FUNCTION_TRUE, FUCNTION_RIGHT_PART,FUNCTION_K, koef):
    Ux = Calculate_X(a,b,U,x,  FUCNTION_RIGHT_PART,FUNCTION_K, koef)
    Real = FUNCTION_TRUE(x)
    return Ux - Real

# Интеграл Симпсона для разностной функции восстановленной
# функций по U и реальной функции
def Integrate_error(t,U,a_ist, b_ist, FUNCTION_TRUE, FUCNTION_RIGHT_PART, FUNCTION_K, koef):
    h1 = (t[1]-t[0])/2
    S = h1/3 * (Error(a_ist, b_ist,U,t[0], FUNCTION_TRUE, FUCNTION_RIGHT_PART,FUNCTION_K, koef) + 4 * Error(a_ist, b_ist,U,(t[0]+t[1])/2, FUNCTION_TRUE,FUCNTION_RIGHT_PART,FUNCTION_K, koef) + Error(a_ist, b_ist,U,t[1], FUNCTION_TRUE,FUCNTION_RIGHT_PART,FUNCTION_K, koef))
    return abs(S)

# Рассчет интеграла адаптивным методом "слева-направо"
def Calculate_integral(a,b,function, p, Epsilon_1, Epsilon_0):
    h = (b-a) / (100)
    res = 0
    a1 = a
    Ih = 0
    Ih2 = 0
    h1 = h
    while a1 < b:
        while True:
            Ih = function((a1,a1+h1))
            Ih2 = (function((a1,a1+h1/2)) + function((a1+h1/2,a1+h1)))
            Current_epsilon = (Ih2 - Ih) / (2**p - 1)
            if abs(Current_epsilon) > Epsilon_1 and abs(Current_epsilon) > Epsilon_0*Ih2:
                h1 /= 2
            else:
                break
        res += Ih
        a1 += h1
        if abs(Current_epsilon)*(2**p) > Epsilon_1 and abs(Current_epsilon)*(2**p) > Epsilon_0*Ih2:
            pass
        else:
            h1 *= 2
        if a1 + h1 > b:
            h1 = b - a1
    return res




def run(starting_msg, test_name,a, b, n, p, koef, Epsilon, Epsilon_1, Epsilon_0, m, koef_lambda, Current_epsilon, plot, FUNCTION_K, FUCNTION_RIGHT_PART, FUNCTION_TRUE, img_save_dir):
    print(test_name)
    print(starting_msg)
    M = Matrix_to_solve(a, b, n, FUNCTION_K=FUNCTION_K)  
    F = Right_part_to_solve(a, b, n, FUCNTION_RIGHT_PART=FUCNTION_RIGHT_PART)   
    U = Solve_system(M, F, n, koef_lambda)    
    num_of_iters = 1   
    while Current_epsilon > Epsilon:   
        U_prev = U.copy()
        n_prev = n
        Current_epsilon_prev = Current_epsilon    
        m *= 2
        n = 3*m + 1
        h = (b-a)/(n-1)
        koef_lambda = 3/8*h*koef   
        M = Matrix_to_solve(a, b, n, FUNCTION_K=FUNCTION_K)   
        F = Right_part_to_solve(a, b, n, FUCNTION_RIGHT_PART=FUCNTION_RIGHT_PART)    
        U = Solve_system(M, F, n, koef_lambda)   
        Current_epsilon = np.sqrt(Calculate_integral (a, b, lambda x: Integrate_difference(x, U, U_prev, a, b,  FUCNTION_RIGHT_PART,FUNCTION_K, koef)**2, p, Epsilon_1, Epsilon_0))   
        if plot:
            X = [a + h*i for i in range(n)]
            Ux = np.array([Calculate_X(a,b,U,x, FUCNTION_RIGHT_PART, FUNCTION_K, koef) for x in X])
            Ux_prev = np.array([Calculate_X(a,b,U_prev,x, FUCNTION_RIGHT_PART=FUCNTION_RIGHT_PART,FUNCTION_K=FUNCTION_K, koef=koef) for x in X])
            plt.figure(figsize=(12, 4))
            plt.plot(X, Ux-Ux_prev, label='Un - U_n-1 от x', lw=2)
            #plt.plot(X, Ux_prev, label='Un-1', linestyle='--', lw=2)
            #plt.ylim(Ux.min() - 1, Ux.max() + 1)
            plt.legend()
            #plt.show()
            plt.savefig(img_save_dir+f"{test_name}_iter_{num_of_iters}.png")      
        print(f"########## Итерация: {num_of_iters} ##########")
        if Current_epsilon > Epsilon:
            print('Текущее Epsilon = {} \nТекущее Epsilon > Epsilon \n{} > {} \nколичество разбиений = {}'.format(Current_epsilon, Current_epsilon, Epsilon, n))
        else:
            print('Текущее Epsilon = {} \nТекущее Epsilon < Epsilon \n{} < {} \nколичество разбиений = {}'.format(Current_epsilon, Current_epsilon, Epsilon, n))  
        print('\n')  
        num_of_iters += 1   
    print('Общее количество итераций = {}'.format(num_of_iters))   
    h = (b-a)/1000
    X = np.array([a + h*i for i in range(1000)])
    real_f = np.array([FUNCTION_TRUE(x) for x in X])
    approx = np.array([Calculate_X(a,b,U,x,  FUCNTION_RIGHT_PART, FUNCTION_K, koef) for x in X])
    plt.figure(figsize=(12,4))
    plt.plot(X, real_f - approx, label = 'График ошибки по координатам', lw = 2)
    #plt.plot(X, approx, label = 'Численная аппроксимация решения', linestyle = '--', lw = 2)
    #if test_variant == 19:
    #	plt.ylim(approx.min() - 1 ,approx.max() + 1)
    plt.legend()
    #plt.show()
    plt.savefig(img_save_dir+f"{test_name}_iter_{num_of_iters}.png")
    print('L2 норма ошибки')
    L2_norm = np.sqrt(Calculate_integral(a, b, lambda x: Integrate_error(x, U, a,b, FUNCTION_TRUE ,FUCNTION_RIGHT_PART=FUCNTION_RIGHT_PART,FUNCTION_K=FUNCTION_K, koef=koef)**2, p, Epsilon_1, Epsilon_0))
    print('L2 норма = {}\n'.format(L2_norm))
    
    print('С норма ошибки')
    C_norm = abs(real_f.reshape(-1,) - approx.reshape(-1,)).max()
    print('C норма на [{},{}] = {}'.format(a,b,C_norm))

