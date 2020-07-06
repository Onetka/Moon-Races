import math
from scipy.integrate import ode
import matplotlib.pyplot as plt

M_Lk = 28000# масса лунного корабля  
M_Lm = 15000# масса лунного модуля
mRN1 = 135000# сухая масса первой ступени
mRN1F = 2145000 - mRN1  # масса топлива в первой ступени
RN1force = 34350000 # сила тяги первой ступени 
Vreact1 = 2580 # скорость истечения продуктов сгорания первой ступени
d1 = 10.1# диаметр первой ступени
mRN2 = 37600# сухая масса второй ступени
mRN2F = 458700 - mRN2 # масса топлива во второй ступени
RN2force = 5115000  # сила тяги второй ступени 
Vreact2 = 4130 # скорость истечения продуктов сгорания второй ступени
d2 = 10.1
mRN3 = 12000 # сухая масса третей ступени
mRN3F = 120000 - mRN3 # масса топлива в третей ступени
RN3force = 1016000 # сила тяги третей ступени
Vreact3 = 4130  # скорость истечения продуктов сгорания третей ступени
d3 = 6.6
M_total = mRN1 + mRN2 + mRN1F + mRN2F + mRN3 + mRN3F + M_Lk + M_Lm
Me = 5.972 * 10 ** 24  # kg  Mass of the Earth
We = 7.29 * 10 ** (-5)
G = 6.67408 * 10 ** (-11)
R = 6375000  # m
h = 350000
vt_ = []
vr_ = []
xs = []
ys = []
xsf = []
ysf = []
h_ = []
t_ = []


def step1(t, y):
    global xs, ys, RN1force, M_total, mRN1F, Vreact1
    xs.append(y[0])
    ys.append(y[2])
    h_.append(math.sqrt(y[0]**2 + y[2]**2) - R)
    vr_.append(- (y[1]*math.cos(math.pi - math.atan(y[2]/y[0])) + y[3]*math.cos(math.pi - math.atan(y[0]/y[2]))))
    vt_.append(math.sqrt(y[1]**2 + y[3]**2 - vr_[-1] **2 ))
    t_.append(t)
    if mRN1F - RN1force/Vreact1 * t <= 0.001:
        M_total -= mRN1 + mRN1F
        print("First stage detached")


def step2(t, y):
    global xs, ys, RN2force, M_total, Vreact2
    xs.append(y[0])
    ys.append(y[2])
    h_.append(math.sqrt(y[0]**2 + y[2]**2) - R)
    vr_.append(- (y[1]*math.cos(math.pi - math.atan(y[2]/y[0])) + y[3]*math.cos(math.pi - math.atan(y[0]/y[2]))))
    vt_.append(math.sqrt(y[1] ** 2 + y[3] ** 2 - vr_[-1] ** 2))
    t_.append(t)
    if mRN2F - RN2force/Vreact2 * t <= 0.001:
        M_total -= mRN2 + mRN2F
        print("Second stage detached")


def step3(t, y):
    global xs, ys, M_total, Vreact3
    xs.append(y[0])
    ys.append(y[2])
    h_.append(math.sqrt(y[0]**2 + y[2]**2) - R)
    vr_.append(- (y[1]*math.cos(math.pi - math.atan(y[2]/y[0])) + y[3]*math.cos(math.pi - math.atan(y[0]/y[2]))))
    vt_.append(math.sqrt(y[1] ** 2 + y[3] ** 2 - vr_[-1] ** 2))
    t_.append(t)
    h = math.sqrt(y[0]**2 + y[2]**2)


def step4(t, y):
    global xs, ys, RN2force, M_total, Vreact2
    xsf.append(y[0])
    ysf.append(y[2])

#
def angle(x, y):
    angle = 90
    h = math.sqrt(x*x + y*y) - R
    H = (10000, 12000, 14000, 16000, 18000, 20000, 22000, 24000, 70000, 75000, 80000, 90000, 100000, 110000, 185000, 198000)
    ang = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85)
    for i in range(16):
        if h <= H[i]:
            angle = ang[i]
            break
    return math.radians(angle)

def ro(x, y):
  Space_lst = [0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 12000, 14000, 16000, 18000, 20000, 24000, 28000, 32000, 36000, 40000, 50000, 60000, 70000]
  ro_lst = [1.2250, 1.1673, 1.1117, 1.0581, 1.0065, 0.9569, 0.9093, 0.8194, 0.7365, 0.6601, 0.59, 0.5258, 0.4671, 0.4135, 0.3119, 0.2279, 0.1665, 0.1216, 0.0889, 0.0469, 0.0251, 0.0136, 0.00726, 0.004, 0.00103, 0.0003, 0]
  h1 = math.sqrt(x*x + y*y)
  i = 0
  while (h1 >= Space_lst[i]):
    i += 1
    if i > 25:
      break
  return ro_lst[i]

def launch1(t, y):  # Запуск первой ступени
    global Vreact1, mRN1, mRN1F, M_total, RN1force, d1
    X, X1, Y, Y1 = y
    ang = angle(X, Y) + math.atan(-Y/X)
    S = math.pi * d1 ** 2 / 4  # Параметры для силы сопротивления
    C = 0.1
    v2 = X1**2 + Y1**2
    D = ro(X, Y)
    dM_fuel = RN1force/Vreact1 * t
    a = RN1force/(M_total - dM_fuel)
    res = C*S*v2*D / 2     # сила сопротивления

    ax = - X * G * Me / ((X*X + Y*Y) ** 1.5) + (a - res) * math.cos(ang)
    ay = - Y * G * Me / ((X*X + Y*Y) ** 1.5) - (a - res) * math.sin(ang)

    return [X1, ax, Y1, ay]


def launch2(t, y):  # Запуск второй ступени
    global Vreact2, mRN2, mRN2F, M_total, RN2force, d2
    X, X1, Y, Y1 = y
    ang = angle(X, Y) + math.atan(-Y/X)
    S = math.pi * d2 ** 2 / 4
    C = 0.1
    v2 = X1**2 + Y1**2
    D = ro(X, Y)

    dM_fuel = RN2force/Vreact2 * t
    a = RN2force/(M_total - dM_fuel)

    res = C*S*v2*D / 2
    ax = - X * G * Me / ((X*X + Y*Y) ** 1.5) + (a - res) * math.cos(ang)
    ay = - Y * G * Me / ((X*X + Y*Y) ** 1.5) - (a - res) * math.sin(ang)

    return [X1, ax, Y1, ay]


def launch3(t, y):  # Запуск третьей ступени
    global Vreact3, mRN3, mRN3F, M_total, RN3force, d3
    X, X1, Y, Y1 = y
    ang = angle(X, Y) + math.atan(-Y/X)
    S = math.pi * d3 ** 2 / 4  # Параметры для силы сопротивления
    C = 0.1
    v2 = X1**2 + Y1**2
    D = ro(X, Y)

    dM_fuel = RN3force/Vreact3 * t
    a = RN3force/(M_total - dM_fuel)
    res = C*S*v2*D / 2     # сила сопротивления

    ax = - X * G * Me / ((X*X + Y*Y) ** 1.5) + (a - res) * math.cos(ang)
    ay = - Y * G * Me / ((X*X + Y*Y) ** 1.5) - (a - res) * math.sin(ang)

    return [X1, ax, Y1, ay]


def free_flight(t, y):
    X, X1, Y, Y1 = y
    ax = - X * G * Me / ((X * X + Y * Y) ** 1.5)
    ay = - Y * G * Me / ((X * X + Y * Y) ** 1.5)
    return [X1, ax, Y1, ay]


y0 = [R, 0, 0, -We * R]  # Начальные условия
t0 = 0
0
f = ode(launch1)
f.set_integrator('dopri5', nsteps=100000)
f.set_initial_value(y0, t0)
f.set_solout(step1) # функция, которая вызывается после каждой итерации

y1 = f.integrate(mRN1F*Vreact1 / RN1force)  # решение системы для первой ступени; аргумент - время, возвращает список: x, x1, y, y1
t1 = mRN1F*Vreact1 / RN1force


f = ode(launch2)
f.set_integrator('dopri5', nsteps=100000)
f.set_initial_value(y1, t1)
f.set_solout(step2) # функция, которая вызывается после каждой итерации

y2 = f.integrate(mRN2F*Vreact2 / RN2force)  # решение для второй ступени
t2 = mRN2F*Vreact2 / RN2force


f = ode(launch3)
f.set_integrator('dopri5', nsteps=100000)
f.set_initial_value(y2, t2)
f.set_solout(step3) # функция, которая вызывается после каждой итерации

# y3 = f.integrate(t2 + (mRN3F * 0.35) * Vreact3 / RN3force) 572.8
y3 = f.integrate(505.5)
t3 = 505.5

print("free flight")
print(y3[0], y3[2], "Height: ",math.sqrt(y3[0] ** 2 + y3[2] ** 2) - R, ", Velocity: ", math.sqrt(y3[1] ** 2 + y3[3] ** 2), ", Fuel left: ", (mRN3F - RN3force/Vreact3 * (t3 - t2)) / mRN3F, sep="")

f = ode(free_flight)
f.set_integrator('dopri5')
f.set_initial_value(y3, t3)
f.set_solout(step4)
f.integrate(10000)


x, vx, y, vy = y3
dm = (mRN3F - RN3force/Vreact3 * (t3 - t2)) / mRN3F
h = math.sqrt(x**2 + y**2)
nx, ny = math.cos(math.atan(-y/x)), math.sin(math.atan(-y/x))
w = math.sqrt((ny * vx + nx * vy)**2) / h
phi = math.atan(-y/x)
f = open("BlastOffFromEarth.txt", 'w')
print(dm, w, h, t3, phi, file=f)
f.close()

xc, yc = [], [] # рисует окружность Земли для графика
for i in range(0, 630):
    xc.append(R*math.cos(i/100))
    yc.append(R*math.sin(i/100))

with plt.style.context("seaborn-white"):
    graph = plt.figure()
    graph.set_facecolor("white")
    graph = plt.gca()
    graph.set_facecolor("white")
    plt.plot(xc, yc, color='black', linewidth="2")
    plt.axis('equal')
    plt.plot(xs, ys, linewidth='1.5', color="#8B0000")
    plt.plot(xsf, ysf, color='#0000FF', alpha=0.8, linewidth='1')
    plt.title("Выход на Низкую Околоземную орбиту", fontsize="15")
    plt.show()
