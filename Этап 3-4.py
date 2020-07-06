import  math as m
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import matplotlib.transforms as trns
import pylab
import numpy as np
import os

G = 6.67e-11 
Mm = 7.36e22  
Rm = 1738000  
m_ship = 6835 #масса сухого лунного модуля
mVSLM = 2315 # ВСЛМ сухая
mPSLMf = 8165 # масса топлива в ПСЛМ
mVSLMf = 2355
PSLMforce = 45040
VSLMforce = 15600
h_0_Moon = 2
Vreact = 3050
dt = 0.01

X, Y, R, V, A, T = [], [], [], [], [], []

f=np.loadtxt("Этап 2 выход.txt", delimiter=' ', dtype=np.float)
x_Moon, y_Moon, x_ship, y_ship, a, b, FullT = f



def change_coordinates (x_Moon, y_Moon, x_ship, y_ship):
    New_x =(x_ship - x_Moon)*(y_Moon/m.sqrt(x_Moon**2 + y_Moon**2)) - (y_ship - y_Moon)*(x_Moon/m.sqrt(x_Moon**2 + y_Moon**2))# пересчет координат из ВСК в систему координат относительно Луны
    New_y =(x_ship - x_Moon)*(x_Moon/m.sqrt(x_Moon**2 + y_Moon**2)) + (y_ship - y_Moon)*(y_Moon/m.sqrt(x_Moon**2 + y_Moon**2))
    R_orbit = m.sqrt(New_x**2 + New_y**2)*1000
    V_orbit = m.sqrt(G * Mm / R_orbit)
    return New_x, New_y, R_orbit, V_orbit

new_x, new_y, R_orbit, V_orbit = change_coordinates (x_Moon, y_Moon, x_ship, y_ship)

def gg(x, y, G, Mm):  
    gx = - (G*Mm*x)/((x**2 + y**2)**1.5)
    gy = - (G*Mm*y)/((x**2 + y**2)**1.5)
    return gx, gy

def current_position(x, y, r):  # местоположение лунного корабля на орбите
    if x <= 0:
        alpha = m.acos(y/r)
    else:
        alpha = m.acos(-y/r) + m.pi
    position = r*alpha
    return position

def newData(x, y, vx, vy, ax, ay):
    x = x + vx*dt + ax*(dt**2)/2
    y = y + vy*dt + ay*(dt**2)/2
    vx = vx + ax*dt
    vy = vy + ay*dt
    return x, y, vx, vy

def free(x, y, vx, vy):  # свободный полет
    while x <= 0:
        ax, ay = gg(x, y, G, Mm)
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
    return x, y, vx, vy

def desired_x(x=0, y=R_orbit, vx=-V_orbit, vy=0):  # высота на которой начинается торможение 
    #speed = vx
    x_i = x
    y_i = y
    vy_i = vy
    height = abs(y)
    while height > Rm + ((y_i - Rm) / 2):
        vx += 0.1
        x = x_i
        y = y_i
        vy = vy_i
        x, y, vx, vy = free(x, y, vx, vy)
        height = abs(y)
    desiredx = (-vx**2 + V_orbit**2)/(2*PSLMforce/(m_ship + mPSLMf))
    return desiredx

#летаем до нужной точки, а потом начинаем тормозить
def after_impulse(x_ship_orbit, y_ship_orbit, mfuel = mPSLMf): # скорость и координаты перед началом свободного полёта
    global X, Y, R, V, A, T 
    x = desired_x()
    y = m.sqrt(R_orbit**2 - x**2)
    decel_startco = current_position(x, y, R_orbit)
    rocket_co = current_position(x_ship_orbit, y_ship_orbit, R_orbit)
    if decel_startco >= rocket_co: # летим до нужных нам координат, с которых начинаем торможение
        t_1stage = (decel_startco - rocket_co) / V_orbit
    else:
        t_1stage = (decel_startco + 2 * m.pi * R_orbit - rocket_co) / V_orbit
    time = 0
    vx = - (V_orbit*y)/R_orbit
    vy = (V_orbit*x)/R_orbit
    ax = PSLMforce/(m_ship + mPSLMf)
    while x > 0:
        gx, ay = gg(x, y,  G, Mm)
        afuel = ax - gx
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
        mfuel -= ((m_ship + mfuel)*afuel/Vreact)*dt
        time += dt
        X.append(x/1000), Y.append(y/1000), R.append(m.sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
        A.append(m.sqrt(ax ** 2 + ay ** 2)), T.append(time)
    return mfuel, vx, vy, x, y, t_1stage, time, decel_startco

mfuel_aftdecel, vx, vy, x, y, t_1stage, time, decel_startco = after_impulse(new_x, new_y)

print('Конец торможения, начало свободного полёта')
print("x =", round(x, 4), "y =", round(-y, 4), "Vx =", round(vx, 4), "Vy =", round(-vy, 4), "масса топлива =",
      round(mfuel_aftdecel, 4))
#функция вывода координат
def writing_coordinates(x, y, vx, vy, mfuel, time, f_max=PSLMforce):  # свободный полёт
    global X, Y, R, V, A, T
    n = 1
    while y > 0:  # Половина траектории
        ax, ay = gg(x, y, G, Mm)
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
        time += dt
        n += 1
        if n % 10 == 0:
            X.append(x/1000), Y.append(y/1000), R.append(m.sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
            A.append(m.sqrt(ax ** 2 + ay ** 2)), T.append(time)
    ax, ay = gg(x, y, G, Mm)
    while (vx ** 2) < abs(2 * x * (- (f_max/(m_ship + mfuel)) + ax)):  # скорость при которой остаётся горизонтальная составляющая V
        ax, ay = gg(x, y, G, Mm)
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
        time += dt
        n += 1
        if n % 10 == 0:
            X.append(x/1000), Y.append(y/1000), R.append(m.sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
            A.append(m.sqrt(ax ** 2 + ay ** 2)), T.append(time)
    return x, y, vx, vy, time


x_1, y_1, vx_1, vy_1, time = writing_coordinates(x, y, vx, vy, mfuel_aftdecel, time)
print("Конец свободного полёта, начало торможения")
print("x =", round(x_1, 4), "y =", round(-y_1, 4), "Vx =", round(vx_1, 4), "Vy =", round(-vy_1, 4))

def stop(x, y, vx, vy, mfuel, time, f_max=PSLMforce):  # торможение
    global X, Y, R, V, A, T
    ax_G, ay_G = gg(x, y,  G, Mm)
    ax = - (PSLMforce/ (mfuel + m_ship)) + ax_G
    while vx > 0:
        ax_G, ay_G = gg(x, y,  G, Mm)
        ax_decel = - ax + ax_G
        ay_decel = - m.sqrt(((f_max / (m_ship + mfuel)) ** 2) - (ax_decel ** 2))
        ay = ay_G + ay_decel
        if mfuel < 7800:
            if vy > 100:
                mfuel -= (f_max /Vreact) * dt
            else:
                ay = ay_G
                mfuel -= (ax_decel * (m_ship + mfuel) / Vreact) * dt
        else:
            mfuel -= (f_max / Vreact) * dt
        time += dt
        x, y, vx, vy = newData(x, y, vx, vy, ax, ay)
        X.append(x/1000), Y.append(y/1000), R.append(m.sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
        A.append(m.sqrt(ax ** 2 + ay ** 2)), T.append(time)
    return x, y, vx, vy, mfuel, time

print("Конец торможения, начало вертикальной посадки")
x_2, y_2, vx_2, vy_2, mfuel_2, time = stop(x_1, y_1, vx_1, vy_1, mfuel_aftdecel,time)
print("x =", round(x_2, 4), "y =", round(-y_2, 4), "Vx =", round(vx_2, 4), "Vy =", round(-vy_2, 4),
      "масса топлива =", round(mfuel_2, 4))


def vertical_stop(x, y, vx, vy, mfuel, time):  # вертикальная посадка
    global X, Y, R, V, A, T
    while vy > 0:
        ax_G, ay_G = gg(x, y,  G, Mm)
        ay = - (vy ** 2 / (2 * abs(abs(y) - Rm)))
        x, y, vx, vy = newData(x, y, vx, vy, ax_G, ay)
        mfuel -= ((ay_G - ay) * (m_ship + mfuel) / Vreact) * dt
        time += dt
        X.append(x/1000), Y.append(y/1000), R.append(m.sqrt(x ** 2 + y ** 2)/1000), V.append(vx ** 2 + vy ** 2)
        A.append(m.sqrt(ax_G ** 2 + ay ** 2)), T.append(time)
    return x, y, vx, vy, mfuel, time

x_3, y_3, vx_3, vy_3, mfuel_3, t_2stage = vertical_stop(x_2, y_2, vx_2, vy_2, mfuel_2, time)
print("Место посадки: x =", round(x_3, 4), "y =", round(-y_3, 4), "Скорость при посадке: v_x =", round(
    vx_3, 4),  "v_y =", round(-vy_3, 4), "Масса топлива =", round(mfuel_3, 4))

t_landing = t_1stage + t_2stage
position_of_lk = (decel_startco + V_orbit * t_2stage)% (2 * m.pi * R_orbit)
y_lk = R_orbit * m.cos(position_of_lk / R_orbit)
x_lk = - R_orbit * m.sin(position_of_lk / R_orbit)  # местоположение ракеты после посадки

print("Местоположения ракеты носителя на орбите: x =", round(x_lk, 4), "y =", round(-y_lk, 4), "Полное время =", round(
    t_landing + FullT, 4))


plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots()
xc,yc=[],[]
for i in range(0, 630):
    xc.append(Rm*0.001*m.cos(i/100))
    yc.append(Rm*0.001*m.sin(i/100))
plt.plot(xc,yc,linewidth=2, c = 'orange')


plt.axis('equal')
plt.plot(list(X), list(Y), marker = "*", c="black", markersize=0.1)
plt.xlabel("Координата x, км")
plt.ylabel("Координата y, км")
ax.set_title("Траектория ракеты вблизи луны")
pylab.xlim(-2000, 2000)
pylab.ylim(-2000, 2000)
plt.grid(False)
plt.savefig("Moon1.png", dpi=300)
plt.show()

fig, ax = plt.subplots()

x, y, Vx, Vy, t_3stage = x_3, -y_3, 0, 0, 0
X, Y, R, V, T = [], [], [], [], []
#все уравнения для взлета
def ready_steady(angle):
    dt = 1
    global x, y, Vx, Vy, mVSLMf, mVSLM, t_3stage, X, Y, R, V, T
    gx = (G * Mm * x) / ((x ** 2 + y ** 2) ** 1.5)
    gy = (G * Mm * y) / ((x ** 2 + y ** 2) ** 1.5)
    m0 = mVSLM + mVSLMf
    m = m0 - (VSLMforce/Vreact) * dt
    Vr = Vreact*m.log(m/m0)
    x = x + Vx * dt + 0.5 * gx * dt ** 2 - m.cos(m.radians(angle)) * ((-m0 / (VSLMforce/Vreact) + dt) * Vr - Vreact*dt)
    y = y + Vy * dt - 0.5 * gy * dt ** 2 - m.sin(m.radians(angle)) * ((-m0 / (VSLMforce/Vreact) + dt) * Vr - Vreact*dt)
    Vx = Vx + gx * dt - m.cos(m.radians(angle)) * Vr
    Vy = Vy - gy * dt - m.sin(m.radians(angle)) * Vr
    mVSLMf = mVSLMf - (VSLMforce/Vreact) * dt
    t_3stage += dt
    X.append(x / 1000), Y.append(y / 1000), R.append(m.sqrt(x ** 2 + y ** 2) / 1000), V.append(m.sqrt(Vx ** 2 + Vy ** 2))
    T.append(t_3stage)


#рассчкт углов
def go(x_0, y_0):
    global x, y, t_3stage, X, Y, R, V, T, mVSLMf
    t_3stage, i = 0, 0
    x = x_0
    y = y_0
    while m.sqrt(x ** 2 + y ** 2)  <= 1740000:
        i += 1
        ready_steady(90 - i)
    while m.sqrt(x ** 2 + y ** 2) <= 1748000:
        ready_steady(38)
    i = 0
    while m.sqrt(x ** 2 + y ** 2) <= 1764100:
        i += 0.25
        ready_steady(38-i)
    while m.sqrt(x ** 2 + y ** 2) <= 1788500:
        ready_steady((m.acos(y/m.sqrt(x ** 2 + y ** 2)))*(180/m.pi))
    i=0
    while i<35:
        i+=1
        ready_steady(270)
    i=0
    while i < 1:
        i += 1
        ready_steady(0)
    V4 = m.sqrt(Vx ** 2 + Vy ** 2)
    H4 = m.sqrt(x ** 2 + y ** 2)
    print("Место стыковки: x =", round(x, 4), "y =", round(y, 4), "Скорость на орбите: v =", round(
    V4, 4), "Высота орбиты =", round(H4, 4), "Масса топлива =", round(mVSLMf, 4))
    Esc_moon = open("Этап 3-4 вывод", 'w')
    Esc_moon.write(str(m.atan(x/y)) + "\n" + str(m.sqrt(Vx ** 2 + Vy ** 2)) + "\n" + str(m.sqrt(x ** 2 + y ** 2))
                   + "\n" + str(t_landing + 3202 + FullT))
    Esc_moon.close()
    drawing(X, Y, R, V, T)
    return (V4, H4, mVSLMf)


def drawing(X, Y, R, V, T ):
    plt.style.use('seaborn-whitegrid')
    fig, ax = plt.subplots()
    xc, yc = [], []
    for i in range(0, 630):
        xc.append(Rm * 0.001 * m.cos(i / 100))
        yc.append(Rm * 0.001 * m.sin(i / 100))
    plt.plot(xc, yc, linewidth=2, c='orange')

    

    plt.plot(list(X), list(Y), marker="*", c="#3589EF", markersize=0.1)
    plt.xlabel("Координата x, км")
    plt.ylabel("Координата y, км")
    ax.set_title("Траектория ракеты вблизи луны")
    pylab.xlim(-75, 365)
    pylab.ylim(1650, 1790)
    plt.grid(False)
    plt.rcParams['examples.directory'] = os.getcwd()
    image_data = cbook.get_sample_data('rocket.png')
    image = plt.imread(image_data)
    im = ax.imshow(image, origin='lower', extent=[327, 335, 1749, 1772])
    trans_data = trns.Affine2D().rotate_deg_around(333, 1761, 75) + ax.transData
    im.set_transform(trans_data)
    plt.savefig("Moon4.png", dpi=300)
    plt.show()

    plt.style.use('seaborn-whitegrid')
    fig, ax = plt.subplots()

t_landing = t_1stage + t_2stage
position_of_lk = (decel_startco + V_orbit * t_2stage) % (2 * m.pi * R_orbit)
y_lk = R_orbit * m.cos(position_of_lk / R_orbit)
x_lk = - R_orbit * m.sin(position_of_lk / R_orbit)


#полет после стыковки
def vzlet():
    position_of_lk1 = (decel_startco + V_orbit * (
    ((R_orbit * (m.asin(-319716.9245 / R_orbit) + 3*m.pi) - decel_startco) / V_orbit))) % (2 * m.pi * R_orbit)
    t_blast = (R_orbit * (m.asin(-319716.9245 / R_orbit) + 3*m.pi) - decel_startco) / V_orbit
    y_lk1 = R_orbit * m.cos(position_of_lk1 / R_orbit)
    x_lk1 = R_orbit * m.sin(position_of_lk1 / R_orbit)
    go(0, Rm)
    print("Ждём на Луне", (t_blast - t_3stage))
    print("Координаты лунного корабля на орбите. Момент стыковки: x =", round(x_lk1, 4), "y =", round(-y_lk1, 4),
          "Полное время =", round(t_landing+t_blast+FullT))

vzlet()

