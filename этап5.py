import numpy as np
import matplotlib.pyplot as plt
import math as m

from scipy.integrate import ode


with open('Этап 4 выход.txt') as f:
    angle_start_0 = float(f.readline())
    w = float(f.readline())
    H = float(f.readline())
    T = float(f.readline())

v = 0                             #начальные условия
G = 6.6743 * 10**(-11)
M_Earth = 5.972 * 10**(24)
R_Earth = 6371000
R_Orbit = 385*10**6
R_Moon = 173800
M_Moon = 7.348*10**22
m_CM = 22500
m_CMf = 17700 
m_KO = 5500
m_dry = m_CM - m_CMf + m_KO
m_fuel = m_CMf - 11242        
M = m_dry + m_fuel

angular_moon_speed = 1023/R_Orbit

H = H - R_Moon

Moon_angle_start = 0
Vorb = (G*M_Moon/(H + R_Moon))**0.5
angle_start = 28.72*m.pi/180          
x_start = R_Orbit + m.cos(angle_start)*(R_Moon + H)
y_start = m.sin(angle_start)*(R_Moon + H)
Vx_start = -Vorb*m.cos(angle_start + 0.5*m.pi)
Vy_start = 1023 - Vorb*m.sin(angle_start + 0.5*m.pi)

ForceCM = 95750
V_reactCM = 3050

Hmin = R_Orbit

def fout1(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    if ForceCM*t/V_reactCM >= m1:
        return -1
    
def fout2(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    if y[0]<0 and y[2]<0 or (y[0]*y[0]+y[2]*y[2])**0.5 < 1.01*R_Earth:        
        return -1
    
def fout3(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    
    if ForceCM*t/V_reactCM + m1 >= m_fuel:
        return -1
    
    
def fout4(t, y):
    ts.append(t)
    ys.append(list(y.copy()))
    if (y[0]*y[0]+y[2]*y[2])**0.5 < R_Earth + 71000:
        return -1
    

def f1(t, y):

    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    
    dm_fuel = ForceCM*t/V_reactCM
    a = ForceCM/(M - dm_fuel)
    angle = angle_start
 
    ax = -1 * G * M_Earth * y1 / (y1**2 + y3**2)**(1.5) - G * M_Moon * (y1m) / (y1m**2 + y3m**2)**(1.5) - a*m.cos(angle + 0.5*m.pi)
    ay = -1 * G * M_Earth * y3 / (y1**2 + y3**2)**(1.5) - G * M_Moon * (y3m) / (y1m**2 + y3m**2)**(1.5) - a*m.sin(angle + 0.5*m.pi)

    return [y2, ax, y4, ay, angular_moon_speed]

def f2(t, y):
    
    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    

    
    ax = -1 * G * M_Earth * y1 / (y1**2 + y3**2)**(1.5) - G * M_Moon * (y1m) / (y1m**2 + y3m**2)**(1.5)
    ay = -1 * G * M_Earth * y3 / (y1**2 + y3**2)**(1.5) - G * M_Moon * (y3m) / (y1m**2 + y3m**2)**(1.5)

    return [y2, ax, y4, ay, angular_moon_speed]

def f3(t, y):
    global Hmin, V, Time
    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    
    if (y1*y1 + y3*y3)**0.5 - R_Earth < Hmin + 1 and (y1*y1 + y3*y3)**0.5 - R_Earth > 70000:
        Hmin = (y1*y1 + y3*y3)**0.5 - R_Earth
        V = (y2*y2 + y4*y4)**0.5
        Time = y5/(2*m.pi)*27.32*24*60*60
    
    dm_fuel = ForceCM*t/V_reactCM + m1
    a = - ForceCM/(M - dm_fuel)
    angle = m.pi#*135/180
    
    ax = -1 * G * M_Earth * y1 / (y1*y1 + y3*y3)**(1.5) - G * M_Moon * (y1m) / (y1m**2 + y3m**2)**(1.5) + a*m.cos(angle + 0.5*m.pi)
    ay = -1 * G * M_Earth * y3 / (y1*y1 + y3*y3)**(1.5) - G * M_Moon * (y3m) / (y1m**2 + y3m**2)**(1.5) + a*m.sin(angle + 0.5*m.pi)

    return [y2, ax, y4, ay, angular_moon_speed]

def f4(t, y):
    global Hmin, V, Time
    
    y1, y2, y3, y4, y5 = y
    y1m, y3m = y1 - m.cos(y5)*R_Orbit, y3 - m.sin(y5)*R_Orbit
    
    if (y1*y1 + y3*y3)**0.5 - R_Earth <= Hmin + 1 and (y1*y1 + y3*y3)**0.5 - R_Earth > 70000:
        Hmin = (y1*y1 + y3*y3)**0.5 - R_Earth
        V = (y2*y2 + y4*y4)**0.5
        Time = y5/(2*m.pi)*27.32*24*60*60
    
    ax = -1 * G * M_Earth * y1 / (y1*y1 + y3*y3)**(1.5) - G * M_Moon * (y1m) / (y1m**2 + y3m**2)**(1.5)
    ay = -1 * G * M_Earth * y3 / (y1*y1 + y3*y3)**(1.5) - G * M_Moon * (y3m) / (y1m**2 + y3m**2)**(1.5)

    return [y2, ax, y4, ay, angular_moon_speed]

def facc(t, y):
   
    dm_fuel = ForceCM*t/V_reactCM
    a = ForceCM/(M - dm_fuel)
    
    return a

time1 = 800000
time2 = 90000 #время
m1 = 4046.5


Time_waiting = (H+R_Moon)/Vorb * (m.pi - angle_start_0 - angle_start)


xc,yc=[],[]                                #Табличная функция для круга Земли
for i in range(0, 630):
    xc.append(R_Earth*m.cos(i/100))
    yc.append(R_Earth*m.sin(i/100))
    

    
#этап 1, ускорение    


y0,t0=[x_start,  Vx_start, y_start, Vy_start, Moon_angle_start], 0 # начальные условия для интегрирования
print(y0)
ODE = ode(f1)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 0.5)
ODE.set_solout(fout1)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time1)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//15 + 1

#graph(Y, period)        #Вывод графиков
print(len(Y))
#этап 2, свободный полет

y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3],Y[-1:,4]], 0 # начальные условия для интегрирования

ODE = ode(f2)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 10)
ODE.set_solout(fout2)
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time1)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//20

#graph(Y, period)      #Вывод графиков 
print(len(Y))
#Этап 3, торможение

y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3],Y[-1:,4]], 0

ODE = ode(f3)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 1)
ODE.set_solout(fout3)
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time1)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//5 + 1

print(len(Y))

y0,t0=[Y[-1:,0],Y[-1:,1],Y[-1:,2],Y[-1:,3],Y[-1:,4]], 0 # начальные условия для интегрирования

ODE = ode(f4)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 2)
ODE.set_solout(fout4)
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time2)      # решение ОДУ
Y=np.array(ys)
period = len(Y)//10   #Период - шаг, с которым проходятся элементы списка Y с проинтегрированными значениями
#Чем меньше период, тем больше графиков выведется, если период равен len(Y), будет один график
def graph(Y, period):                                            #Функция вывода графиков
    for k in range(1, len(Y), period):                       #Период - шаг, с которым проходятся элементы списка Y с проинтегрированными значениями
        if (Y[k, 0]**2 + Y[k, 2]**2)**0.5 < 70000000:               #Одно из трех условий для вывода разных видов графиков
            
            Y1 = Y[30000:k]
            plt.plot(Y1[:,0],Y1[:,2],linewidth=2)
            plt.axis('equal')
            plt.plot(xc,yc,'black',linewidth=1)
        
            plt.title("Подлет к Земле\n ")
            plt.grid(True)
            plt.show()
        
        elif ((Y[k, 0])**2 + (Y[k, 2])**2)**0.5 > 330000000:
            if k < 50000:
                Y1 = Y[0:k]
            else:
                Y1 = Y[k-50000:k]
            
            Y2 = []
            Y3 = []
            
            for i in range(len(Y1)):
                Y2.append(m.cos(Y1[i,4])*R_Orbit)
                Y3.append(m.sin(Y1[i,4])*R_Orbit)
            plt.plot(Y2[:], Y3[:], '#3589EF', linewidth = 1)
            plt.plot(Y1[:,0],Y1[:,2],linewidth=2)
            plt.axis('equal')
            Moon_x, Moon_y = m.cos(Y1[-1,4])*R_Orbit, m.sin(Y1[-1,4])*R_Orbit
            
            xn,yn = [], []
            for i in range(0, 630):
                xn.append(1737000*m.cos(i/100) + Moon_x)
                yn.append(1737000*m.sin(i/100) + Moon_y)
            plt.plot(xn,yn,'black',linewidth=2)
            
            plt.title("Flight near Moon\n ")
            
            plt.grid(True)
            
            plt.show()
        else:
            Y1 = Y[:k]
            plt.plot(Y1[:,0],Y1[:,2],linewidth=2)
            
            Y2 = []
            Y3 = []
            
            for i in range(len(Y1)):
                Y2.append(m.cos(Y1[i,4])*R_Orbit)
                Y3.append(m.sin(Y1[i,4])*R_Orbit)
            plt.plot(Y2[:], Y3[:], '#3589EF', linewidth = 1)
            plt.axis('equal')
            plt.plot(xc,yc,linewidth=1)
            Moon_x, Moon_y = m.cos(Y1[-1,4])*R_Orbit, m.sin(Y1[-1,4])*R_Orbit
            
            xn,yn = [], []
            for i in range(0, 630):
                xn.append(1737000*m.cos(i/100) + Moon_x)
                yn.append(1737000*m.sin(i/100) + Moon_y)
            plt.plot(xn,yn,linewidth=2)
              
            plt.title("Полет к Земле\n ")
            plt.grid(True)
            plt.show()
graph(Y, period) 
print(len(Y))

#расчет скорости отлета

y0,t0=Vorb, 0 # начальные условия для интегрирования
ODE = ode(facc)
ODE.set_integrator('dopri5', nsteps=200000, max_step = 0.5)
ODE.set_solout(fout1)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений
ODE.integrate(time1)
Y=np.array(ys)      # решение ОДУ
v = Y[-1]
print('Скорость отлета = ', v)
print('Hmin (km) = ', Hmin/1000)
print('V = ', V, Time+Time_waiting)


            
f = open('Этап 5 выход.txt', 'w')
f.write(str(Hmin)+'\n')
f.write(str(V)+'\n')
f.write(str(T + Time + Time_waiting))
f.close()