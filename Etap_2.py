# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 22:45:55 2020

@author: Пользователь
"""

# -*- coding: utf-8 -*-
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
import math
f=np.loadtxt("Этап 1 выход.txt", delimiter=' ', dtype=np.float)
out=open('Этап 2 выход.txt','w')
dm, W, h, t_previous, phi = f
print(phi)

Ne = 39844e10 #гравитационный параметр Земли
Rem = 384405e3 #радиус орбиты Луны
w0 = 265e-8 #угловая скорость Луны
Nm = 489210064e4 #гравитационный параметр Луны
Rm = 1738000 #радиус Луны
mCMF = 17700 #СМ: масса топлива и окислителя
mRN3 = 12000
mLK = 28000
mLM = 15000
LKforce = 95750
alphaM=0.91-0.2827 #начальный угол луны с осью Оу в новых координатах

t0=(math.pi-phi+w0*t_previous+(math.pi/2)-alphaM)/(W-w0) #время свободного вращения до попадания в точку с координатами (х0,у0) в новой системе координат
delta=(math.pi/2)-alphaM+w0*t_previous #угол между прямой, соединяющей центры Земли и Луны, и осью х в ВСК
print(delta)
print(t0)
YYY=[]
#РС
q0=246
Vreact3=4130
#ЛК
q1=31.4
VtLK=3050
q=q0
Vt=Vreact3
F=q*Vt
ON=1 #двигатель включен  
y0=0
x0=-h
uy0=W*h
ux0=0
fuel=dm*102
fuel=fuel*1000+mCMF-1838
'''
print(h)
print(str((fuel-mCMF)/1000)+'%')
'''
mass=fuel+mRN3+25300
YY=[]    
ts=[]
ys,xs=[],[]

def fout2(t, y):# обработчик шага(перелет) 
        ts.append(t)
        YY.append(list(y.copy()))
        y1, y2, y3, y4, y5 = y
        if np.sqrt((Rem*math.sin(y5)-y1)**2+(Rem*math.cos(y5)-y2)**2) <= 4000000:
            return -1
#Обработчики(торможение)
#y1 - координата по x, y2 - координата по y, y3 - скорость по x, y4 - скорость по y, y5 - угол между прямой, соединяющей центр Луны и корабль, и прямой, соединяющей центры Луны и Земли
def fout3Stop(t,y):
        YY.append(list(y.copy()))
        y1, y2, y3, y4, y5 = y
        y_m=y2-Rem*math.cos(y5)
        x_m=y1-Rem*math.sin(y5)
        ux_m=y3-w0*Rem*math.cos(y5)
        uy_m=y4+w0*Rem*math.sin(y5)
        YYY.append([x_m,y_m,ux_m,uy_m])
        
def fout4Stop(t,y):
        ts.append(t)
        YY.append(list(y.copy()))
        y1, y2, y3, y4, y5 = y
        y_m=y2-Rem*math.cos(y5)
        x_m=y1-Rem*math.sin(y5)
        ux_m=y3-w0*Rem*math.cos(y5)
        uy_m=y4+w0*Rem*math.sin(y5)
        r=np.sqrt(x_m**2+y_m**2)
        arr.append(r)
        YYY.append([x_m,y_m,ux_m,uy_m])

    
def fun(t,y): #ODE Земля и Луна 
    y1, y2, y3, y4, y5 = y
    u=np.sqrt(y4**2+y3**2)
    return [y3,y4,((-1)*Ne*y1)/((np.sqrt(y1**2+y2**2))**3)+
            (F*ON*(y3/u))/(mass-q*t)-((Nm*(y1-Rem*math.sin(y5)))
            /((np.sqrt((y1-Rem*math.sin(y5))**2+(y2-Rem*math.cos(y5))**2))**3)),
    ((-1)*Ne*y2)/((np.sqrt(y1**2+y2**2))**3)+
    (F*ON*(y4/u))/(mass-q*t)-((Nm*(y2-Rem*math.cos(y5)))
            /((np.sqrt((y1-Rem*math.sin(y5))**2+(y2-Rem*math.cos(y5))**2))**3)),w0]

#Разгон 
t1=(fuel-mCMF)/q
Y0, t0=[x0, y0, ux0, uy0, alphaM], 0
ODE=ode(fun)
ODE.set_integrator('dopri5',nsteps=10000)#, max_step=0.01)
ODE.set_solout(fout2)
ODE.set_initial_value(Y0, t0) 
ODE.integrate(t1)      
print()
Y=np.array(YY)
ys,xs=list(Y[:,1]),list(Y[:,0])

#Сброс РС
mass=mLM + mLK
fuel=mCMF
u=VtLK
q=q1
F=LKforce
q=F/VtLK
x0,y0,ux0,uy0,alphaM=list(Y[-1])

#Свободный полет к луне
t2=450000
ON=0 # выключен двигатель
ym,xm=[],[]
YY=[]
ts=[]
Y0, t0=[x0,y0,ux0,uy0,alphaM], 0
ODE=ode(fun)
ODE.set_integrator('dopri5')#, max_step=0.01)
ODE.set_solout(fout2)
ODE.set_initial_value(Y0, t0) 
ODE.integrate(t2)     
t2=ts[-1]
print()
Y=np.array(YY)

#до торможения
for i in range(len(ts)):
    ym.append(Rem*math.cos(Y[i][4]))
    xm.append(Rem*math.sin(Y[i][4]))

for i in range(len(ts)):
    ys.append(Y[i][1])
    xs.append(Y[i][0])

def circle(x,r):
    return np.array(np.sqrt(r**2-x**2))

plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots()

xc,yc=[],[]
for i in range(0, 630):
    xc.append(6371000*math.cos(i/100))
    yc.append(6371000*math.sin(i/100))
plt.plot(xc,yc,linewidth=1, c = 'green')

plt.axis('equal')
plt.plot(xs,ys, marker="*", c="#3589EF", markersize=0.5)
plt.plot(xm,ym, marker="*", c="#FF00FF", markersize=0.5, linestyle = '--')
plt.xlabel("X, м")
plt.ylabel("Y, м")
ax.set_title("Траектория от Земли до Луны")
plt.grid(False)
plt.show()

t3=1308 #Время до начала торможения
t4=357 #Время торможения
YY==[]
YYY=[]
arr=[]
#свободный полет до начала торможения
x,y,ux,uy,alphaM=list(Y[-1])       
Y0,t0=[x,y,ux,uy,alphaM],0
ODE.set_solout(fout3Stop)
ODE.set_initial_value(Y0, t0) 
ODE.integrate(t3)

#Торможение
x,y,ux,uy,alphaM=list(YY[-1])
Y0,t0=[x,y,ux,uy,alphaM],0
ON=-1
ODE.set_solout(fout3Stop)
ODE.set_initial_value(Y0, t0) 
ODE.integrate(t4)
mass += -q*(t4) # потраченное на торможение топливо 

#Притяжение ЛК во время свободного полета вокругу Луны для снижения высоты 
ON=0
t5=12915
ts=[]
x,y,ux,uy,alphaM=list(YY[-1])
Y0,t0=[x,y,ux,uy,alphaM],0
ODE.set_solout(fout4Stop)
ODE.set_initial_value(Y0, t0) 
ODE.integrate(t5)

#торможение до первой космической
x,y,ux,uy,alphaM=list(YY[-1])
Y0,t0=[x,y,ux,uy,alphaM],0
ON=-1
t6=2.01
ODE.set_solout(fout3Stop)
ODE.set_initial_value(Y0, t0) 
ODE.integrate(t6)
mass+=-q*t6

#Свободный полет по орбите Луны
ON=0
t7=8000
arr=[]
x,y,ux,uy,alphaM=list(YY[-1])
Y0,t0=[x,y,ux,uy,alphaM],0
ODE.set_solout(fout4Stop)
ODE.set_initial_value(Y0, t0) 
ODE.integrate(t7)
a=max(arr)
b=min(arr)

Y=np.array(YY)

Y1=np.array(YYY)    
xs1,ys1=list(Y1[:,0]),list(Y1[:,1])
xb1=np.linspace(-Rm,Rm,500)
yb1=circle(xb1,Rm)

#Траектория вблизи Луны
plt.style.use('seaborn-whitegrid')
fig, ax = plt.subplots()

plt.axis('equal')
plt.plot(xs1,ys1, c="#3589EF", linewidth=0.5)
plt.plot(xb1,yb1, c="orange", linewidth=1)
plt.plot(xb1,-yb1, c="orange", linewidth=1)
plt.xlabel("X, м")
plt.ylabel("Y, м")
ax.set_title("Траектория вблизи Луны")
plt.grid(False)
plt.show()


x_m=Rem*math.sin(Y[-1][4]) #проекция радиуса Луны в визуализируемую систему координат(ВСК)
y_m=Rem*math.cos(Y[-1][4])
x,y=Y[-1][0],Y[-1][1]
print(b-Rm,a-Rm)
V=np.sqrt((Y[-1][2]-w0*Rem*math.cos(Y[-1][4]))**2+(Y[-1][3]+w0*Rem*math.sin(Y[-1][4]))**2) #скорость в ВСК
print("Velosity:",V)
print("Orbital velosity:",np.sqrt(Nm*2/(b+a)))
t_spent=t0+t1+t2+t3+t4+t5+t6+t7+t_previous 
print("Spent Time :",t_spent)
H=np.sqrt(xs1[-1]**2+ys1[-1]**2)-Rm
print("Height: ",H)

#Переход в ВСК
X_m=x_m*math.cos(delta)+y_m*math.sin(delta)
Y_m=(-x_m*math.sin(delta)+y_m*math.cos(delta))*(-1)
X=x*math.cos(delta)+y*math.sin(delta)
Y=(-x*math.sin(delta)+y*math.cos(delta))*(-1)
print('Final coordinates:',X_m,Y_m,X,Y)
F_spent=(t4+t6)*q
print('sent Fuel:',F_spent)

W=(-V/(Rm+H)) #положительное направление вращения-направление вращения луны
out.write(str(X)+' '+str(Y)+' '+str(X_m)+' '+str(Y_m)+' '+str(W)+' '+str(F_spent)+' '+str(t_spent))
out.close()
