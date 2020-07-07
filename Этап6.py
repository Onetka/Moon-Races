import numpy as np
import matplotlib.pyplot as plt
import math as m
import random

from scipy.integrate import ode

arr = np.loadtxt("Этап 5 выход.txt", delimiter='\n', dtype=np.float)
h = arr[0] #высота над землей
vx0 = arr[1] #сорбитальная скорость
vy0 = 0
#print(h,vx0,vy0)

max_overload = 0

totalTotalTime = arr[2]
totalTime = 0
overload = []
g_e = 9.81
R_e = 6375*10**3 
S = 0.25*m.pi*3.9**2  
p0 = 1.584 
mu = 398600e9
Cx = 0.85
Cy0 = 0.134*0.85
vx0 = 10000
TotalMass = 5500
F = True
qwerty = 20000
result = open('Этап 6 выход.txt', 'w')
result.close()

def Safety(F1, F2, F):
    np = (F1**2 + F2**2)**0.5/g_e  
    if np>=10:
        print ("Пилот разбился ;C ", np)
        print (F1, F2)
        F = False
    overload.append(np)
    return F


def fout1(t, y):
        ts.append(t)
        ys.append(list(y.copy()))
        y1, y2, y3, y4 = y
        if m.sqrt(y2*y2+y4*y4)<=300 and (m.sqrt(y1*y1 +y3*y3) - R_e)<10000:
          return -1

def density(x,y): #по барометрической формуле рассчитываем плотность воздуха 
    h = m.sqrt(x*x +y*y) - R_e
    g = g_e*m.sqrt(R_e/(R_e+h))
    p = 1.225*m.exp(-0.029*g*h/(8.31*273))
    return p

def tang_angle(h):
    global  qwerty
    angle_lst = [0]
    Cy_lst = np.zeros(qwerty)
    Cy_lst[0] = 0.134*0.85*m.cos(0*m.pi/180)
    for i in range(1,qwerty):
        if (h-R_e)<70000:
                angle_lst.append(angle_lst[i-1] + 1)
        else:
            angle_lst.append(0)
        Cy_lst[i] = Cy_lst[i-1]*m.cos(angle_lst[i]*m.pi/180)
            
    angle0 = angle_lst [-1]
    Cy_lst0 = Cy_lst[-1]
    return angle0*m.pi/180, Cy_lst0
        
# функция правых частей системы ОДУ
def function1(t, y):
         global TotalMass, max_overload, F, angle, Cy    
         if F:
             y1, y2, y3, y4 = y
             result = open('Этап 6 выход.txt', 'a')         
             vv=y2*y2+y4*y4
             v = m.sqrt(vv)
             r = (y1**2+y3**2)**0.5
             angle = tang_angle(r)[0]
             Cy =tang_angle(r)[1]
             p = density(y1,y3)
             Q = (Cx*S*p*(v**2))/(2*TotalMass)
             N = (Cy*S*p*(v**2))/(2*TotalMass)
             Py = N*m.cos(angle)
             F = Safety(Q, N, F)
             ax = -Q*(y2/v) - Py*(y4/v) - (mu/r**2)*y1/r
             ay = -Q*(y4/v) + Py*(y2/v) - (mu/r**2)*y3/r 
             result.write(str(r) + '\t' + str(m.sqrt(vv)) + '\t' +
                          str(m.sqrt(ax**2 +ay**2)) + '\t' + str(t) + 
                          str(round(angle*180/m.pi)) + '\t'+ str(round(Cy)) + 
                          '\n')
             result.close()
             return [y2, ax, y4,ay] 
     
x_start = 0
Vx_start = vx0
y_start = h + R_e
Vy_start = -vy0

xc,yc=[],[]
for i in range(0, 630):
    xc.append(0.001*R_e*m.cos(i/100))
    yc.append(0.001*R_e*m.sin(i/100))
tmax = 15000

y0,t0=[x_start,  Vx_start, y_start, Vy_start], 0 # начальные условия 

ODE=ode(function1)
ODE.set_integrator('dopri5')
ODE.set_solout(fout1)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) # задание начальных значений  
ODE.integrate(tmax)      # решение ОДУ
Y=np.array(ys)
xgraph = Y[:,0]/1000
ygraph = Y[:,2]/1000


print()
print('Конечная высота ',"%.3f" % (m.sqrt(Y[-1:,0]**2 + Y[-1:,2]**2) - R_e),'м')
print('Конечная скорость ', "%.3f" % (m.sqrt(Y[-1:,1]*Y[-1:,1] + Y[-1:,3]*Y[-1:,3])),'м/с')
totalTime+=ts[-1]
totalTotalTime +=totalTime
print('Полное время снижения', "%.0f" %totalTime,'с')
max_overload = max(overload)
print ('Макс перегрузка ', round(max_overload), 'g')
print ('Полное время гонки ', round(totalTotalTime/(60*60*24)),' дн')

end = open('Финал.txt', 'w')
end.write('Конечная высота '+ str(round(m.sqrt(Y[-1:,0]**2 + Y[-1:,2]**2) - R_e)) +' м \n' +
          'Конечная скорость '+ str(round(m.sqrt(Y[-1:,1]*Y[-1:,1] + Y[-1:,3]*Y[-1:,3]))) +' м/с \n'+
          'Полное время снижения '+ str(round(totalTime)) +' с \n' +
          'Макс перегрузка ' + str(round(max_overload,1)) + ' g \n' +
          'Полное время гонки ' + str(round(totalTotalTime/(60*60*24))) + ' дн')
end.close()


plt.style.use('seaborn-white')
plt.title("Тормозим \n ", size=35)       
    
plt.plot(x, y, marker ="*", c="white", linestyle=" ")
plt.plot(xgraph, ygraph, linewidth=2, color='b')
plt.axis('equal')
#plt.text (1000, 5000, u'the Earth', size=40, color="bisque")
plt.plot(xc, yc, linewidth=3, color='green')
#plt.xlim(0, 3500)
#plt.ylim(5000, 7000)
plt.xlabel('X, км')
plt.ylabel('Y, км')
plt.grid(False)
plt.show()
plt.plot(xgraph, ygraph, linewidth=2, color='b')
plt.axis('equal')
#plt.text (1000, 5000, u'the Earth', size=40, color="bisque")
plt.plot(xc, yc, linewidth=3, color='green')
plt.xlim(0, 2000)
plt.ylim(6000, 7000)
plt.xlabel('X, км')
plt.ylabel('Y, км')
plt.grid(False)
plt.show()

def visual(data):
    for i in range(0,len(data)):
        y = data[i]
        plt.plot(i, y, c='b',marker = ".", markersize=5, linestyle="--")
        plt.pause(0.001)
    

dataa = np.zeros(len(Y[:,2]))
