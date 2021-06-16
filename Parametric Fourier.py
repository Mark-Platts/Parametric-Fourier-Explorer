#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
from numpy import trapz
import math
import numpy as np
pi = math.pi

'''
READ ME
This script uses Fourier series equations as parametric equations to plot 2D shapes made out of straight lines. 
To draw the shape, fill the 'points' variable below with 2D coordinates for the ends of the lines in order of the parametrisation
Examples are provided below.
The 'order' variable sets the accuracy of the five outcomes. The number represents the order to which the Fourier series sums are evaluated.
A higher number will give better accuracy.
'''

points = [[1,1],[0,0],[1,-1],[-1,-1],[0,0],[-1,1]]
orders = [2,4,6,10,100]

#EXAMPLES

#triangle
#points = [[0,0],[0,1],[1,0]]

#square
#points = [[1,1],[1,-1],[-1,-1],[-1,1]]

#hourglass
#points = [[1,1],[-1,-1],[1,-1],[-1,1]]

#proper hourglass
#points = [[1,1],[0,0],[1,-1],[-1,-1],[0,0],[-1,1]]

#generic
#points = [[0,0],[1,1],[2,1],[3,-2],[0,-3],[-3,2]]

#Encircled Pentagon
'''
points = [[0,1],
          [math.cos((pi/2)+2*pi/5),math.sin((pi/2)+2*pi/5)],
          [math.cos((pi/2)+2*2*pi/5),math.sin((pi/2)+2*2*pi/5)],
          [math.cos((pi/2)+3*2*pi/5),math.sin((pi/2)+3*2*pi/5)],
          [math.cos((pi/2)+4*2*pi/5),math.sin((pi/2)+4*2*pi/5)],
          [0,1],
          [math.cos((pi/2)+2*2*pi/5),math.sin((pi/2)+2*2*pi/5)],
          [math.cos((pi/2)+4*2*pi/5),math.sin((pi/2)+4*2*pi/5)],
          [math.cos((pi/2)+2*pi/5),math.sin((pi/2)+2*pi/5)],
          [math.cos((pi/2)+3*2*pi/5),math.sin((pi/2)+3*2*pi/5)],
         ]
'''

#Customisations for the outcome plots
view_draw_dots = 1
view_draw_grid = 1
view_four_grid = 1


#LOGIC

#Gets highest order in 'orders'
hig_ord = 0
for i in orders:
    if i > hig_ord:
        hig_ord = i
        
#Gets the needed mathematical data for the lines that make the shape
lines = []
for i in range(len(points)-1):
    a = []
    a.append(points[i])
    a.append(points[i+1])
    lines.append(a)
a = []
a.append(points[-1])
a.append(points[0])
lines.append(a)

#Finds the period of the Fourier series
p = len(lines)

#Given initial coordinate I, final coordinate F and line number L, returns contribution to a0 coefficient
def a0_segment(I,F,L):
    return((F+I)/2)

#Given initial coordinate I, final coordinate F,line number L,and n, returns contribution to an coefficient
def an_segment(I,F,L,n):
    N = 2*pi*n/p
    D = F-I
    term1 = ((I-L*D)/N)*(math.sin(N*(L+1))-math.sin(N*L))
    term2 = (D/N**2)*(N*(L+1)*math.sin(N*(L+1))-N*L*math.sin(N*L)+math.cos(N*(L+1))-math.cos(N*L))
    return(term1+term2)

#Given initial coordinate I, final coordinate F,line number L,and n, returns contribution to bn coefficient
def bn_segment(I,F,L,n):
    N = 2*pi*n/p
    D = F-I
    term1 = (-(I-L*D)/N)*(math.cos(N*(L+1))-math.cos(N*L))
    term2 = (D/N**2)*(-N*(L+1)*math.cos(N*(L+1))+N*L*math.cos(N*L)+math.sin(N*(L+1))-math.sin(N*L))
    return(term1+term2)

#c defines which coordinate to get coefficient for: 0=x, 1=y
def get_a0(c):
    a = 0
    for i in range(len(lines)):
        a += a0_segment(lines[i][0][c],lines[i][1][c],i)
        a0 = a/p
    return(a0)

def get_an(c,n):
    a = 0
    for i in range(len(lines)):
        a += an_segment(lines[i][0][c],lines[i][1][c],i,n)
        an = 2*a/p
    return(an)

def get_bn(c,n):
    b = 0
    for i in range(len(lines)):
        b += bn_segment(lines[i][0][c],lines[i][1][c],i,n)
        bn = 2*b/p
    return(bn)



def sine_part(n,T,x):
    return(math.sin(2*pi*n*x/T))

def cosine_part(n,T,x):
    return(math.cos(2*pi*n*x/T))


#Puts together the final Fourier series and returns y for given x
def end_fun(a0,a,b,T,x,o):
    y = a0
    for i in range(o):
        y += a[i]*cosine_part(i+1,T,x)
        y += b[i]*sine_part(i+1,T,x)
    return(y)


#Finds the coefficients for both coordinates
ax0 = get_a0(0)
ax = []
bx = []
ay0 = get_a0(1)
ay = []
by = []
for i in range(hig_ord):
    n = i+1
    ax.append(get_an(0,n))
    bx.append(get_bn(0,n))
    ay.append(get_an(1,n))
    by.append(get_bn(1,n))


#Prepare storage for plot data
fx0=[]
fy0=[]
fx1=[]
fy1=[]
fx2=[]
fy2=[]
fx3=[]
fy3=[]
fx4=[]
fy4=[]

t = np.arange(0, p, p/1000)

#Get plot data
for i in t:
    fx0.append(end_fun(ax0,ax,bx,p,i,orders[0]))
    fy0.append(end_fun(ay0,ay,by,p,i,orders[0]))
    fx1.append(end_fun(ax0,ax,bx,p,i,orders[1]))
    fy1.append(end_fun(ay0,ay,by,p,i,orders[1]))
    fx2.append(end_fun(ax0,ax,bx,p,i,orders[2]))
    fy2.append(end_fun(ay0,ay,by,p,i,orders[2]))
    fx3.append(end_fun(ax0,ax,bx,p,i,orders[3]))
    fy3.append(end_fun(ay0,ay,by,p,i,orders[3]))
    fx4.append(end_fun(ax0,ax,bx,p,i,orders[4]))
    fy4.append(end_fun(ay0,ay,by,p,i,orders[4]))

#Restates coordinates of each point so that they can be graphed
lx = []
ly = []
for i in points:
    lx.append(i[0])
    ly.append(i[1])
lx.append(points[0][0])
ly.append(points[0][1])

#Make plots and show
fig, axs = plt.subplots(3, 2)
axs[0, 0].plot(lx, ly, 'r-')
if len(points) < 30 and view_draw_dots != 0:
    axs[0, 0].plot(lx, ly, 'ro')
if view_draw_grid != 0:
    axs[0, 0].grid()
axs[0, 0].set_title('Original Points')
axs[0, 0].axis('equal')
axs[0, 1].plot(fx0, fy0, 'b-')
axs[0, 1].set_title('n = '+str(orders[0]))
axs[0, 1].axis('equal')
axs[1, 0].plot(fx1, fy1, 'b-')
axs[1, 0].set_title('n = '+str(orders[1]))
axs[1, 0].axis('equal')
axs[1, 1].plot(fx2, fy2, 'b-')
axs[1, 1].set_title('n = '+str(orders[2]))
axs[1, 1].axis('equal')
axs[2, 0].plot(fx3, fy3, 'b-')
axs[2, 0].set_title('n = '+str(orders[3]))
axs[2, 0].axis('equal')
axs[2, 1].plot(fx4, fy4, 'b-')
axs[2, 1].set_title('n = '+str(orders[4]))
axs[2, 1].axis('equal')
if view_four_grid != 0:
    axs[0, 1].grid()
    axs[1, 0].grid()
    axs[1, 1].grid()
    axs[2, 0].grid()
    axs[2, 1].grid()
fig.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.5,
                    wspace=0.35)
fig.show()


# In[ ]:




