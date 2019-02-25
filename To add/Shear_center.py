import math
import numpy as np
import random
#t_1=t_sk
#t_2=t_sp
#spar_height=h

#S_y/I_xzz
I_xx=1565571000 #mm^4

spar_height=205
large_pitch=81.33511390749766
small_pitch=8.512094985784858
stiffener_pitch=89.84720889328251
stiffener_area=40.56
t_1=1.1
t_2=2.8

boomarea=[]
    
angle=math.degrees(0.2012196329874848)
#print (angle)
z_pos= [560.0763955533588, 560.0763955533588, 470.22918666007627, 470.22918666007627, 380.38197776679374, 380.38197776679374, 290.5347688735112, 290.5347688735112, 200.68755998022868, 200.68755998022868, 110.84035108694617, 110.84035108694617, 102.5, 102.5, 36.920485267356455, 36.920485267356455, 0]
#print(z_pos)
y_pos= [9.163521305036266, -9.163521305036266, 27.49056391510882, -27.49056391510882, 45.817606525181375, -45.817606525181375, 64.14464913525393, -64.14464913525393, 82.47169174532648, -82.47169174532648, 100.79873435539903, -100.79873435539903, 102.5, -102.5, 78.7754863357313, -78.7754863357313, 0]

#print(y_pos)

boomarea2 = []
#boomarea 1 and 2
b1 = (t_1*stiffener_pitch)/6*(2+y_pos[1]/y_pos[0]) + stiffener_area + (t_1*stiffener_pitch)/6*(2+y_pos[2]/y_pos[0])
boomarea2.append(b1)
boomarea2.append(b1)
#boomarea 3 to 10
for i in range(0,4):
    j = 2+2*i
    b2 = ((t_1*stiffener_pitch)/6)*(2+y_pos[j-2]/y_pos[j]) + stiffener_area + ((t_1*stiffener_pitch)/6)*(2+y_pos[j+2]/y_pos[j])
    boomarea2.append(b2)
    boomarea2.append(b2)  
#boomarea 11 and 12
b11=(t_1*small_pitch)/6*(2+y_pos[12]/y_pos[10]) + stiffener_area + (t_1*stiffener_pitch)/6*(2+y_pos[8]/y_pos[10])
boomarea2.append(b11)
boomarea2.append(b11)
#boomarea 13 and 14
b13=(t_2*spar_height)/6*(2+y_pos[13]/y_pos[12]) + (t_1*small_pitch)/6*(2+y_pos[10]/y_pos[12]) + (t_1*large_pitch)/6*(2+y_pos[14]/y_pos[12])
boomarea2.append(b13)
boomarea2.append(b13)
#boomarea 15 and 16
b15=(t_1*large_pitch)/6*(2+y_pos[12]/y_pos[14]) + stiffener_area + (t_1*stiffener_pitch)/6*(2+y_pos[16]/y_pos[14])
boomarea2.append(b15)
boomarea2.append(b15)
#boomarea 17
b17=((t_1*stiffener_pitch)/6)*2 + stiffener_area
boomarea2.append(b17)
boomarea=boomarea2
q_bottom=[0]
q_prev = float(0)
for i in range(0,6):
    j = 11-(2*i)
    q=(-boomarea[j]/I_xx)*y_pos[j] + q_prev
    q_prev = q
    q_bottom.append(q)
#print(q_bottom)
#print (q_bottom[-1])

q_top=[]
q_prev = q_bottom[-1]
for i in range(0,5):
    j = (2*i)
    q=(-boomarea[j]/I_xx)*y_pos[j] + q_prev
    q_prev = q
    q_top.append(q)
q_top.append(0)
#print(q_top)

q_spar=[]
q=(-boomarea[13]/I_xx)*y_pos[13]
q_spar.append(q)
#print(q_spar)

q_curve=[0]

q16=(-boomarea[15]/I_xx)*y_pos[15]
q_curve.append(q16)
q17=(-boomarea[16]/I_xx)*y_pos[16]+q16
q_curve.append(q17)
#q15=(-boomarea[14]/I_xx)*boomy[14]+q17
#q_curve.append(q15)
q_curve.append(0)
# print(q_curve)

q_right=q_bottom+q_top

a_1=((q_curve[1]+q_curve[2])*stiffener_pitch)/t_1
b_1=-(q_spar[0])*spar_height/t_2

m_1=((2*stiffener_pitch+2*large_pitch)/t_1)+(spar_height/t_2)
m_2=(-spar_height/t_2)

a_2=(sum(q_right)*stiffener_pitch)/t_1
b_2=(q_spar[0])*spar_height/t_2

n_1=(-spar_height/t_2)
n_2=((11*stiffener_pitch+2*small_pitch)/t_1)+(spar_height/t_2)
#print(m_1,m_2,a_1,b_1,n_1,n_2,a_2,b_2)
A=np.array([[m_1,m_2],[n_1,n_2]])
B=np.array([[-a_1-b_1],[-a_2-b_2]])
#print (np.linalg.solve(A, B ))
q_s0=np.linalg.solve(A,B)
q_s01=q_s0[0][0]
q_s02=q_s0[1][0]
d=spar_height*math.sin(math.radians(90-angle))
#M around 5 by upper part
M_1=(((-q_s02*5+sum(q_top[0:5]))*stiffener_pitch)+(-q_s02*small_pitch))*d
#M caused by q between 1 and 2
M_2=(q_bottom[-1]-q_s02)*(y_pos[0]-y_pos[1])*(z_pos[0]-z_pos[14])
#M caused by curve
M_3= -((q_curve[1]+q_s01)*(2*(y_pos[14]))*(z_pos[13]-z_pos[14]))
M_4=-((q_s01*(z_pos[13]-z_pos[14])*(y_pos[14]-y_pos[13]))+(q_s01*(y_pos[12]-y_pos[14])*(z_pos[13]-z_pos[14])))

M=M_1+M_2+M_3+M_4
shearcenter_u=M+(spar_height/2)
print(shearcenter_u)
