import numpy as np
from matplotlib import pyplot as plt 
import matplotlib.cm as cm 
from numpy.linalg import inv
import time

f1 = open('lcparam_DS17f.txt', 'r')
l = f1.readlines()[1:]


#print(l)

error = []
redc = []
redl = []
mb = []

for line in l:
    #print(line,'\n\n')
    error.append(line.strip().split(" ")[5])
    redc.append(line.strip().split(" ")[1])
    mb.append(line.strip().split(" ")[4])



print('aparente \n', mb)
i = 0
while i < len(redc):
    redl.append(redc[39-i])
    i=i+1
#print(redc)
#print(redl)
errornp = np.array(error)
dim = len(errornp)

Ce = np.zeros([dim,dim])

i=0
while i < dim:
    Ce[i,i] = errornp[i]
    i = i +1

C = Ce**2
#print('C\n ', C)

f1.close()

f2 = open('syserror.txt', 'r')

Csys = np.zeros([40,40])
r = f2.readlines()
errorsys = []

for line in r:
    errorsys.append(line.strip().split(" ")[0])


errornpsys = np.array(errorsys)
dimsys = len(errornpsys)


i=0
j=0
while i < 40:
    while j < 40:
        Csys[i,j] = errornpsys[j + i*40]
        j= j+1
    j=0
    i = i+1     

np.set_printoptions(threshold=1600)
#print('Csys\n ',Csys)

#print('linha', Csys[0])

Ctot = C + Csys


#print('Ctot\n ',Ctot)


Corr = np.zeros([dim,dim])
i=0
j=0
while i < 40:
    while j< 40:
        Corr[i,j] = Ctot[i,j]/((Ctot[i,i]*Ctot[j,j])**(0.5)) 
        j=j+1
    j=0
    i=i+1

#print('\n\nCorr\n', Corr)

inicio = time.time()
Ctotinv = inv(Ctot)
fimzao = time.time()

tempao = fimzao - inicio

print('Inversa\n',Ctotinv)

print('\n\n\n', tempao, '\n\n\n')

cross = np.matmul(Ctot,Ctotinv)
print('\n', cross)
'''
x = np.array(redc)    
y = np.array(redl) 
X, Y = np.meshgrid(x, y) 

Z = X + Y 

plt.pcolormesh(X, Y, Z, cmap = cm.gray) 
plt.show()
'''
'''
n = 256 
x = np.array([float(i) for i in redc]) 
y = np.array([float(i) for i in redc])

X, Y = np.meshgrid(x, y) 

#print(np.meshgrid(x,y))

Z = Ctot  
W = Corr  



plt.xlim(0.014, 1.6123)
plt.ylim(1.6123, 0.014)
plt.pcolormesh(X, Y, Z, cmap = cm.gray)
plt.colorbar(label='scale') 
plt.clim(0.0015, -0.0015)
plt.show()
plt.xlim(0.014, 1.6123)
plt.ylim(1.6123, 0.014)
plt.pcolormesh(X, Y, W, cmap = cm.gray)
plt.colorbar(label='scale') 
plt.clim(1, -1)
plt.show()
 '''
