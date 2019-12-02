import numpy as np
from matplotlib import pyplot as plt 
import matplotlib.cm as cm 

f1 = open('lcparam_DS17f.txt', 'r')
bytes = f1.readlines()[1:]


#print(bytes)

error = []
redc = []
redl = []

for line in bytes:
    #print(line,'\n\n')
    error.append(line.strip().split(" ")[5])
    redc.append(line.strip().split(" ")[1])

i = 0
while i < len(redc):
    redl.append(redc[39-i])
    i=i+1
#print(redc)
#print(redl)
errornp = np.array(error)
dim = len(errornp)

C = np.zeros([dim,dim])

i=0
while i < dim:
    C[i,i] = errornp[i]
    i = i +1
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

print('Ctot\n ',Ctot)


'''

x = np.array(redc)    
y = np.array(redl) 
X, Y = np.meshgrid(x, y) 

Z = X + Y 

plt.pcolormesh(X, Y, Z, cmap = cm.gray) 
plt.show()
'''

n = 256 
x = np.array([float(i) for i in redc]) 
y = np.array([float(i) for i in redc])

X, Y = np.meshgrid(x, y) 

#print(np.meshgrid(x,y))

Z = Ctot  

plt.xlim(0.0, 1.6)
plt.ylim(1.6, 0.0)
plt.pcolormesh(X, Y, Z, cmap = cm.gray)
plt.colorbar(label='scale') 
plt.clim(0.001, -0.001)
plt.show()
 
