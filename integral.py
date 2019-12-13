import scipy.integrate
import time
import numpy as np
import math
import pylab as pl
import CDM 


coloumn0 = []
coloumn1 = []

with open(r"DAT.txt", "r+") as f:
    data = f.readlines()
   
    for line in data:
        coloumn1.append(line.strip().split("   ")[1])
        #strip to remove the \n
        #split at every interval of comma
        #second element is indexed 1
    for line in data:
        coloumn0.append(line.strip().split("   ")[0])
        #strip to remove the \n
        #split at every interval of comma
        #second element is indexed 1

data1 = np.zeros(len(coloumn0)) #redshift coletado
data2 = np.zeros(len(coloumn1)) #magnitude aparente coletada
i=0
while i < len(coloumn0):
    data1[i] = np.array([float(coloumn0[i])])
    i=i+1

i=0
while i < len(coloumn0):
    data2[i] = np.array([float(coloumn1[i])])
    i=i+1


Dc =np.zeros(300)
DA =np.zeros(300)
DL0 =np.zeros(300)
DL1 =np.zeros(300)
DL2 =np.zeros(300)
MI =np.zeros(300)
Z = np.linspace(7.9385996e-02,9.5739698e-01,300)



M = CDM.LCDModel(70.0,299792.4580, 0.0, 0.27, 0.72, 0.00, -1, 0.002)

'''i=0
ini = time.time()
while i < len(Z):

    Dc[i] = M.r(Z[i])
    i=i+1    
fim = time.time()
t = fim - ini
print('o tempo para calcular todos os pontos da distância comóvel ====>  ',t)
'''


'''i=0
ini = time.time()
while i < len(Z):

    DA[i] = M.DA(Z[i])
    i=i+1    
fim = time.time()
t = fim - ini
print('o tempo para calcular todos os pontos da distância Angular ====>  ',t)
'''

i=0
ini = time.time()
while i < len(Z):

    DL0[i] = M.DL(Z[i])
    print('DL0: ', DL0[i], '\n')
    i=i+1    
fim = time.time()
t = fim - ini
print('o tempo para calcular todos os pontos da distância luminosa ====>  ',t)


M = CDM.LCDModel(70.0,299792.4580, 0.0, 0.27, 0.72, 0.1, -1, 0.002)


i=0
ini = time.time()
while i < len(Z):

    DL1[i] = M.DL(Z[i])
    print('DL1: ', DL1[i], '\n')
    i=i+1    
fim = time.time()
t = fim - ini
print('o tempo para calcular todos os pontos da distância luminosa ====>  ',t)


M = CDM.LCDModel(70.0,299792.4580, 0.0, 0.27, 0.72, -0.1, -1, 0.002)


i=0
ini = time.time()
while i < len(Z):

    DL2[i] = M.DL(Z[i])
    print('DL2: ', DL2[i], '\n')
    i=i+1    
fim = time.time()
t = fim - ini
print('o tempo para calcular todos os pontos da distância luminosa ====>  ',t)

'''
i=0
ini = time.time()
while i < len(Z):

    MI[i] = M.MI(Z[i])
    i=i+1    
fim = time.time()
t = fim - ini
print('o tempo para calcular todos os pontos do modulo de distância  ====>  ',t)

M0 = -19.3
c = data2 - M0
Dif = c - MI

i = 0
while i < 300:
    print(f'dado {c[i]} e teorico {MI[i]}')
    print(f'\n onde a diferença: {Dif}')
    i=i+1
'''

#pl.plot(Z, Dc, 'b.', label='DC')
#pl.plot(Z, DA, 'r.', label='DA')
#pl.plot(Z, DL0, 'y.', label='DL0')
#pl.plot(Z, DL1, 'y.', label='DL1')
#pl.plot(Z, DL2, 'y.', label='DL2')
pl.plot(Z, MI, 'g.', label='MI')
pl.plot(data1, c, 'r.', label='MI')

pl.legend(loc='upper left')
pl.show()

