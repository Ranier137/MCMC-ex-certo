import numpy as np

f1 = open('lcparam_DS17f.txt', 'r')
bytes = f1.readlines()[1:]



error = []

for line in bytes:
    error.append(line.strip().split(" ")[5])

errornp = np.array(error)
dim = len(errornp)

C = np.zeros([dim,dim])

i=0
while i < dim:
    C[i,i] = errornp[i]
    i = i +1
print('C\n ', C)

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
print('Csys\n ',Csys)

#print('linha', Csys[0])

Ctot = C + Csys

print('Ctot\n ',Ctot)






 
