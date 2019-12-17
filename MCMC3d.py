import math
import numpy as np
import CDM 
import matplotlib.pyplot as plt
import time
from numpy.linalg import inv

#o problema e gerar uma distribuiçao no espaço dos parametros OM_m(materia) e w(dark energy state equation) dado o modelo LCDM (usando o pacote que escrevi CDM). Os parametros sao denotados por x = [OM_m, w].

#aqui defino a distribuiçao de transiçao. Pega os dois parametros, guardados no vetor x, e devolve outros dois parametros de acordo com a distribuiçao normal, centrada em cada um dos parametros iniciais:

#os passos de cada parâmetro
pass0 = 0.05
pass1 = 0.5
pass2 = 0.05



Transition = lambda x: [np.random.normal(x[0], pass0), np.random.normal(x[1], pass1), np.random.normal(x[2], pass2)]

#Definir a distribuiçao prior para cada calor dos parametros armazenado no vetor x
    
def prior(x):
# coloquei aqui a informaçao de que a probabilidade deve se anular caso algum dos parametros se anule ou se torne negativo.
    if x[0]<= 0.: 
        return 0.0
    if x[1] > -1./3.:
        return 0.0
    if  abs(x[2]) < 18.7 or abs(x[2]) > 19.8:
        return 0.0 
#usei prior uniforme. coloquei apenas 1.0 pq nao importa a constante de normalizaçao que deve aparecer, nao influenciara no metodo.	
    else:
        return 1.0

# aqui devo defini o logaritmo da likelihood*prior
#alem dos parmetros x, as entradas sao sig = 0.4 desvio padrao (assumindo covariancia diagonal) e dadaos representados por data1 = coluna dos redshifts e data2 = coluna das respectivas magnitude aparente.


def LnLike(x,data1, data2, Cinv): 
    d = len(data1) 		#número de dados coletados
    M = CDM.LCDModel(72., 299792.4580, 0.0, x[0], 1.0-x[0], 0.0, x[1], 0.0005) #é criado uma instância do modelo LCDM, fixando constante de hubble atualmente H0 = 72, c = 299792.458km/s, OM_r = 0.
    i=0
    M0 = x[2]	#assumindo que todas as SNIa possuem magnitude absoluta iguais
    Mbo = np.zeros(len(data2))
    deltaMb = np.zeros(len(data2))	
    while i < len(data2):
        Mbo[i] = data2[i] #magnitude aparente observado
        deltaMb[i] = Mbo[i] - M.Mb(data1[i], M0) #diferença entre o magnitude aparentes, observado e o teórico (dado os possíveis valores dos parâmbetros).        
        i = i+1                         
    Qui2 = np.matmul(np.matmul(deltaMb, Cinv),deltaMb)    
    result = np.log(prior(x))-Qui2/2. # logaritmo da posterior ln(post) = ln(L*prior) = ln(L) + ln(prior)
    print(f'Ln(post) =  {result} \n')
    return result

# comparar as probabilidades posterior do ponto atual e possível próximo ponto no espaço de parâmetros
def Passo(xi,xp, data1, data2, Cinv, LnLike, Lnmax):
    print('este é o xi atualmente: ',xi, '\n') 
    LNi = LnLike(xi, data1, data2, Cinv) #LnLike do ponto inicial/atual xi
    print('este é o xp sorteado: ',xp, '\n') 
    LNp = LnLike(xp, data1, data2, Cinv) #LnLike do possível ponto posterior xp
    alpha = np.random.uniform(0.,1.) #sorteio uniforme de um número entre 0 e 1.
    D = LNp-LNi
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> diferença ', D) 
    if LNp > LNi: #caso a probabilidade favorece o próximo ponto xp.
        print('Lnp > Lni \n')
        Lnmax = LNp
        return 1 #retorna 1 caso aceito o próximo ponto xp
    else: #caso a probabilidade não favoreça xp, é necessário utilizar o sorteio de alpha.
        r = np.exp(D) #razão entre posterio de xp e posterior de xi
        print(f'Essa é a razão: {r}            Este o alpha: {alpha}\n') 
        if alpha < r: #critério de aceite do ponto xp
            print('alpha < r \n')
            return 1 # aceito
        else:
            print('alpha > r \n')
            return 0 #recusado


    



#IMPLEMENTAÇÃO DO MÉTODO MCMC



#______________________________________________________________________________________________________________________________________________________

#abrir arquivo com os dados
f1 = open('lcparam_DS17f.txt', 'r')
l = f1.readlines()[1:] # cada linha do arquivo é l
error = [] #armezenar erro diagonal
red = [] # armazenar redshift do primeiro até o ultimo
mb = [] #armazenar as magnitude aparente

for line in l:
    #print(line,'\n\n')
    error.append(line.strip().split(" ")[5])
    red.append(line.strip().split(" ")[1])
    mb.append(line.strip().split(" ")[4])



#print('aparente \n', mb)
print('redshift', red)
errornp = np.array(error)

Ce = np.zeros([len(errornp),len(errornp)])

i=0
while i < len(errornp):
    Ce[i,i] = errornp[i]
    i = i +1

C = Ce**2
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

#np.set_printoptions(threshold=1600)
#print('Csys\n ',Csys)

#print('linha', Csys[0])

Ctot = C + Csys


print('Ctot\n ',Ctot)


inicio = time.time()
Ctotinv = inv(Ctot)
fimzao = time.time()

tempao = fimzao - inicio

print('Inversa\n',Ctotinv)

print('\n\n\n', tempao, '\n\n\n')

cross = np.matmul(Ctot,Ctotinv)
print('\n', cross)

#________________________________________________________________________________________________________________________________________________
# Criando data1 e data2, ainda como listas vazias.

#Transformando as listas em np.arrays
data1 = np.zeros(len(red)) #redshift coletado
data2 = np.zeros(len(mb)) #magnitude aparente coletada
print(len(mb), len(red))
i=0
while i < len(red):
    data1[i] = np.array([float(red[i])])
    i=i+1

i=0
while i < len(mb):
    data2[i] = np.array([float(mb[i])])
    i=i+1


#chute inicial dos parâmetros para a cadeia
param_init0 = float(input(print('digite um chute inicial para Omega_m: \n'))) 
param_init1 = float(input(print('digite um chute inicial para w: \n')))
param_init2 = float(input(print('digite um chute inicial para M0: \n')))
iterations = int(input(print('digite o numero de iterações: \n')))    

sig2 = Ctotinv




xi = np.array([param_init0, param_init1, param_init2])    #ponto inicial da cadeia
#lista com os valores aceitos dos parâmetros. Inicialmente há apenas um valor
chain_m = [xi[0]]	
chain_w = [xi[1]]
chain_M0 = [xi[2]]
#lista de valores sorteados para os parâmetros porém rejeitados 
chain_mrej = []
chain_wrej = []
chain_M0rej = []
accepted = [0] #numeros de sorteios aceitados
rejected = [0] #numeros de sorteios rejeitados

ini = time.time()
for i in range(iterations): #iteração i
    print('_______________________________________________________________________________Faltam Iteração ', iterations - i, '\n') #monitoramento da contagem   
    Lnmax = LnLike(xi,data1, data2, sig2)
    xp = Transition(xi) #sorteio de um novo ponto, xp, a partir do ponto atual xi
    ret = Passo(xi,xp, data1, data2, sig2, LnLike, Lnmax) # verificar a condição de aceite de xp. Note que assumirá valores 1 (aceite) ou 0 (recusa). 
    print('retrun do passo: ', ret, '\n')
    if ret == 1: #se for aceito
        xi = xp #atualiza o ponto atual para o ponto sorteiado xp.
        accepted[0] = accepted[0] + 1 # atualiza quantos foram aceitos. 
    else: #não aceito
        rejected[0] = rejected[0] + 1 # atualiza quantos foram rejeitados. continua com o ponto xi.
	
	# insere na lista dos recusados xp
        chain_mrej.append(xp[0]) 
        chain_wrej.append(xp[1])
        chain_M0rej.append(xp[2])

    print(' true: ', accepted, '\n\n')
    print(' false: ', rejected, '\n\n')

#acrescenta os valores atualizados à cadeia.
    chain_m.append(xi[0])
    chain_w.append(xi[1])
    chain_M0.append(xi[2])
    
#após todas as iterações cálculo do valor médio    


dimension = len(chain_m)


C_m = np.array(chain_m)
C_w = np.array(chain_w)
C_M0 = np.array(chain_M0)
burnin = 0.8
burn = int(burnin*dimension)
burnper = burnin*100.
del chain_m[0:burn]
del chain_w[0:burn]
del chain_M0[0:burn]
n = dimension - burn
ValorMM = sum(chain_m)/n
ValorMw = sum(chain_w)/n
ValorMM0 = sum(chain_M0)/n

#montando a matrix covariancia dos parâmetros

cov = np.zeros([3,3])
a = cov[0,0] = sum((np.array(chain_m)-ValorMM)**2)/(n-1)
b = cov[0,1]= cov[1,0] = sum((np.array(chain_m)-ValorMM)*(np.array(chain_w)-ValorMw))/(n-1)
c = cov[1,1] = sum((np.array(chain_w)-ValorMw)**2)/(n-1)

d = cov[0,2] = cov[2,0] = sum((np.array(chain_m)-ValorMM)*(np.array(chain_M0)-ValorMM0))/(n-1)
e = cov[1,2] = cov[1,2] = sum((np.array(chain_w)-ValorMw)*(np.array(chain_M0)-ValorMM0))/(n-1)

f = cov[2,2] = sum((np.array(chain_M0)-ValorMM0)**2)/(n-1)



#definindo desvio padrão
sigx = np.sqrt(cov[0,0])
sigy = np.sqrt(cov[1,1])

#Os eixos da elipse de confianca
Lx = (a+c)/2.0 + np.sqrt(((a-c)/2)**2 + b**2)
Ly = (a+c)/2.0 - np.sqrt(((a-c)/2)**2 + b**2)
Theta = np.arctan((Lx-a)/b) 
print('esse é o angulo', Theta, 'em radianos \n\n' )

#95% de confianca
s1 = -2*np.log(1.0-0.95)
rx1 = np.sqrt(s1*Lx)
ry1 = np.sqrt(s1*Ly)

#68% de confianca
s2 = -2*np.log(1.0-0.68)
rx2 = np.sqrt(s2*Lx)
ry2 = np.sqrt(s2*Ly)




print(len(C_m), '\n\n\n')
print(len(chain_m))


fim = time.time()
tempoT = (fim - ini)/60. 
print(f'\n o valor esperado dos parametros sao: \n OM_m = {ValorMM} +/- {(a)**(0.5)} \n w = {ValorMw} +/- {(c)**(0.5)} \n M0 = {ValorMM0} +/- {(f)**(0.5)}\n\n Burnin = {burnper}% \n O tempo de execução total das iterações {tempoT} min \n os parmatros iniciais foram [{param_init0}, {param_init1}, {param_init2}] \n Ln máximo {Lnmax} \n Passos [{pass0}, {pass1}, {pass2}]')

#plotar os histogramas
plt.hist(C_m, bins = 100, label= 'Matter', color= 'blue')
plt.hist(C_w, bins = 100, label= 'w', color= 'red')
plt.hist(C_M0, bins = 100, label= 'M0', color= 'green')
plt.legend()
plt.show()



#pontos percorridos no espaço de parâmetros


plt.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca()
plt.scatter(chain_mrej, chain_wrej, marker='x',color='y',linewidths=1) # recusados em amarelo (yellow)
plt.scatter(C_m, C_w,marker='.',color='b',linewidths=3)# aceitos em azul (blue)
plt.scatter(ValorMM, ValorMw,marker='.',color='r',linewidths=3)# valor medio em vermelho 
phi = np.linspace(0.0, 2*np.pi, 100)
#elipse 95
x1 = rx1*np.cos(Theta)*np.cos(phi) - ry1*np.sin(Theta)*np.sin(phi) + ValorMM
y1 = rx1*np.sin(Theta)*np.cos(phi) + ry1*np.cos(Theta)*np.sin(phi) + ValorMw

#elipse 68
x2 = rx2*np.cos(Theta)*np.cos(phi) - ry2*np.sin(Theta)*np.sin(phi) + ValorMM
y2 = rx2*np.sin(Theta)*np.cos(phi) + ry2*np.cos(Theta)*np.sin(phi) + ValorMw

ax.plot(x1, y1, color = 'blue', label = '95%')
ax.plot(x2, y2, color='red', label = '68%')

plt.xlabel('OM_m')
plt.ylabel('w')


plt.show()
plt.title('Espaço de Parâmetros')
		

    
'''
def LnLike(x,data1, data2, Cinv): 
    d = len(data1) 		#número de dados coletados
    deltax = np.zeros(d)
    M = CDM.LCDModel(72., 299792.4580, 0.0, x[0], 1.0-x[0], 0.0, x[1], 0.0005) #é criado uma instância do modelo LCDM, fixando constante de hubble atualmente H0 = 72, c = 299792.458km/s, OM_r = 0.
    i=0
    M0 = x[2]	#assumindo que todas as SNIa possuem magnitude absoluta iguais
    Mbo = np.zeros(len(data2))	
    while i < len(data2):
        Mbo[i] = data2[i] #magnitude aparente observado
        deltaMb[i] = Mbo[i] - M.Mb(data1[i], M0) #diferença entre o magnitude aparentes, observado e o teórico (dado os possíveis valores dos parâmbetros).        
        i = i+1                         
    Qui2 = np.matmul(np.matmul(deltaMb, Cinv),deltaMb)    
    result = np.log(prior(x))-Qui2/2. # logaritmo da posterior ln(post) = ln(L*prior) = ln(L) + ln(prior)
    print(f'Ln(post) =  {result} \n')
    return result

# comparar as probabilidades posterior do ponto atual e possível próximo ponto no espaço de parâmetros
def Passo(xi,xp, data1, data2, Cinv, LnLike):
    print('este é o xi atualmente: ',xi, '\n') 
    LNi = LnLike(xi, data1, data2, Cinv) #LnLike do ponto inicial/atual xi
    print('este é o xp sorteado: ',xp, '\n') 
    LNp = LnLike(xp, data1, data2, Cinv) #LnLike do possível ponto posterior xp
    alpha = np.random.uniform(0.,1.) #sorteio uniforme de um número entre 0 e 1.
    if LNp > LNi: #caso a probabilidade favorece o próximo ponto xp.
        print('Lnp > Lni \n')
        return 1 #retorna 1 caso aceito o próximo ponto xp
    else: #caso a probabilidade não favoreça xp, é necessário utilizar o sorteio de alpha.
        r = np.exp(Decimal(LNp))/np.exp(Decimal(LNi)) #razão entre posterio de xp e posterior de xi
        print(f'Essa é a razão: {r}            Este o alpha: {alpha}\n') 
        if alpha < r: #critério de aceite do ponto xp
            print('alpha < r \n')
            return 1 # aceito
        else:
            print('alpha > r \n')
            return 0 #recusado









def LnLike(x,data1, data2, sig): 
    d = len(data1) 		#número de dados coletados
    deltax = np.zeros(d)
    M = CDM.LCDModel(72., 299792.4580, 0.0, x[0], 1.0-x[0], 0.0, x[1], 0.0005) #é criado uma instância do modelo LCDM, fixando constante de hubble atualmente H0 = 72, c = 299792.458km/s, OM_r = 0.
    i=0
    M0 = x[2]	#assumindo que todas as SNIa possuem magnitude absoluta iguais
    Mbo = np.zeros(len(data2))	
    deltaMb = np.zeros(len(data2))	
    while i < len(data2):
        Mbo[i] = data2[i] #magnitude aparente observado
        deltaMb[i] = Mbo[i] - M.Mb(data1[i], M0) #diferença entre o magnitude aparentes, observado e o teórico (dado os possíveis valores dos parâmbetros).        
        i = i+1                         
    Qui2 = (deltaMb/sig)**2    
    result = np.log(prior(x))-np.sum(Qui2/2.) # logaritmo da posterior ln(post) = ln(L*prior) = ln(L) + ln(prior)
    print(f'Ln(post) =  {result} \n')
    return result

# comparar as probabilidades posterior do ponto atual e possível próximo ponto no espaço de parâmetros
def Passo(xi,xp, data1, data2, sig, LnLike):
    print('este é o xi atualmente: ',xi, '\n') 
    LNi = LnLike(xi, data1, data2, sig) #LnLike do ponto inicial/atual xi
    print('este é o xp sorteado: ',xp, '\n') 
    LNp = LnLike(xp, data1, data2, sig) #LnLike do possível ponto posterior xp
    alpha = np.random.uniform(0.,1.) #sorteio uniforme de um número entre 0 e 1.
    D = LNp-LNi
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>> diferença ', D) 
    if LNp > LNi: #caso a probabilidade favorece o próximo ponto xp.
        print('Lnp > Lni \n')
        return 1 #retorna 1 caso aceito o próximo ponto xp
    else: #caso a probabilidade não favoreça xp, é necessário utilizar o sorteio de alpha.
        r = np.exp(D) #razão entre posterio de xp e posterior de xi
        print(f'Essa é a razão: {r}            Este o alpha: {alpha}\n') 
        if alpha < r: #critério de aceite do ponto xp
            print('alpha < r \n')
            return 1 # aceito
        else:
            print('alpha > r \n')
            return 0 #recusado
'''

