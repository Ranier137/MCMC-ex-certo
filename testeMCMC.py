import math
import numpy as np
import CDM 
import matplotlib.pyplot as plt
import time
import pylab as pl

#o problema e gerar uma distribuiçao no espaço dos parametros OM_m(materia) e w(dark energy state equation) dado o modelo LCDM (usando o pacote que escrevi CDM). Os parametros sao denotados por x = [OM_m, w].

#aqui defino a distribuiçao de transiçao. Pega os dois parametros, guardados no vetor x, e devolve outros dois parametros de acordo com a distribuiçao normal, centrada em cada um dos parametros iniciais:

Transition = lambda x: [np.random.normal(x[0], 0.05), np.random.normal(x[1], 0.05)]

#Definir a distribuiçao prior para cada calor dos parametros armazenado no vetor x
    
def prior(x):
# coloquei aqui a informaçao de que a probabilidade deve se anular caso algum dos parametros se anule ou se torne negativo.
	if x[0]< 0.: 
		return 0.0
	if x[1] > -1./3.:
		return 0.0
#usei prior uniforme. coloquei apenas 1.0 pq nao importa a constante de normalizaçao que deve aparecer, nao influenciara no metodo.	
	else:
		return 1.0

# aqui devo defini o logaritmo da likelihood*prior
#alem dos parmetros x, as entradas sao sig = 0.4 desvio padrao (assumindo covariancia diagonal) e dadaos representados por data1 = coluna dos redshifts e data2 = coluna das respectivas magnitude aparente.
def LnLike(x,data1, data2,sig = 0.4): 
    d = len(data1) 		#número de dados coletados
    deltax = np.zeros(d)
    M = CDM.LCDModel(72., 299792.4580, 0.0, x[0], 1.0-x[0], 0.0, x[1]) #é criado uma instância do modelo LCDM, fixando constante de hubble atualmente H0 = 72, c = 299792.458km/s, OM_r = 0.
    i=0
    Mo = -19.3	#assumindo que todas as SNIa possuem magnitude absoluta iguais
    MI2 = np.zeros(len(data2))	
    while i < len(data2):
        MI2[i] = data2[i] - Mo #módulo de distância observado
        deltax[i] = MI2[i] - M.MI(data1[i]) #diferença entre o modulo de distância observado e o teórico (dado os possíveis valores dos parâmetros).        
        i = i+1                         
    result = np.log(prior(x))-np.sum((deltax**2)/(2*sig**2)) # logaritmo da posterior ln(post) = ln(L*prior) = ln(L) + ln(prior)
    print(f'Ln(post) =  {result} \n')
    return result

# comparar as probabilidades posterior do ponto atual e possível próximo ponto no espaço de parâmetros
def Passo(xi,xp, data1, data2, sig, LnLike):
    print('este é o xi atualmente: ',xi, '\n') 
    LNi = LnLike(xi, data1, data2, sig) #LnLike do ponto inicial/atual xi
    print('este é o xp sorteado: ',xp, '\n') 
    LNp = LnLike(xp, data1, data2, sig) #LnLike do possível ponto posterior xp
    alpha = np.random.uniform(0.,1.) #sorteio uniforme de um número entre 0 e 1.
    if LNp > LNi: #caso a probabilidade favorece o próximo ponto xp.
        f = open('Ln.txt', 'a')
        f.write(str(LNp) + ' ' + str(xp[1]) + ' ' + str(xp[0]) + '\n')
        f.close()
        print('Lnp > Lni \n')
        return 1 #retorna 1 caso aceito o próximo ponto xp
    else: #caso a probabilidade não favoreça xp, é necessário utilizar o sorteio de alpha.
        r = np.exp(LNp-LNi) #razão entre posterio de xp e posterior de xi
        print(f'Essa é a razão: {r}            Este o alpha: {alpha}\n') 
        if alpha < r: #critério de aceite do ponto xp
            print('alpha < r \n')
            f = open('Ln.txt', 'a')
            f.write(str(LNp) + ' ' + str(xp[1]) + ' '+ str(xp[0]) +'\n')
            f.close()
            return 1 # aceito
        else:
            print('alpha > r \n')
            f = open('Ln.txt', 'a')
            f.write(str(LNi) + ' ' + str(xi[1]) + ' '+ str(xi[0]) +'\n')
            f.close()
            return 0 #recusado
    



######################################################



# Criando data1 e data2, ainda como listas vazias.
coloumn0 = []
coloumn1 = []

#abrir o arquivo com os dados
with open(r"DAT.txt", "r+") as f:
    data = f.readlines() #ler linha por linha


# No arquivo DAT.txt os dados estão separados em colunas. essas colunas estão separados por dois espaços. Para organizar e inserir os dados da primeira coluna na lista column0 é feito:

    for line in data:
        coloumn0.append(line.strip().split("   ")[0])
        #strip para remover  \n
        #split todo intervalo entre colunas
        #second element is indexed 1

#fazendo o mesmo para a segunda coluna, acrescentando seus valores em coloumn1

    for line in data:
        coloumn1.append(line.strip().split("   ")[1])
        
#Transformando as listas em np.arrays
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

####################################################
m = []    
w = []
nm = 100
nw = 100
dm = 0.9/nm
dw = 29./(3.*nw)

i=0
while i < nm+1:
    m.append(0.1 + (i)*dm)
    i=i+1    
j=0 
while j< nw+1:
    w.append(-10 + (j)*dw)
    j=j+1

n = 256 
x = np.array([float(k) for k in m]) 
y = np.array([float(k) for k in w])

    
Z = np.zeros([nm,nw])



 
    


X, Y = np.meshgrid(x, y) 

ini = time.time()
i=0
j=0
while i < nm:
    while j < nw:    
        Z[i,j] = LnLike([x[i],y[j]],data1, data2,sig = 0.4)
        print(nw*nm - (i+1)*(j+1))
        j=j+1
    j=0
    i=i+1

fim = time.time()

T = fim - ini

print(Z)

print(T/60., 'min')
#print(np.meshgrid(x,y))



plt.xlim(0.1, 1.)
plt.ylim(-10, -1./3.)
plt.pcolormesh(X, Y, Z)
plt.colorbar(label='scale') 
plt.clim(-300, -100)
plt.show()

