import math
import numpy as np
import CDM 
import matplotlib.pyplot as plt
from decimal import Decimal
import time

#o problema e gerar uma distribuiçao no espaço dos parametros OM_m(materia) e w(dark energy state equation) dado o modelo LCDM (usando o pacote que escrevi CDM). Os parametros sao denotados por x = [OM_m, w].

#aqui defino a distribuiçao de transiçao. Pega os dois parametros, guardados no vetor x, e devolve outros dois parametros de acordo com a distribuiçao normal, centrada em cada um dos parametros iniciais:

Transition = lambda x: [np.random.normal(x[0], 0.05), np.random.normal(x[1], 0.05)]

#Definir a distribuiçao prior para cada calor dos parametros armazenado no vetor x
    
def prior(x):
# coloquei aqui a informaçao de que a probabilidade deve se anular caso algum dos parametros se anule ou se torne negativo.
	if x[0]<= 0.: 
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
    M = CDM.LCDModel(68., 299792.4580, 0.0, x[0], 1.0-x[0], 0.0, x[1]) #é criado uma instância do modelo LCDM, fixando constante de hubble atualmente H0 = 72, c = 299792.458km/s, OM_r = 0.
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
    



#IMPLEMENTAÇÃO DO MÉTODO MCMC


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


#chute inicial dos parâmetros para a cadeia
param_init0 = float(input(print('digite um chute inicial para Omega_m: \n'))) 
param_init1 = float(input(print('digite um chute inicial para w: \n')))
iterations = int(input(print('digite o numero de iterações: \n')))    
sig = 0.4



xi = np.array([param_init0, param_init1])    #ponto inicial da cadeia
#lista com os valores aceitos dos parâmetros. Inicialmente há apenas um valor
chain_m = [xi[0]]	
chain_w = [xi[1]]
#lista de valores sorteados para os parâmetros porém rejeitados 
chain_mrej = []
chain_wrej = []
accepted = [0] #numeros de sorteios aceitados
rejected = [0] #numeros de sorteios rejeitados

ini = time.time()
for i in range(iterations): #iteração i
    print('iteração ', i, '\n') #monitoramento da contagem   
    xp = Transition(xi) #sorteio de um novo ponto, xp, a partir do ponto atual xi
    ret = Passo(xi,xp, data1, data2, sig, LnLike) # verificar a condição de aceite de xp. Note que assumirá valores 1 (aceite) ou 0 (recusa). 
    print('retrun do passo: ', ret, '\n')
    if ret == 1: #se for aceito
        xi = xp #atualiza o ponto atual para o ponto sorteiado xp.
        accepted[0] = accepted[0] + 1 # atualiza quantos foram aceitos. 
    else: #não aceito
        rejected[0] = rejected[0] + 1 # atualiza quantos foram rejeitados. continua com o ponto xi.
	
	# insere na lista dos recusados xp
        chain_mrej.append(xp[0]) 
        chain_wrej.append(xp[1])

    print(' true: ', accepted, '\n\n')
    print(' false: ', rejected, '\n\n')

#acrescenta os valores atualizados à cadeia.
    chain_m.append(xi[0])
    chain_w.append(xi[1])
    
#após todas as iterações cálculo do valor médio    
n = len(chain_m)
ValorMM = sum(chain_m)/n
ValorMw = sum(chain_w)/n

#montando a matrix covariancia dos parâmetros
cov = np.zeros([2,2])
a = cov[0,0] = sum((chain_m-ValorMM)**2)/(n-1)
b = cov[0,1]= cov[1,0] = sum((chain_m-ValorMM)*(chain_w-ValorMw))/(n-1)
c = cov[1,1] = sum((chain_w-ValorMw)**2)/(n-1)

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







fim = time.time()
tempoT = (fim - ini)/60. 
print(f'\n o valor esperado dos parametros sao OM_m = {ValorMM}, w = {ValorMw} e o tempo de execução total das iterações {tempoT} min')

#plotar os histogramas
plt.hist(chain_m, bins = 100, label= 'Matter', color= 'blue')
plt.hist(chain_w, bins = 100, label= 'w', color= 'red')
plt.legend()
plt.show()



#pontos percorridos no espaço de parâmetros


plt.rcParams['legend.fontsize'] = 10
fig = plt.figure()
ax = fig.gca()
plt.scatter(chain_mrej, chain_wrej, marker='x',color='y',linewidths=1) # recusados em amarelo (yellow)
plt.scatter(chain_m, chain_w,marker='.',color='b',linewidths=3)# aceitos em azul (blue)
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
talvez seja necessário para a prior:

if x[0] + x[1] > 1.0: 
		return 0.0
'''

