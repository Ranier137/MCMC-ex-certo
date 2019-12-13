import math
import scipy.integrate
import numpy as np

#modelo LambdaCDM

class LCDModel:
    
    #construtor das instâncias
    def __init__(self, H0, c, OM_r, OM_m, OM_de, OM_k, w, tol = 0.0):
        self.H0 = H0
        self.c = c
        self.OM_m = abs(OM_m)
        self.OM_de = abs(OM_de)
        self.OM_r = abs(OM_r)
        self.OM_k = OM_k		
        self.w = w
        self.tol = abs(tol)
    #parâmetro de Hubble/H0
 
    def e(self, z):
        #OM_k = 1.0 - self.OM_r - self.OM_m - self.OM_de      
        return (self.OM_r*(1.0+z)**4 + self.OM_m*(1.0+z)**3 + self.OM_k*(1+z)**2 +self.OM_de*(1+z)**(3*(1 + self.w)))**(0.5)
    
	
    def I(self,z):
        f = lambda z1: 1.0/self.e(z1)
        return scipy.integrate.quad(f,0,z)[0]  

    #comoving distance
    def r(self, z):
        return (self.c/self.H0)*self.I(z)  

    #Angular Diameter Distance
	#def DA(self, z):
     #   	return self.r(z)/(1.0 + z)


    #Luminosity distance
    def DL(self, z):
        if self.OM_k > self.tol:
            #print('integra +')
            return (1+z)*(self.c/self.H0)*np.sinh(self.I(z)*((abs(self.OM_k))**(0.5)))/((abs(self.OM_k))**(0.5))
        if self.OM_k < -self.tol:
            #print('integra -')
            return (1+z)*(self.c/self.H0)*np.sin(self.I(z)*((abs(self.OM_k))**(0.5)))/((abs(self.OM_k))**(0.5))
        if abs(self.OM_k) <= self.tol:
            #print('integra 0')
            return (1+z)*(self.c/self.H0)*self.I(z)

    #módulo de distância (supernovas), isto é mb - M
    def MI(self, z):
        	return 5*math.log(self.DL(z),10) + 25.0

    def Mb(self,z, M0):
            return self.MI(z) + M0
    
   



