import math
import scipy.integrate
import numpy as np

#modelo LambdaCDM

class LCDModel:
    
    #construtor das instâncias
	def __init__(self, H0, c, OM_r, OM_m, OM_de, w):
		self.H0 = H0
		self.c = c
		self.OM_m = abs(OM_m)
		self.OM_de = abs(OM_de)
		self.OM_r = abs(OM_r)
		self.w = w
    #parâmetro de Hubble/H0
	def e(self, z):
		#OM_k = 1.0 - self.OM_r - self.OM_m - self.OM_de      
		return math.sqrt(self.OM_r*(1.0+z)**4 + self.OM_m*(1.0+z)**3 + self.OM_de*(1+z)**(3*(1 + self.w)))
    
	#comoving distance
	def r(self, z):
        	f = lambda z1: 1.0/self.e(z)
        	I = (self.c/self.H0)*scipy.integrate.quad(f,0,z)[0]        
        	return I

    #Angular Diameter Distance
	def DA(self, z):
        	return self.r(z)/(1.0 + z)


    #Luminosity distance
	def DL(self, z):
        	return (1.0 + z)*self.r(z)


    #módulo de distância (supernovas), isto é mb - M
	def MI(self, z):
        	return 5*math.log(self.DL(z),10) + 25.0
    
   



