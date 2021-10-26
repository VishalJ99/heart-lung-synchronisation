#Libraries
import numpy as np 


class HH_utils():
		
	'''
	Utility package for the Hodgkin-Huxley class
	'''

	#Membrane Capicitance [uF/cm^2]
	mem_C  =   1.0

	#Sodium (Na) maximum conductance [mS/cm^2]
	g_Na = 120.0

	#Potassium (K) maximumum conductance [mS/cm^2]
	g_K  =  36.0

	#Leak maximum conductance [mS/cm^2]
	g_L  =   0.3

	#Sodium (Na) Nernst reversal potentials [mV]
	E_Na =  50.0
   
	#Postassium (K) Nernst reversal potentials [mV]
	E_K  = -77.0

	#Leak Nernst reversal potentials [mV]
	E_L  = -54.387
	
	'''
	Channel Gating kinetic coefficients
	'''

	@staticmethod
	def alpha_m(voltage):
		return  0.1*(voltage + 40)  / ( 1 - np.exp(- (voltage + 40) / 10))

	@staticmethod	
	def beta_m(voltage):
		return 4*np.exp(- (voltage + 65)/18)

	@staticmethod	
	def alpha_n(voltage):
		return 0.01*(voltage + 55) / (1 - np.exp(-(voltage + 55)/10))

	@staticmethod	
	def beta_n(voltage):
		return 0.125*np.exp(-(voltage + 65) / 80)

	@staticmethod	
	def alpha_h(voltage):
		return 0.07*np.exp(-(voltage + 65) / 20)

	@staticmethod	
	def beta_h(voltage):
		return 1 / 	(1 + np.exp(-(voltage + 35)/10))



