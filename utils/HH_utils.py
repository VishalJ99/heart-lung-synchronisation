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

	@staticmethod
	def get_inp_curr(t,t_base,i_base,i_max):
		'''Returns the current based on the period--current input based on a step function'''
		
		if t//t_base%2: return i_max
		else: return i_base

	@staticmethod
	def find_max(array):
	    '''
	    Finds the turning points within an 1D array and returns the indices of the
	    maximum turning points in a list.
	    '''
	    idx_max = []

	    NEUTRAL, RISING, FALLING = range(3)
	    def get_state(a, b):
	        if a < b: return RISING 
	        elif a > b: return FALLING
	        return NEUTRAL

	    point_state = get_state(array[0], array[1])
	    begin = 1
	    for i in range(2, len(array)):
	        state = get_state(array[i - 1], array[i])
	        if state != NEUTRAL:
	            if point_state != NEUTRAL and point_state != state:
	                if state == FALLING and array[i] > 0: 
	                    idx_max.append((begin + i - 1) // 2)
	            begin = i
	            point_state = state
	    return idx_max	



