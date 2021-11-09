#Libraries
import numpy as np 


class HH_utils():
		
	'''
	Utility package for the Hodgkin-Huxley class
	'''

	#Membrane Capicitance [uF/cm^2]
	mem_C  =   1.0

	#Sodium (Na) maximum conductance [mS/cm^2]
	g_Na = 40

	#Potassium (K) maximumum conductance [mS/cm^2]
	g_K  =  35

	#Leak maximum conductance [mS/cm^2]
	g_L  =   0.3

	#Sodium (Na) Nernst reversal potentials [mV]
	E_Na =  55
   
	#Postassium (K) Nernst reversal potentials [mV]
	E_K  = -77.0

	#Leak Nernst reversal potentials [mV]
	E_L  = -65
	
	
	'''
	Channel Gating kinetic coefficients
	'''

	@staticmethod	
	def alpha_n(voltage):
		return 0.02*(voltage - 25) / (1 - np.exp(-(voltage - 25)/9))

	@staticmethod	
	def beta_n(voltage):
		return -0.002*(voltage - 25) / (1 - np.exp((voltage -25) / 9))

	@staticmethod
	def alpha_m(voltage):
		return  0.182*(voltage + 35)  / ( 1 - np.exp(- (voltage + 35) / 9))

	@staticmethod	
	def beta_m(voltage):
		return -0.124*(voltage + 35) / (1 - np.exp((voltage + 35)/9))

	@staticmethod	
	def alpha_h(voltage):
		return 0.25*np.exp(-(voltage + 90) / 12)

	@staticmethod	
	def beta_h(voltage):
		return 0.25* np.exp((voltage + 62)/6) / (np.exp((voltage + 90)/12))

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
	    #Fn to see if the line is rising/falling/neutral
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

	            	#Only interested in maximas above zero
	                if state == FALLING and array[i] > 0: 
	                    idx_max.append((begin + i - 1) // 2)

	            begin = i
	            point_state = state
	    return idx_max	



