#Libraries
import scipy as sp
import pylab as plt
import numpy as np
from scipy.integrate import odeint
from HH_utils import HH_utils as utils

class HodgkinHuxley_model():

	'''
	Starting Paramters 
	'''

	m0 = 0.05
	n0 = 0.32
	h0 = 0.6
	v0 = -65
	t = np.arange(0.0, 450.0, 0.01)	

	'''
	Membrane Currents (i_element) [uA/cm^2]
	'''

	#Sodium (Na) 
	@staticmethod	
	def i_Na(voltage, m, h):

		conductance = utils.g_Na
		eqlm_voltage = utils.V_Na

		return conductance * (m**3) * h * (voltage	- eqlm_voltage) 

	#Potassium (K)
	@staticmethod	
	def i_K(voltage, n):
		
		conductance = utils.g_K
		eqlm_voltage = utils.V_K

		return conductance * (n**4) * (voltage - eqlm_voltage)

	#Leak (L) 
	@staticmethod	
	def i_L(voltage):

		conductance	= utils.g_L
		eqlm_voltage = utils.V_L

		return conductance * (voltage - eqlm_voltage)

	@staticmethod
	def input_current(t):
		'''
		For a given time (t) the fn returns the input current
		'''

		temp_mag = 10

		return temp_mag*(t>100) - temp_mag*(t>200)

	def dAlldt(self, init_conditions_list, t):
		'''
		Input Parameters:
		init_conditions_list = [v0, m0, h0, n0]
		t = list containing times of interest
		
		Returns:
		HH differential equations
		'''

		v, m, h, n = init_conditions_list
		mem_C = utils.mem_C	


		dVdt = (self.input_current(t) - self.i_Na(v, m, h) - self.i_K(v, n) - self.i_L(v)) / mem_C	
		dmdt = utils.alpha_m(v) * (1 - m) - utils.beta_m(v)*m
		dhdt = utils.alpha_h(v) * (1 - h) - utils.beta_h(v)*h
		dndt = utils.alpha_n(v) * (1 - n) - utils.beta_n(v)*n
		

		return dVdt, dmdt, dhdt, dndt 


	def Main(self):
		
		init_conditions_list = [self.v0, self.m0, self.h0, self.n0]

		hh_integral_matrix = odeint(self.dAlldt, init_conditions_list, self.t)

		#Unpacks the voltage and m,h,n values at each time step
		v = hh_integral_matrix[:,0]
		m = hh_integral_matrix[:,1]
		h = hh_integral_matrix[:,2]
		n = hh_integral_matrix[:,3]

		curr_Na = self.i_Na(v, m ,h)
		curr_K = self.i_K(v, n)
		curr_L = self.i_L(v)

		plt.plot(self.t, v, 'k')
		plt.show()



model = HodgkinHuxley_model()
model.Main()
#model.i_Na()	