#Libraries
import scipy as sp
import pylab as plt
import numpy as np
from scipy.integrate import odeint
from utils.HH_utils import HH_utils as utils

class HodgkinHuxley_model():


	m0 = 0.05
	n0 = 0.32
	h0 = 0.6
	v0 = -65
	t = np.arange(0.0, 450.0, 0.01)	

	
	@staticmethod	
	def get_i_Na(voltage, m, h):
		'''Returns Sodium (Na) Current'''

		return utils.g_Na * (m**3) * h * (voltage	- utils.E_Na) 

	
	@staticmethod	
	def get_i_K(voltage, n):
		'''Returns Potassium Current'''

		return utils.g_K * (n**4) * (voltage - utils.E_K)

	@staticmethod	
	def get_i_L(voltage):
		'''Returns Leak Current'''

		return utils.g_L * (voltage - utils.E_L)

	@staticmethod
	def get_input_current(t):
		'''For a given time (t) the fn returns the input current'''

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


		dVdt = (self.get_input_current(t) - self.get_i_Na(v, m, h) - self.get_i_K(v, n) - self.get_i_L(v)) / mem_C	
		dmdt = utils.alpha_m(v) * (1 - m) - utils.beta_m(v)*m
		dhdt = utils.alpha_h(v) * (1 - h) - utils.beta_h(v)*h
		dndt = utils.alpha_n(v) * (1 - n) - utils.beta_n(v)*n
		

		return dVdt, dmdt, dhdt, dndt 

	def get_HH_vars(self):

		init_conditions_list = [self.v0, self.m0, self.h0, self.n0]

		hh_integral_matrix = odeint(self.dAlldt, init_conditions_list, self.t)
		return hh_integral_matrix



	def plot_vars(self):

		hh_integral_matrix = self.get_HH_vars()

		#Unpacks the voltage and m,h,n values at each time step
		v = hh_integral_matrix[:,0]
		m = hh_integral_matrix[:,1]
		h = hh_integral_matrix[:,2]
		n = hh_integral_matrix[:,3]

		curr_Na = self.get_i_Na(v, m ,h)
		curr_K = self.get_i_K(v, n)
		curr_L = self.get_i_L(v)

		fig, ax = plt.subplots(5, 1, figsize = (15,35), sharex = True )

		ax[0].plot(self.t, v, color = 'black', label = 'Action potential')
		ax[0].set_title('Action potential')

		ax[1].set_title('Variable')
		ax[1].plot(self.t, m, color = 'red', label = 'm variable')
		ax[1].plot(self.t, h, color = 'blue', label = 'h variable')
		ax[1].plot(self.t, n, color = 'green', label = 'n variable')
		ax[1].legend(fontsize = 17)

		ax[2].plot(self.t, curr_Na, color = 'black')
		ax[2].set_title('Na Current')

		ax[3].plot(self.t, curr_K, color = 'black')
		ax[3].set_title('K Current')
		
		ax[4].plot(self.t, curr_L, color = 'black')
		ax[4].set_title('Leak Current')

		#plt.savefig('Inital_results.png')
		plt.show()


model = HodgkinHuxley_model()
model.plot_vars()
#model.i_Na()	


