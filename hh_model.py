#Libraries
import scipy as sp
import pylab as plt
import numpy as np
from scipy.integrate import odeint
from utils.HH_utils import HH_utils as utils

class HodgkinHuxley_model():

	def __init__(self):
		self.m0 = 0.05
		self.n0 = 0.32
		self.h0 = 0.6
		self.v0 = -65
		self.set_time()
		self.set_curr_bounds()
		#self.t = np.arange(0.0, 450.0, 0.01)	

	

	def set_time(self, t_resp = 400, n_resps = 2, t_points = 100000):
		'''Creates an array of time periods'''
		
		self.t_resp = t_resp
		max_t = t_resp*n_resps
		self.t = np.linspace(0, max_t, t_points )

	def set_curr_bounds(self, i_base = 5, i_max = 10):
		'''Sets the base and max input current for the currrent step function'''

		self.i_base = i_base
		self.i_max = i_max


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


		dVdt = (utils.get_inp_curr(t, self.t_resp, self.i_base, self.i_max) - self.get_i_Na(v, m, h) - self.get_i_K(v, n) - self.get_i_L(v)) / mem_C	
		dmdt = utils.alpha_m(v) * (1 - m) - utils.beta_m(v)*m
		dhdt = utils.alpha_h(v) * (1 - h) - utils.beta_h(v)*h
		dndt = utils.alpha_n(v) * (1 - n) - utils.beta_n(v)*n
		

		return dVdt, dmdt, dhdt, dndt 

	def get_HH_vars(self):
		'''Returns a matrix of all the integrated hh variables'''

		init_conditions_list = [self.v0, self.m0, self.h0, self.n0]

		hh_integral_matrix = odeint(self.dAlldt, init_conditions_list, self.t)
		return hh_integral_matrix

	def get_curr_vals(self, hh_integral_matrix):
			'''Returns a matrix of all the currents over time'''

			#Unpacks the voltage and m,h,n values at each time step
			v = hh_integral_matrix[:,0]
			m = hh_integral_matrix[:,1]
			h = hh_integral_matrix[:,2]
			n = hh_integral_matrix[:,3]

			#Calcs the currents
			curr_Na = self.get_i_Na(v, m ,h)
			curr_K = self.get_i_K(v, n)
			curr_L = self.get_i_L(v)
			input_curr = [utils.get_inp_curr(time, self.t_resp, self.i_base, self.i_max) for time in self.t]

			return curr_Na, curr_K, curr_L, input_curr



model = HodgkinHuxley_model()
#model.plot_vars()

# model.set_time(t_resp = 500, n_resps = 10, t_points = 4000)
# model.plot_vars()
#model.i_Na()	


