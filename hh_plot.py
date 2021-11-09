#Libraries
import scipy as sp
import pylab as plt
import numpy as np
from scipy.integrate import odeint
from utils.HH_utils import HH_utils as utils
from hh_model import HodgkinHuxley_model as hh_model

class plot_hh():

	def plot_vars(self):


		hh_integral_matrix = hh_model.get_HH_vars()

		#Unpacks the voltage and m,h,n values at each time step
		v = hh_integral_matrix[:,0]
		m = hh_integral_matrix[:,1]
		h = hh_integral_matrix[:,2]
		n = hh_integral_matrix[:,3]

		curr_Na = hh_model.get_i_Na(v, m ,h)
		curr_K = hh_model.get_i_K(v, n)
		curr_L = hh_model.get_i_L(v)

		fig, ax = plt.subplots(5, 1, figsize = (15,35), sharex = True )

		ax[0].plot(hh_model.t, v, color = 'black', label = 'Action potential')
		ax[0].set_title('Action potential')

		ax[1].set_title('Variable')
		ax[1].plot(hh_model.t, m, color = 'red', label = 'm variable')
		ax[1].plot(hh_model.t, h, color = 'blue', label = 'h variable')
		ax[1].plot(hh_model.t, n, color = 'green', label = 'n variable')
		ax[1].legend(fontsize = 17)

		ax[2].plot(hh_model.t, curr_Na, color = 'black')
		ax[2].set_title('Na Current')

		ax[3].plot(hh_model.t, curr_K, color = 'black')
		ax[3].set_title('K Current')
		
		ax[4].plot(hh_model.t, curr_L, color = 'black')
		ax[4].set_title('Leak Current')

		#plt.savefig('Inital_results.png')
		plt.show()




