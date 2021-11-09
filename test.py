#Libraries
import scipy as sp
import pylab as plt
import numpy as np
from scipy.integrate import odeint
from utils.HH_utils import HH_utils as utils
from hh_model import HodgkinHuxley_model as hh_model



model = hh_model()

for current in range(10, 13):
	#model.set_time(t_resp = 500, n_resps = 10, t_points = 4000)
	model.set_curr_bounds(i_base = 7, i_max = current )
	model.plot_vars()