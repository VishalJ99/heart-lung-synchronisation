import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import sys
import matplotlib.pyplot as plt

def find_firing_freq(I):
	try:
		X = np.load('filtered_Is.npy')[:-1]
		Y = np.load('filtered_Fs.npy')[:-1]
		f = interp1d(X,Y,'cubic')

	except FileNotFoundError:
		print('No npy files found in working dir to load, func will not load')
		sys.exit(1)
 
	return f(I)








