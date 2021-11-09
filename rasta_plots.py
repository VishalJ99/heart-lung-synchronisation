#Libraries
import scipy as sp
import pylab as plt
import numpy as np
from scipy.integrate import odeint
from utils.HH_utils import HH_utils as utils
from hh_model import HodgkinHuxley_model as hh_model

class plots:

    def __init__(self, hh_model):
        self.hh_model = hh_model



        
    def get_vals(self):

        #Gets HH variables
        hh_int_vars = self.hh_model.get_HH_vars()
        #Unpacks HH_int_vars
        self.v = hh_int_vars[:,0]
        self.m = hh_int_vars[:,1]
        self.h = hh_int_vars[:,2]
        self.n = hh_int_vars[:,3]

        #Gets and unpacks currents currents  
        self.curr_Na, self.curr_K, self.curr_L, self.input_curr = self.hh_model.get_curr_vals(hh_int_vars)

        #Gets time array
        self.t = self.hh_model.t

        #Gets Index positions of the peaks
        self.peak_idx = utils.find_max(self.v)

        #Gets voltage and Time peak vals
        self.peak_v, self.peak_t = self.v[self.peak_idx], self.t[self.peak_idx]

        #print(self.peak_v)

    def plot_action(self, curr_max = 5.9):   

        fig, ax = plt.subplots(2, 1, figsize = (15,35), sharex = True )  

        self.hh_model.set_curr_bounds(i_base = 0, i_max = curr_max)
        self.get_vals()

        #Plot action Potential
        ax[0].plot(self.t, self.v, color = 'black', label = 'Action potential')
        ax[0].set_title('Action potential --%s'%curr_max)

        #Plot maximas
        ax[0].scatter(self.peak_t, self.peak_v)

        #Plot Input current
        ax[1].plot(self.t, self.input_curr)
        plt.show()


    def plot_all(self, curr_max = 6.2):

        

        fig, ax = plt.subplots(4,1, figsize = (15,28))

        self.hh_model.set_curr_bounds(i_base = 0, i_max = curr_max)
        self.get_vals()
        print(f'm_pre = {self.m[0]}, m_post = {self.m[-1]}, n_pre = {self.n[0]}, n_post = {self.n[-1]}, h_pre = {self.h[0]}, h_post = {self.h[-1]}')
        #Plot action potential
        ax[0].plot(self.t, self.v, color = 'black', label = 'Action potential')
        ax[0].set_title('Action potential --%s'%curr_max)
        ax[0].legend()

        #Plot coeffs
        ax[1].set_title('Coeffs')
        ax[1].plot(self.t, self.m, color = 'blue', label = 'm')
        ax[1].plot(self.t, self.n, color = 'red', label = 'n')
        ax[1].plot(self.t, self.h, color = 'green', label = 'h')
        ax[1].legend()

        #Plot Ion currents
        ax[2].plot(self.t, self.curr_Na, color = 'blue', label = 'Na current')
        ax[2].plot(self.t, self.curr_K, color = 'red', label = 'K current')
        ax[2].plot(self.t, self.curr_L, color = 'green', label = 'L current')
        ax[2].legend()

        #Plot input current
        ax[3].plot(self.t, self.input_curr, label = 'Input Current')
        ax[3].legend()

        plt.show()

    def plot_cut_off(self):

        freq_list, curr_list = [] , []

        for curr in list(np.arange(5, 7, 0.01)):
            #print(curr)
            self.hh_model.set_curr_bounds(i_base = 0, i_max = curr)
            self.get_vals()

            #Checks if action potential peaks were prdoduced
            if len(self.peak_t) > 2:

                #Checks action potentials were generated across all of the step current
                time_diff = self.peak_t[-1] - self.peak_t[0]
                if time_diff >= 0.8 * self.hh_model.t_resp:

                    #Removes initial action potentials--caused by initial conditions
                    peak_times = self.peak_t[4:]

                    if len(peak_times) > 1:

                        #Calcs avg interval between peaks
                        avg_interval = (peak_times[-1] - peak_times[0]) / len(peak_times)
                        freq = 1 / avg_interval

            else:
                freq = 0

            freq_list.append(freq)
            curr_list.append(curr)

        plt.plot(curr_list,freq_list)
        plt.show()  


    def plot_rasta(self,i_min, i_max, num_vals):

            # self.hh_model.set_curr_bounds(i_base = 7, i_max = 10)
            # self.get_vals()
            peak_matrix = []


            scale = 0.1
            for curr in list(np.arange(i_min, i_max, num_vals)):

                self.hh_model.set_curr_bounds(i_base = 8, i_max = curr)

                self.get_vals()

                y = np.ones(len(self.peak_idx))
                y_plot = y * scale
                plt.scatter(self.peak_t, y_plot, color ='black', marker = 'o', s = 1)

                scale += 0.1
       

            plt.show()

model = hh_model()
plot = plots(model)
#plot.plot_action(curr_max = 6.265)
plot.plot_all()
# plot.plot_cut_off()
#plot.plot(10,20,.2)