import scipy as sp
import pylab as plt
from scipy.integrate import odeint
import numpy as np
from scipy.signal import find_peaks
from find_neuron_firing_freq import find_firing_freq
import pandas as pd 
import sys
from icecream import ic
import os
class HodgkinHuxleyModel():
    '''
    
    Hodgekin Huxley song bird (?) neuron 

    ...

    Attributes
    ----------
   
    C_m : float
        membrane capacitance (uF/cm^2)
    
    g_Na : float
        maximum sodium channel conductance (mS/cm^2)
    g_K : float
        maximum potassium conductances (mS/cm^2)
    g_L : float
        maximum leak channel conductances (mS/cm^2)
    
    E_Na : float
        Nernst reversal potential for sodium channel (mV)
    E_K : float
        Nernst reversal potential for potassium channel (mV) 
    E_L : float
        Nernst reversal potential for leak channel (mV)

    f_base : float
        Base firing rate of the neuron corresponding to the base current value (mHz)     
    
    t_resp : float
        period of current (ms?)
    i_base : float
        base value for current
    rsa : float
        rsa value which defines i_max by relation i_max = (1 + rsa) * i_base
    dc : float
        Duty cycle of signal
    
    t : np.ndarray 
        Time array 
    t_max : float
        maximum time during which HH equations are being integrated across
    num_time_steps: int 
        defines number of values in time array



    Methods
    -------
    alpha_m(V):
        Calculates close to open transition rate for sodium activation gating variable at a specific voltage
    alpha_n(V):
        Calculates close to open transition rate for potassium activation gating variable at a specific voltage
    alpha_h(V):
        Calculates close to open transition rate for sodium inactivation gating variable at a specific voltage
    
    beta_m(V):
        Calculates open to close transition rate for sodium activation gating variable at a specific voltage
    beta_n(V):
        Calculates open to close transition rate for potassium activation gating variable at a specific voltage
    beta_h(V):
        Calculates open to close transition rate for sodium inactivation gating variable at a specific voltage
            
    I_in(t):
        Calculates the injected current at time t
    I_Na(V, m, h):
        Calculates the sodium channel current at time t 
    I_K(V, n):
        Calculates the potassium channel current at time t
    I_L(V):
        Calculates the leak channel current at time t

      
    dALLdt(X, t):
        Calculates values for the HH ODEs at time t

    forward():
        Integrates the HH ODEs using scipy odeint, returns integrated dictionary with keys = time and values = integrated voltage values.   
    '''
    def __init__(self,t_resp,i_base=2.32,rsa=0.3,dc = 0.5,t_max=1000,num_time_steps=10000 ):
        """
        Sets the necessary attributes for the HH neuron

        Parameters
        ----------
            t_resp : float
                period of current (ms)
            i_base : float
                base value for current (nA?)
            rsa : float
                rsa value which defines i_max by relation i_max = (1 + rsa) * i_base 
            dc : float 
                Duty cycle for injected current signal
            t_max : float
                maximum time during which HH equations are being integrated across
            num_time_steps: int 
                defines number of values in time array
        """
        # membrane capacitance, in uF/cm^2
        self.C_m = 1.0
        # maximum conductances, in mS/cm^2
        self.g_Na = 69.0
        self.g_K = 6.9
        self.g_L =   0.165
        # Nernst reversal potentials, in mV
        self.E_Na =  41.0
        self.E_K  = -100.0
        self.E_L  = -65
        # set time step array
        self.t = np.linspace(0,t_max,num_time_steps)
        # properties of respiratory step signal 
        self.t_resp = t_resp
        self.i_base = i_base
        self.f_base = find_firing_freq(i_base) / 10**3 # divide by 10**3 since freq returned in hz and working time unit is ms
        self.rsa = rsa
        self.dc = dc

    def alpha_m(self, V):
        """        
        Prints the person's name and age.

        If the argument 'additional' is passed, then it is appended after the main info.

        Parameters
        ----------
        additional : str, optional
            More info to be displayed (default is None)

        Returns
        -------
        None
        
        """
        V_t = -39.92
        dV = 10
        dVt = 23.39
        t0 = 0.143
        taueps = 1.099
        
        thetai = (V - V_t) / dV 
        thetait = (V - V_t) /dVt

        tauj = t0 + taueps * (1 - np.tanh(thetait)**2)
        
        return 0.5*(1 + np.tanh(thetai))/tauj

    def beta_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_t = -39.92
        dV = 10
        dVt = 23.39
        t0 = 0.143
        taueps = 1.099
        
        thetai = (V - V_t) / dV 
        thetait = (V - V_t) /dVt

        tauj = t0 + taueps * (1 - np.tanh(thetait)**2)
        
        return 0.5*(1 - np.tanh(thetai))/tauj

    def alpha_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_t = -65.37
        dV = -17.65
        dVt = 27.22
        t0 = 0.701
        taueps = 12.9
        
        thetai = (V - V_t) / dV 
        thetait = (V - V_t) /dVt

        tauj = t0 + taueps * (1 - np.tanh(thetait)**2)
        
        return 0.5*(1 + np.tanh(thetai))/tauj

    def beta_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_t = -65.37
        dV = -17.65
        dVt = 27.22
        t0 = 0.701
        taueps = 12.9
        
        thetai = (V - V_t) / dV 
        thetait = (V - V_t) /dVt

        tauj = t0 + taueps * (1 - np.tanh(thetait)**2)
        
        return 0.5*(1 - np.tanh(thetai))/tauj

    def alpha_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_t = -34.58
        dV = 22.17
        dVt = 23.58
        t0 = 1.291
        taueps = 4.314
        
        thetai = (V - V_t) / dV 
        thetait = (V - V_t) /dVt

        tauj = t0 + taueps * (1 - np.tanh(thetait)**2)
        
        return 0.5*(1 + np.tanh(thetai))/tauj

    def beta_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_t = -34.58
        dV = 22.17
        dVt = 23.58
        t0 = 1.291
        taueps = 4.314
        
        thetai = (V - V_t) / dV 
        thetait = (V - V_t) /dVt

        tauj = t0 + taueps * (1 - np.tanh(thetait)**2)
        
        return 0.5*(1 - np.tanh(thetai))/tauj

    def I_Na(self, V, m, h):
        """
        Sodium membrane current (in uA/cm^2)
        """
        return self.g_Na * m**3 * h * (V - self.E_Na)

    def I_K(self, V, n):
        """
        Potassium membrane current (in uA/cm^2)
        """
        return self.g_K  * n**4 * (V - self.E_K)
    #  Leak
    def I_L(self, V):
        """
        Leak membrane current (in uA/cm^2)
        """
        return self.g_L * (V - self.E_L)

    def I_in(self,t):
        ''' Returns the current at time t
            Current is modeled of a periodic step function with a duty cycle (dc) 
        '''
        if (t/self.t_resp)%1 <= self.dc:
            return self.i_base+self.rsa
        else:
            return self.i_base

    @staticmethod
    def dALLdt(X, t, self):
        """
        Integrate HH model ODEs 
        |returns at membrane potential & activation variables at timestep t
        """
        V, m, h, n = X

        dVdt = (self.I_in(t) - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_L(V)) / self.C_m
        dmdt = self.alpha_m(V)*(1.0-m) - self.beta_m(V)*m
        dhdt = self.alpha_h(V)*(1.0-h) - self.beta_h(V)*h
        dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n
        return dVdt, dmdt, dhdt, dndt
    
    def forward(self):
        """
        Integrate HH ODE's
        """
        X = odeint(self.dALLdt, [-65.01324, 0.00657, 0.4899, 0.06034], self.t, args=(self,))
        self.V = X[:,0]
        self.m = X[:,1]
        self.h = X[:,2]
        self.n = X[:,3]

        return self.V,self.t, self.m, self.h, self.n


class HodgkinHuxleyModelPlots(HodgkinHuxleyModel):
    '''
    Class for generating various plots
    
    __init__ parameters: 


    Methods:
    forward method integrates eqs and returns voltage (mV) list and time (ms?) list. 
    '''
    def gen_freq_response_curve(self, I_min=0,I_max=12.5,step = 0.05,savefig=False):
        '''
        Generate a plot of the frequency response curve of the neuron, I_min and I_max are currents in nA, step defines the increment step size
        if plot to be saved, specify savefig arguement as a string to the save path
        '''
        # current and frequency arrays
        Is = []
        fs = []
        # set rsa to 0 to generate constant currents
        self.rsa = 0
        # loop through current range specified
        for i_in in np.arange(I_min,I_max,step):
            print('Current I_in:',i_in)
            self.i_base = i_in
            V,t = self.forward()
            peaks,_ = find_peaks(V,height=0)
            if len(peaks)!=0:
                times = t[peaks]
                T_heart = (times[-1] - times[0]) / len(times)
                fs.append(10**3/T_heart)
            else:
                fs.append(0)
            Is.append(i_in)
        plt.xlabel('Input current (nA)')
        plt.ylabel('Firing frequency (Hz)')
        plt.plot(Is,fs)
        # plt.scatter(X,Y,marker='o',s=1,c='k')
        if savefig: plt.savefig(savefig)
        plt.show()

        
    def gen_action_potential_plot(self, show_gating_and_channel_currs=False, save_fig=False):
        ''' returns a plot of the action potential against time along with the gating variables and injected current
        if plot to be saved, specify savefig argument as the save path string'''
        # integrate equations
        self.forward()
        ina = self.I_Na(self.V, self.m, self.h)
        ik = self.I_K(self.V, self.n)
        il = self.I_L(self.V)

        plt.figure()
        if show_gating_and_channel_currs: plt.subplot(4,1,1)
        else: plt.subplot(2,1,1)

        plt.title('Hodgkin-Huxley Neuron')
        plt.plot(self.t, self.V, 'k')
        plt.ylabel('V (mV)')

        # plt.subplot(4,1,2)
        # peak_v_indices,_ = find_peaks(self.V,height=0)
        # peak_times = self.t[peak_v_indices]
        # peak_time_diff_list = np.roll(peak_times,-1) - peak_times 
        # # set last instantaneous period to the second last value in array as it cant be calculated
        # peak_time_diff_list[-1] = peak_time_diff_list[-2] 
        # plt.plot(peak_time_diff_list,'k')
        # plt.ylabel('T_heart (ms)')

        # plt.subplot(4,1,3)
        # peak_time_diff_diff_list = np.roll(peak_time_diff_list,-1) - peak_time_diff_list
        # # set last instantaneous period to the second last value in array as it cant be calculated
        # peak_time_diff_diff_list[-1] = peak_time_diff_diff_list[-2]
        # plt.plot(peak_time_diff_diff_list,'k')
        # plt.ylabel('delta T_heart (ms)')

        if show_gating_and_channel_currs:
            plt.subplot(4,1,2)
            plt.plot(self.t, ina, 'c', label='$I_{Na}$')
            plt.plot(self.t, ik, 'y', label='$I_{K}$')
            plt.plot(self.t, il, 'm', label='$I_{L}$')
            plt.ylabel('Current (nA)')
            plt.legend()
            
            plt.subplot(4,1,3)
            plt.plot(self.t, self.m, 'r', label='m')
            plt.plot(self.t, self.h, 'g', label='h')
            plt.plot(self.t, self.n, 'b', label='n')
            plt.ylabel('Gating Variables')
            plt.legend()

        plt.subplot(2,1,2)

        i_inj_values = [self.I_in(t,) for t in self.t]
        plt.plot(self.t, i_inj_values, 'k')
        plt.xlabel('t (ms)')
        plt.ylabel('$I_{inj}$ (nA)')
        plt.ylim(np.amin(i_inj_values)-1, np.amax(i_inj_values)+1)
        

        if save_fig: plt.savefig(save_fig)
        plt.show()

    def gen_arnold_tongue_plot_data(self,steady_state_time = 0, t_step_size=500,rel_rsa_list = [0.3],dc_list = [0.5],t_resp_list=None,save_csv=True):

        if t_resp_list is None:
            t_resp_list = np.linspace(1/(2*self.f_base),4/self.f_base,t_step_size)
            # t_resp_list = np.linspace(self.f_base/2,8*self.f_base,f_step_size)

        '''
        time to frequency conversion
        2T0 -> T0/8
        f0/2 -> 8f0
        '''
        f_heart_max = find_firing_freq(self.i_base+self.rsa) / 10**3 

        # loop over duty cycles
        for dc in dc_list:
            self.dc = dc
            # loop over relative amplitudes (coupling strengths)
            for rel_rsa in rel_rsa_list:
                X1,X2,X3 = [],[],[]
                Y1,Y2,Y3 = [],[],[]
                print('Current rel_RSA:',rel_rsa)
                self.rsa = rel_rsa * self.i_base
                # vary period of stimulus current
                for t_resp in t_resp_list:
                    print('Current Tresp/T0:',t_resp*self.f_base)
                    self.t_resp = t_resp
                    '''
                    Frequency list!
                    '''
                    V,t = self.forward()
                    peak_indices,_ = find_peaks(V,height=0)
                    
                    assert len(peak_indices)!=0, '[ERROR] no peaks detected, check i_base and make sure it is set to a value greater than 2 nA'
                    assert np.amin(t) <= steady_state_time <= np.amax(t), 'Steady state time not within t grid, check time grid and steady state vals'

                    peak_times = t[peak_indices]
                    steady_state_peaks = peak_times[peak_times>=steady_state_time]
                    t_heart = (steady_state_peaks[-1] - steady_state_peaks[0]) / len(steady_state_peaks)
                    
                    X1.extend(peak_times)
                    X2.append(self.t_resp*self.f_base)
                    X3.append(2*np.pi/self.t_resp)
                    
                    
                    Y1.extend([self.t_resp/t_heart] * len(peak_times))
                    Y2.append(t_heart*self.f_base)
                    Y3.append(self.t_resp/t_heart)
                    
                
                pd.DataFrame(np.asarray([X1,Y1])).to_csv(f"raster_plot_META_i_base_{self.i_base}_rsa_{self.rsa}_t_heart_base_{1/self.f_base}_t_heart_max_{1/f_heart_max}_dc_{dc}_steady_state_t_{steady_state_time}_t_resp_min_{t_resp_list[0]*self.f_base}_t_resp_max_{t_resp_list[-1]*self.f_base}_t_max_{self.t[-1]}.csv",header=None, index= None)
            
                pd.DataFrame(np.asarray([X2,Y2])).to_csv(f"lung_heart_period_plot_META_i_base_{self.i_base}_rsa_{self.rsa}_t_heart_base_{1/self.f_base}_dc_{dc}_steady_state_t_{steady_state_time}_t_resp_min_{t_resp_list[0]*self.f_base}_t_resp_max_{t_resp_list[-1]*self.f_base}_t_max_{self.t[-1]}.csv",header=None, index= None)

                pd.DataFrame(np.asarray([X3,Y3])).to_csv(f"lung_heart_freq_ratio_plot_META_i_base_{self.i_base}_rsa_{self.rsa}_t_heart_base_{1/self.f_base}_dc_{dc}_steady_state_t_{steady_state_time}_t_resp_min_{t_resp_list[0]*self.f_base}_t_resp_max_{t_resp_list[-1]*self.f_base}_t_max_{self.t[-1]}.csv",header=None, index= None)


    def plot_from_csvs(self,csv_path_list,order_by_rsa = True, save_fig=True):
        '''
        takes in a list of csv files and plots them
        x and y labels set from fname 
        
        gray scale used to differentiate between different datasets 
        '''
        color_scale_list = np.zeros((len(csv_path_list)))
        ic(color_scale_list)
        type_str = csv_path_list[0].split('\\')[-1].split('META')[0][:-1]
        # meta_str = 
        # dc
        # steady_state_time
        if order_by_rsa:
            rsa_list = []
            for csv in csv_path_list:
                csv_type,csv_meta_str = csv.split('\\')[-1].split('META')
                rsa = float(csv_meta_str.split('rsa_')[1].split('_')[0])
                rsa_list.append(rsa)
            enumerated_csv_path_list = sorted(enumerate(csv_path_list),key = lambda x:rsa_list[x[0]],reverse=True)
            
            ic(rsa_list)

        fig,ax = plt.subplots(2,2)
        ax = ax.flatten()    
        for idx,csv in enumerated_csv_path_list:
            ax_=ax[idx]
            csv_type,csv_meta_str = csv.split('\\')[-1].split('META')
            ic(csv_meta_str)
            csv_type = csv_type[:-1]
            ic(csv_type)

            # assert csv_type[:-1] == type_str, '[ERROR] csv paths given contain data for different types of plots!'
            rsa = float(csv_meta_str.split('rsa_')[1].split('_')[0])
            X,Y = np.asarray(pd.read_csv(csv,header=None))
            t_heart_0 = float(csv_meta_str.split('t_heart_base_')[1].split('_')[0])
            w_heart_0 = 2*np.pi / t_heart_0
            rel_rsa =  np.round((rsa) / 2.32,2)  

            # plt.title(csv_meta_str)
            if csv_type == 'raster_plot': 
                plt.xlabel(r'$time (ms)$',fontsize=12)
                plt.ylabel(r'$T_{resp}$/$T_{heart}$',fontsize=12)
                                
            elif csv_type == 'lung_heart_period_plot':
                ax_.set_xlabel(r'$t_{resp}/t_{heart_{0}}$',fontsize=12)
                ax_.set_ylabel(r'$t_{heart}/t_{heart_{0}}$',fontsize=12)
                ax_.set_title(f'rel rsa = {rel_rsa}')
                ax_.set_ylim(0.4,1)
            else:
                X/=w_heart_0
                plt.xlabel(r'$ω_{lung}$/$ω_{heart_{0}}$',fontsize=12)
                plt.ylabel(r'$ω_{heart} / ω_{lung}$',fontsize=12)
            
            #plt.scatter(X,Y,marker='o',s=1,c=f'{color_scale_list[idx]}')
            ic(rsa)
            ax_.plot(X,Y,ms=1,c=f'{color_scale_list[-(idx+1)]}', label = f'rel rsa = {rel_rsa}')
            # plt.legend()
        plt.show()
        plt.clf()
        if save_fig: plt.savefig(f"{csv[:-4]}.png")

                
    

        
if __name__ == '__main__':
    ds = 0.5
    i_base = 2.32 # nA
    t_base = 10**3/find_firing_freq(i_base)
    t_resp = 20*t_base
    HHplots = HodgkinHuxleyModelPlots(t_resp=t_resp,i_base=2.32,rsa=1.2,dc = 0.5,t_max=500)
    csv_data = [f for f in os.listdir() if f.endswith('csv')]
    rel_rsa = [1.2,1.5]
    # HHplots.plot_from_csvs(csv_data)
    HHplots.gen_action_potential_plot()
    # HHplots.gen_arnold_tongue_plot_data(steady_state_time=0,rel_rsa_list = rel_rsa)

    # HHplots.gen_action_potential_plot()
    # plot t resp grid on arnold tongue raster plot
    # how does theart change wrt to t_resp
    # how does t_heart inhale change wrt t_heart base as a function of t_resp 
    # plot consecutive t_hearts and see how it changes during a mode locked period vs a normal period

    # identify mode locked regimes delta t curves vs unlocked and see how this varies as a function of n:m tongue and DC 
             




    
        
    
    
    