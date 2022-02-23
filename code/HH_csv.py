import numpy as np
import pandas as pd 
from scipy.signal import find_peaks
import os
import pylab as plt
class hh_processesor():
	'''
	Class to process the raw data produced by the HH C code

	Attributes:
	-----------
	None
	'''

	def __init__(self): 
		return


	def TI_hh_processor(self, bin_folder_name):
		'''
		Function that processes all raw time independent hh model CSV file data in a specific directory 
		and saves the processed data to a seperate directory for use

		Attributes:
		-----------
		bin_foler_name : String
			Name of the folder in the bin that stores the information-
			Format of folder name should be TD_DC_baseCurrent_relRSA

		'''

		#Gets the path that the data is stored in
		curr_dir = os.getcwd()
		bin_dir = os.path.join(curr_dir, 'Bin')
		working_dir = os.path.join(bin_dir, str(bin_folder_name))


		#Initialises a list that contains all the csv files of interest
		files = [f for f in os.listdir(working_dir) if f.startswith('sigma')]

		#Initialise a data frames to store the processed data
		stair_case_df = pd.DataFrame(columns = ['w_lung', 'w_heart'])
		rasta_df = pd.DataFrame	(columns = ['time_arr', 't_resp/t_heart_arr', 't_resp/t_heart'])
		lung_heart_df = pd.DataFrame(columns = ['t_resp', 't_heart'])


		for i, file in enumerate(files):
			print(file)
			#Gets path for file of interest
			file_path = os.path.join(working_dir, str(file))

			#Creates a date frame
			df = pd.read_csv(file_path, names = ['t','I','V','m','h','n']) 

			t= df.loc[:,'t'].values
			V =  df.loc[:,'V'].values
			print(type(t))
			#Finds t_resp (stored in file name)
			t_resp = float(file.split('_')[-2])
			# print(t_resp)

			#Gets index position of each peak
			peak_indices,_ = find_peaks(V,height=0) 

			# Gets average t_heart 
			peak_times = t[peak_indices]
			t_heart = (peak_times[-1] - peak_times[0]) / len(peak_times)

			#Appends a row of processed data to df
			stair_case_df.loc[i] = [2*np.pi/t_resp, t_resp/t_heart]
			rasta_df.loc[i] = [peak_times, (t_resp/t_heart)*np.ones(len(peak_times)), t_resp/t_heart]
			lung_heart_df.loc[i] = [t_resp,t_heart]

		#Sorts Dataframes
		stair_case_df.sort_values('w_lung')
		rasta_df.sort_values('t_resp/t_heart')
		lung_heart_df.sort_values('t_resp')


		#Gets directories to save the data in
		processed_data_dir = os.path.join(curr_dir, 'Processed_data')

		stair_case_dir = os.path.join(processed_data_dir, 'Stair_case_data')
		rasta_dir = os.path.join(processed_data_dir, 'Rasta_data')
		lung_heart_dir = os.path.join(processed_data_dir,'Lung_heart_data')

		#Gets path to save the CSVs
		stair_case_path = os.path.join(stair_case_dir, bin_folder_name) +'.csv'
		rasta_path = os.path.join(rasta_dir, bin_folder_name) + '.csv'
		lung_heart_path = os.path.join(lung_heart_dir, bin_folder_name) + '.csv'

		#Saves df as CSVs
		stair_case_df.to_csv(stair_case_path, index = False)
		rasta_df.to_csv	(rasta_path, index = False)
		lung_heart_df.to_csv(lung_heart_path, index = False)

	def process_triangles(self, filename):

		#Gets path data stored as a csv
		curr_dir = os.getcwd()
		processed_data_dir = os.path.join(curr_dir, 'Processed_data')
		stair_case_dir = os.path.join(processed_data_dir, 'Stair_case_data')
		path = os.path.join(stair_case_dir, str(filename))




class plotter():

	def plot_stair_case(self, filename):

		#Gets path data stored as a csv
		curr_dir = os.getcwd()
		processed_data_dir = os.path.join(curr_dir, 'Processed_data')
		stair_case_dir = os.path.join(processed_data_dir, 'Stair_case_data')
		path = os.path.join(stair_case_dir, str(filename))

		df = pd.read_csv(path)

		x = df.loc[:,'w_lung'].values
		y = df.loc[:, 'w_heart'].values
		# print(type(x[1]))
		plt.scatter(x,y)
		plt.show()
	



	def plot_rasta(self,filename):


		#Gets path data stored as a csv
		curr_dir = os.getcwd()
		processed_data_dir = os.path.join(curr_dir, 'Processed_data')
		stair_case_dir = os.path.join(processed_data_dir, 'Rasta_data')
		path = os.path.join(stair_case_dir, str(filename))

		df = pd.read_csv(path)
		# print(df)

		x_matrix = df.loc[:,'time_arr']
		y_matrix = df.loc[:,'t_resp/t_heart_arr']
		# x_data = df['time_arr'].to_numpy()
		# x_data

		# x_data = [[float(i) for i in e] for e in x_data[1:]]
		# print(type(x_data[0]))
		print(x_matrix[0])
		for i in range(len(x_matrix)):

			x_data = x_matrix[i].strip('[').strip(']').replace('\n',''  ).replace('  ', ' ').split(' ')
			y_data = y_matrix[i].strip('[').strip(']').replace('\n',''  ).replace('  ', ' ').split(' ')
			

			# split('[')[1].split(']')[0].strip('\n').split('  ')[1:]
			# y_data = y_matrix[i].split('[')[1].split(']')[0].strip('\n').split('  ')[1:]
			for x in x_data:
				print(x)
			print(len(x_data))
			print(len(y_data))
				
			x = [float(j) for j in x_data if j != '']
			y = [float(j) for j in y_data if j != '']

			plt.scatter(x,y)


		plt.show()






processor = hh_processesor()
# processor.TI_hh_processor('sigma_0.000000_tdeIdx_1_iBase_0.002400_relRsa_0.300000_dc_0.500000')
plot = plotter()
# plot.plot_rasta('sigma_0.000000_tdeIdx_1_iBase_0.002400_relRsa_0.300000_dc_0.500000.csv')
plot.plot_stair_case('sigma_0.000000_tdeIdx_1_iBase_0.002400_relRsa_0.300000_dc_0.500000.csv')
