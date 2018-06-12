import numpy as np
import math
from scipy.linalg import logm, expm
import re
import json
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import gc
import math as m
import cmath
import _pickle as pickle

class MM():
	
	
	

	def __init__(self):
		self.data={}
		self.muellerNames =['mm11', 'mm12', 'mm13', 'mm14',
			'mm21', 'mm22', 'mm23', 'mm24',
			'mm31', 'mm32', 'mm33', 'mm34',
			'mm41', 'mm42', 'mm43', 'mm44']
		



#*****************************************************************************************************************#		
	def ReadData(self,directory,filename,wvlcutofmin,wvlcutofmax):
		'''
		Reads data from a single file with the filename from the given directory, in range between wvlcutofmin,wvlcutofmax.
		'''
		
		print ("Processing File")
		if filename+'.p' in os.listdir(directory):
			with open(directory+'\\'+filename+'.p','rb') as f:
				self.data=pickle.load(f)
				print('Loaded from pickle')
				print(type(self.data))
		else:		
			with open(directory+'/'+filename+ '.dat') as f:
				for lines in f:				
					if lines[0]=='m':								   # skip first non mm-values				 
						lines = re.split(r'\t+', lines)					   #formating the lines
					
						AOI=lines[2]
						if AOI not in self.data.keys():
							self.data[AOI]=dict()
						
						try:		
							AZI=float(lines[7])
						except IndexError:
							AZI=0.0
						finally:
							if AZI not in self.data[AOI].keys():
								self.data[AOI][AZI]=dict() 
								
								
						wvl=float(lines[1])									 
						if wvl>=wvlcutofmin and wvl<=wvlcutofmax and wvl not in self.data[AOI][AZI].keys():
							self.data[AOI][AZI][wvl]=list()
						
						
						
						if lines[0] in self.muellerNames and wvl<=wvlcutofmax and wvl>=wvlcutofmin:
							if lines[3]=='Infinity':
								self.data[AOI][AZI][wvl].append(float(1))
							else:						 
								self.data[AOI][AZI][wvl].append(float(lines[3]))
						
						
			os.chdir(directory)
			pickle.dump( self.data, open( filename+".p", "wb" ) )	
		for AOI in sorted(self.data.keys()):
			for azi in sorted(self.data[AOI].keys()):
				for wl in sorted(self.data[AOI][azi].keys()):
					self.data[AOI][azi][wl] = np.reshape(np.asarray(self.data[AOI][azi][wl]), (4,4))
		
		
						
		
		
#*****************************************************************************************************************#
	def ADDA_read(self,directory,Norm):
		'''
		Reads data from a single file with the filename from the given directory, in range between wvlcutofmin,wvlcutofmax.
		
		Provides normalized data if  Norm is True
		'''
		
		print ("Processing File")
		if 'result'+'.p' in os.listdir(directory):
			with open(directory+'\\'+'result.p','rb') as f:
				self.data=pickle.load(f)
				print('Loaded from pickle')
				print(type(self.data))
		else:
			folders=os.listdir(directory)
			#print(folders)
			for folder in folders:
				print(folder)
				with open(directory+'/'+folder+'/'+'mueller') as f:
					for lines in f:
						if lines.startswith('theta'):
							continue
						else:				 
							lines = re.split(r'\ +', lines)					   #formating the lines
					
							AOI=str(0.000)
							if AOI not in self.data.keys():
								self.data[AOI]=dict()
						
								
							AZI=float(float(lines[0]))
							if AZI not in self.data[AOI].keys():
								self.data[AOI][AZI]=dict()
						
								
								
							wvl=float(folder)
							if wvl not in self.data[AOI][AZI].keys():
								self.data[AOI][AZI][wvl]=list()
						
						
						
							for i in range(1,17):
								self.data[AOI][AZI][wvl].append(float(lines[i]))
						
						
				os.chdir(directory)
				pickle.dump( self.data, open( "result.p", "wb" ) )	
		
		
		for AOI in sorted(self.data.keys()):
			for azi in sorted(self.data[AOI].keys()):
				for wl in sorted(self.data[AOI][azi].keys()):
					if Norm==True:
						m11=self.data[AOI][azi][wl][0]
						for i in range(0,16):
							self.data[AOI][azi][wl][i]=self.data[AOI][azi][wl][i]/m11
					self.data[AOI][azi][wl] = np.reshape(np.asarray(self.data[AOI][azi][wl]), (4,4))
		
		print(self.data)
						
		
		
#*****************************************************************************************************************#
	def Plot(self,what,angle,title,mult):
		'''
		Plot data specified as what parameter. The title is the prefix (type of data) of the titel.
		'''
		if what =='data':
			datas=self.data
		elif what =='Lm':
			datas=self.Lm
		elif what =='Lu':
			datas=self.Lu	
		# generate dictionary for each AOI to plot
		muellerGraphs = dict()
		for AOI in datas.keys():
			muellerGraphs[AOI] = dict()
	
		# generate dictionary entry for each mm-element
		for AOI in muellerGraphs.keys():
			for element in self.muellerNames:
				muellerGraphs[AOI][element] = list()
	
		# now insert all mm-values of each element
		for AOI in sorted(datas.keys()):
			for AZI in sorted(datas[AOI].keys()):
				for wvl in sorted(datas[AOI][AZI].keys()):
					currMM = datas[AOI][AZI][wvl].reshape(-1)
					for element, i in zip(sorted(self.muellerNames), range(16)):
						muellerGraphs[AOI][element].append(currMM[i])
	
		AOIkeys = list(datas.keys())
		#print(AOIkeys)
		aziKeys = list(datas[AOIkeys[0]].keys())
		
		#print(aziKeys)
		
		azimuth = np.radians(sorted(np.asarray(list(datas[AOIkeys[0]].keys()))))
		zenith	= np.asarray(sorted(np.asarray(list(datas[AOIkeys[0]][aziKeys[0]].keys()))))
		
		#print(zenith)
		# adjust plot settings
		mpl.rcParams['axes.labelsize']	= 6
		mpl.rcParams['xtick.labelsize'] = 10
		mpl.rcParams['ytick.labelsize'] = 15
		mpl.rcParams['legend.fontsize'] = 15
		
		datax4=['mm12','mm21','mm13','mm31','mm24','mm42']
		datax10=['mm32','mm23',]
		datax100=['mm14','mm41']
		diag=['mm11','mm22','mm33','mm44']
		
		#saveDirectory = 'D:/Ievgen Ell/Measurements/NanoDiscs/MM_transmition/OD5/result_MM_plot/'
		#if not os.path.isdir(saveDirectory):
		#os.makedirs(saveDirectory)
		print('Building plot')
		for AOI in sorted(muellerGraphs.keys()):
			if AOI==angle:
				fig, axs = plt.subplots(4, 4, subplot_kw=dict(projection='polar'))	
				for element, ax in zip(sorted(muellerGraphs[AOI].keys()), iter(axs.reshape(-1))):
					data = np.reshape(muellerGraphs[AOI][element], (azimuth.size, zenith.size)).T
					if element in datax4:
						data=data*4*mult
						ax.set_title('x{}'.format(mult*4))
					elif element in datax10:
						data=data*10*mult
						ax.set_title('x{}'.format(mult*10))
					elif element in datax100:
						data=data*100*mult
						ax.set_title('x{}'.format(mult*100))
					#elif element in diag:
						#data=data*10
						#ax.set_title('x10')
						
					
					cs = ax.contourf(azimuth, zenith, data, 50, vmin=-1,vmax=1,cmap='jet')
					ax.grid(False)
					#ax.set_title('%s'%(element),fontsize=8, horizontalalignment='left', verticalalignment='bottom')		   
					ax.yaxis.set_major_formatter(plt.NullFormatter())
					thetaticks = np.arange(0,360,180)
					ax.set_thetagrids(thetaticks, frac=1.3)
					minmax = np.asarray([np.min(data), np.max(data)])
				
				#plt.tight_layout()
				#print('Saving L_u %s plot at %s'%(element, AOI))
				#plt.savefig(saveDirectory + 'L_u_AOI_%03.1f.png'%float(AOI), dpi=600, bbox_inches='tight')
				cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
				cNorm = mpl.colors.Normalize(vmin=-1, vmax=1)
				cb1 = mpl.colorbar.ColorbarBase(cax, norm=cNorm,cmap='jet')
				plt.suptitle("%s plot at AOI-%s from %s to %s wvl."%(title,AOI,min(zenith),max(zenith)),fontsize=18)
				plt.show()
				
			else:
				continue
		print('Finished!!!!!!!!!')

#*****************************************************************************************************************#



	def Cloude(self,AOI,AZI,Temp):
		'''
		Performs a Cloude decomposition. The outcome is the plot of eigenvalues at given AOI and AZI.
		Returns the plot  of eigenvlues versus wavelength.
		'''
		m1=[]
		m2=[]
		m3=[]
		m4=[]
		ev1=[]
		ev2=[]
		ev3=[]
		ev4=[]
		
		for wvl in self.data[AOI][AZI].keys():
			M = np.copy(self.data[AOI][AZI][wvl])
			# calculate coherency matrix H of Mueller Matrix M
			M[0][0] = 1.0
		
			H = np.reshape(np.zeros(16, dtype='complex64'), (4, 4))
			H[0][0] = M[0][0]	 + M[1][1]	  + M[2][2]	   + M[3][3]
			H[0][1] = M[0][1]	 + M[1][0]	  - 1j*M[2][3] + 1j*M[3][2]
			H[0][2] = M[0][2]	 + M[2][0]	  + 1j*M[1][3] - 1j*M[3][1]	 
			H[0][3] = M[0][3]	 - 1j*M[1][2] + 1j*M[2][1] + M[3][0]  
			
			H[1][0] = M[0][1]	 + M[1][0]	  + 1j*M[2][3] - 1j*M[3][2]
			H[1][1] = M[0][0]	 + M[1][1]	  - M[2][2]	   - M[3][3]
			H[1][2] = 1j*M[0][3] + M[1][2]	  + M[2][1]	   - 1j*M[3][0]
			H[1][3] = 1j*M[2][0] - 1j*M[0][2] + M[1][3]	   + M[3][1]
			
			H[2][0] = M[0][2]	 + M[2][0]	  - 1j*M[1][3] + 1j*M[3][1]
			H[2][1] = M[1][2]	 - 1j*M[0][3] + M[2][1]	   + 1j*M[3][0]
			H[2][2] = M[0][0]	 - M[1][1]	  + M[2][2]	   -  M[3][3]
			H[2][3] = 1j*M[0][1] - 1j*M[1][0] + M[2][3]	   + M[3][2]
			
			H[3][0] = M[0][3]	 + 1j*M[1][2] - 1j*M[2][1] + M[3][0]
			H[3][1] = 1j*M[0][2] - 1j*M[2][0] + M[1][3]	   + M[3][1]
			H[3][2] = 1j*M[1][0] - 1j*M[0][1] + M[2][3]	   + M[3][2]
			H[3][3] = M[0][0]	 - M[1][1]	  - M[2][2]	   + M[3][3]
		
			H = np.divide(H, 4.0)
		
			eigValH, eigVecH = np.linalg.eig(H)
		
			# Eigenvectors of H must be real (H is hermitian).
			# Therefore we neglect the imaginary parts, since
			# they come from numerical uncertainties.
			# Real part of eigValH is in the range 1E-1 ... 1E-2
			# Imaginary part of eigValH is in the range 1E-29 ... 1E-31 
			eigValH = np.real(eigValH)
			eigVecH = np.real(eigVecH)
			
			# sort Eigenvalues and Eigenvectors descending by Eigenvalues
			idx = eigValH.argsort()[::-1]
			eigValH = eigValH[idx]
			eigVecH = eigVecH[:,idx]
			
			#------------------------------------------------------------------------
			# calculate Jones matrices (M_Ji, i=0...3) from the obtained Eigenvectors
			#------------------------------------------------------------------------
			
			# create empty jones and mueller matrices
			jonesMatrix = dict()
			cloudeDecom = dict()
			
			JtoMtransformer = np.array([[1.0, 0.0, 0.0,	 1.0],
										[1.0, 0.0, 0.0, -1.0],
										[0.0, 1.0, 1.0,	 0.0],
										[0.0, 1j , -1j,	 0.0]])
			invJtoMtransformer = np.array([[ 0.5,  0.5,	 0.0,  0.0	 ],
										[ 0.0,	0.0,  0.5, -0.5*1j],
										[ 0.0,	0.0,  0.5,	0.5*1j],
										[ 0.5, -0.5,  0.0,	0.0	  ]])
			
			for i in range(4):
				jonesMatrix[i] = np.reshape(np.zeros(4, dtype='complex64'), (2,2))
				cloudeDecom[i] = np.reshape(np.zeros(16, dtype='float32'), (4,4))
			
			for i in range(4):
				# assign jones matrix values from eigenvectors
				jonesMatrix[i][0][0] = eigVecH[i][0] + eigVecH[i][1]		# = r_pp
				jonesMatrix[i][0][1] = eigVecH[i][2] - 1j*eigVecH[i][3]		# = r_ps
				jonesMatrix[i][1][0] = eigVecH[i][2] + 1j*eigVecH[i][3]		# = r_sp
				jonesMatrix[i][1][1] = eigVecH[i][0] - eigVecH[i][1]		# = r_ss
		
				# make mueller matrix from jones matrix
				J = np.kron(jonesMatrix[i], np.conjugate(jonesMatrix[i]))
				J = np.dot(J, invJtoMtransformer)
				MJ = np.dot(JtoMtransformer, J)
				# Mueller matrix must be real, so remove the imaginary parts
				# since they are here always zero.
				cloudeDecom[i] = np.copy(np.real(MJ))
			
			# At this point we got the complete decomposited mueller matrix.
			# Now we delete all matrices which have a negative eigenvalue from
			# above as pre-factor and after that we divide each remaining matrix
			# by the sum of its eigenvalues.
			# Finally, summing up all remaining matrices delivers a Mueller matrix,
			# which is free from non-physical realizable parts and which is the
			# closest approximation to the original measured Mueller matrix.
			
			# This will be the final mueller matrix to be returned
			# with removed non-physical realizable part
			MM = np.reshape(np.zeros(16, dtype='float32'), (4,4))
			
			#for i in sorted(cloudeDecom.keys()):
				## divide matrix by the sum of its eigenvalues
				## if the decomposition factor obtained from H is > 0.0
				#if eigValH[i] > 0.0:
					##eValCD, eVecCD = np.linalg.eig(cloudeDecom[i])
					##eValCD = np.real(eValCD) # just to be sure it will be real ;-)
					##cloudeDecom[i] = np.multiply(cloudeDecom[i], (eigValH[i]/np.sum(eValCD)))
					#cloudeDecom[i] = eigValH[i] * cloudeDecom[i]
					#MM = np.add(MM, cloudeDecom[i])
					
			#eigValMM, eigVecMM = np.linalg.eig(MM)
			#MM /= np.sum(eigValMM)
			
			MM1 = eigValH[0]*cloudeDecom[0]
			MM2 = eigValH[1]*cloudeDecom[1]
			MM3 = eigValH[2]*cloudeDecom[2]
			MM4 = eigValH[3]*cloudeDecom[3]
			m1.append(MM1)
			m2.append(MM2)
			m3.append(MM3)
			m4.append(MM4)
			
			ev1.append(eigValH[0])
			ev2.append(eigValH[1])
			ev3.append(eigValH[2])
			ev4.append(eigValH[3])
			
			
		plt.plot(sorted(list(self.data[AOI][AZI].keys())),ev1,label='Lambda1')
		plt.plot(sorted(list(self.data[AOI][AZI].keys())),ev2,label='Lambda2')
		plt.plot(sorted(list(self.data[AOI][AZI].keys())),ev3,label='Lambda3')
		plt.plot(sorted(list(self.data[AOI][AZI].keys())),ev4,label='Lambda4')
		plt.suptitle('Eigenvalues of Cloude decomposition at AOI-{}, AZI-{}'.format(AOI,AZI))
		plt.xlim([min(self.data[AOI][AZI].keys()),max(self.data[AOI][AZI].keys())])
		plt.legend()
		plt.show()
	
#*****************************************************************************************************************#		
	def TdepCloude(self,AOI,AZI,Temp):
		'''
		Performs a Cloude decomposition on the temperature dependent dataset. The outcome is the plot of eigenvalues at given AOI and AZI.
		Returns the plot  of eigenvlues versus wavelength together with the eigenvalues of the hermitian matrix.
		'''
		m1=[]
		m2=[]
		m3=[]
		m4=[]
		ev1=[]
		ev2=[]
		ev3=[]
		ev4=[]
		EV={}
		for wvl in self.data[Temp][AOI][AZI].keys():
			M = np.copy(self.data[Temp][AOI][AZI][wvl])
			# calculate coherency matrix H of Mueller Matrix M
			M[0][0] = 1.0
		
			H = np.reshape(np.zeros(16, dtype='complex64'), (4, 4))
			H[0][0] = M[0][0]	 + M[1][1]	  + M[2][2]	   + M[3][3]
			H[0][1] = M[0][1]	 + M[1][0]	  - 1j*M[2][3] + 1j*M[3][2]
			H[0][2] = M[0][2]	 + M[2][0]	  + 1j*M[1][3] - 1j*M[3][1]	 
			H[0][3] = M[0][3]	 - 1j*M[1][2] + 1j*M[2][1] + M[3][0]  
			
			H[1][0] = M[0][1]	 + M[1][0]	  + 1j*M[2][3] - 1j*M[3][2]
			H[1][1] = M[0][0]	 + M[1][1]	  - M[2][2]	   - M[3][3]
			H[1][2] = 1j*M[0][3] + M[1][2]	  + M[2][1]	   - 1j*M[3][0]
			H[1][3] = 1j*M[2][0] - 1j*M[0][2] + M[1][3]	   + M[3][1]
			
			H[2][0] = M[0][2]	 + M[2][0]	  - 1j*M[1][3] + 1j*M[3][1]
			H[2][1] = M[1][2]	 - 1j*M[0][3] + M[2][1]	   + 1j*M[3][0]
			H[2][2] = M[0][0]	 - M[1][1]	  + M[2][2]	   -  M[3][3]
			H[2][3] = 1j*M[0][1] - 1j*M[1][0] + M[2][3]	   + M[3][2]
			
			H[3][0] = M[0][3]	 + 1j*M[1][2] - 1j*M[2][1] + M[3][0]
			H[3][1] = 1j*M[0][2] - 1j*M[2][0] + M[1][3]	   + M[3][1]
			H[3][2] = 1j*M[1][0] - 1j*M[0][1] + M[2][3]	   + M[3][2]
			H[3][3] = M[0][0]	 - M[1][1]	  - M[2][2]	   + M[3][3]
		
			H = np.divide(H, 4.0)
		
			eigValH, eigVecH = np.linalg.eig(H)
		
			# Eigenvectors of H must be real (H is hermitian).
			# Therefore we neglect the imaginary parts, since
			# they come from numerical uncertainties.
			# Real part of eigValH is in the range 1E-1 ... 1E-2
			# Imaginary part of eigValH is in the range 1E-29 ... 1E-31 
			eigValH = np.real(eigValH)
			eigVecH = np.real(eigVecH)
			
			# sort Eigenvalues and Eigenvectors descending by Eigenvalues
			idx = eigValH.argsort()[::-1]
			eigValH = eigValH[idx]
			eigVecH = eigVecH[:,idx]
			EV[wvl]=eigValH
			#------------------------------------------------------------------------
			# calculate Jones matrices (M_Ji, i=0...3) from the obtained Eigenvectors
			#------------------------------------------------------------------------
			
			# create empty jones and mueller matrices
			jonesMatrix = dict()
			cloudeDecom = dict()
			
			JtoMtransformer = np.array([[1.0, 0.0, 0.0,	 1.0],
										[1.0, 0.0, 0.0, -1.0],
										[0.0, 1.0, 1.0,	 0.0],
										[0.0, 1j , -1j,	 0.0]])
			invJtoMtransformer = np.array([[ 0.5,  0.5,	 0.0,  0.0	 ],
										[ 0.0,	0.0,  0.5, -0.5*1j],
										[ 0.0,	0.0,  0.5,	0.5*1j],
										[ 0.5, -0.5,  0.0,	0.0	  ]])
			
			for i in range(4):
				jonesMatrix[i] = np.reshape(np.zeros(4, dtype='complex64'), (2,2))
				cloudeDecom[i] = np.reshape(np.zeros(16, dtype='float32'), (4,4))
			
			for i in range(4):
				# assign jones matrix values from eigenvectors
				jonesMatrix[i][0][0] = eigVecH[i][0] + eigVecH[i][1]		# = r_pp
				jonesMatrix[i][0][1] = eigVecH[i][2] - 1j*eigVecH[i][3]		# = r_ps
				jonesMatrix[i][1][0] = eigVecH[i][2] + 1j*eigVecH[i][3]		# = r_sp
				jonesMatrix[i][1][1] = eigVecH[i][0] - eigVecH[i][1]		# = r_ss
		
				# make mueller matrix from jones matrix
				J = np.kron(jonesMatrix[i], np.conjugate(jonesMatrix[i]))
				J = np.dot(J, invJtoMtransformer)
				MJ = np.dot(JtoMtransformer, J)
				# Mueller matrix must be real, so remove the imaginary parts
				# since they are here always zero.
				cloudeDecom[i] = np.copy(np.real(MJ))
			
			# At this point we got the complete decomposited mueller matrix.
			# Now we delete all matrices which have a negative eigenvalue from
			# above as pre-factor and after that we divide each remaining matrix
			# by the sum of its eigenvalues.
			# Finally, summing up all remaining matrices delivers a Mueller matrix,
			# which is free from non-physical realizable parts and which is the
			# closest approximation to the original measured Mueller matrix.
			
			# This will be the final mueller matrix to be returned
			# with removed non-physical realizable part
			MM = np.reshape(np.zeros(16, dtype='float32'), (4,4))
			
			#for i in sorted(cloudeDecom.keys()):
				## divide matrix by the sum of its eigenvalues
				## if the decomposition factor obtained from H is > 0.0
				#if eigValH[i] > 0.0:
					##eValCD, eVecCD = np.linalg.eig(cloudeDecom[i])
					##eValCD = np.real(eValCD) # just to be sure it will be real ;-)
					##cloudeDecom[i] = np.multiply(cloudeDecom[i], (eigValH[i]/np.sum(eValCD)))
					#cloudeDecom[i] = eigValH[i] * cloudeDecom[i]
					#MM = np.add(MM, cloudeDecom[i])
					
			#eigValMM, eigVecMM = np.linalg.eig(MM)
			#MM /= np.sum(eigValMM)
			
			MM1 = eigValH[0]*cloudeDecom[0]
			MM2 = eigValH[1]*cloudeDecom[1]
			MM3 = eigValH[2]*cloudeDecom[2]
			MM4 = eigValH[3]*cloudeDecom[3]
			m1.append(MM1)
			m2.append(MM2)
			m3.append(MM3)
			m4.append(MM4)
			
			ev1.append(eigValH[0])
			ev2.append(eigValH[1])
			ev3.append(eigValH[2])
			ev4.append(eigValH[3])
			
			
		#plt.plot(sorted(list(self.data[Temp][AOI][AZI].keys())),ev1,label='Lambda1')
		#plt.plot(sorted(list(self.data[Temp][AOI][AZI].keys())),ev2,label='Lambda2')
		#plt.plot(sorted(list(self.data[Temp][AOI][AZI].keys())),ev3,label='Lambda3')
		#plt.plot(sorted(list(self.data[Temp][AOI][AZI].keys())),ev4,label='Lambda4')
		#plt.suptitle('Eigenvalues of Cloude decomposition at AOI-{}, AZI-{}'.format(AOI,AZI))
		#plt.xlim([min(self.data[Temp][AOI][AZI].keys()),max(self.data[Temp][AOI][AZI].keys())])
		#plt.legend()
		#plt.show()
		return EV
#*****************************************************************************************************************#	

	def diffDecomp(self):
		'''
		Differential decomposition is performed on a given data. The Lm and Lu are stored in a self.Lm and self.Lu
		'''
		print("Differential decomposition")
		g=[1,0,0,0,
		0,-1,0,0,
		0,0,-1,0,
		0,0,0,-1]
		G=np.reshape(g,(4,4))
		L=dict()
		L_u=dict()
		L_m=dict()
		
		
		for AOI in sorted(self.data.keys()):
			if AOI not in L.keys():
				L[AOI]=dict()
				L_m[AOI]=dict()
				L_u[AOI]=dict()
			for AZI in sorted(self.data[AOI].keys()):
				if AOI not in L[AOI].keys():
					L[AOI][AZI]=dict()
					L_m[AOI][AZI]=dict()
					L_u[AOI][AZI]=dict()			
				for wvl in sorted(self.data[AOI][AZI].keys()):
					if wvl not in L[AOI][AZI].keys():
						L[AOI][AZI][wvl]=np.empty([4,4],dtype=complex)
						L_m[AOI][AZI][wvl]=np.empty([4,4],dtype=complex)
						L_u[AOI][AZI][wvl]=np.empty([4,4],dtype=complex)
					if AOI=='0.000' or AOI=='0.0':
						mm=np.copy(self.data[AOI][AZI][wvl])
						mm[0][0]=1.0
							
						L[AOI][AZI][wvl] = logm(mm)	   # Right option 
						L_m[AOI][AZI][wvl]=1/2*(L[AOI][AZI][wvl]-(G@np.transpose(L[AOI][AZI][wvl]))@G)
						L_u[AOI][AZI][wvl]=1/2*(L[AOI][AZI][wvl]+(G@np.transpose(L[AOI][AZI][wvl]))@G)
					else:
						print('Differential decomposition is only valid for 0 AOI')
		
				
		self.Lm=L_m
		self.Lu=L_u
		
		
#*****************************************************************************************************************#
	def slicePlot(self,what,AOI, AZI,Temp,wvlmin,wvlmax,elements,normfac):
		'''
		Plots the linecut of chosen elements of type what, at AOI and AZI from wvlmin to wvlmax.
		If the elements is empty, plots all elements.
		Returns a 2D plot of the required values versus wavelength.
		'''
		plt.rcParams['axes.facecolor'] = 'white'
		s=0
		if Temp==[]:
			print('Building the plot')
			if what =='data':
				datas=self.data
			elif what =='Lm':
				datas=self.Lm
			elif what =='Lu':
				datas=self.Lu
			
			mpl.rcParams['axes.labelsize']	= 6
			mpl.rcParams['xtick.labelsize'] = 15
			mpl.rcParams['ytick.labelsize'] = 15
			mpl.rcParams['legend.fontsize'] = 15
			fig=plt.figure()
			colors=['k','r','brown','darkred','orange','darkorange','y','darksalmon','g','darkseagreen','c','darkcyan','b','darkslateblue','violet','darkviolet']
			
			for Azi in AZI:
				for aoi in sorted(AOI):
					if elements==[]:
						for i in range(0,16):
							plt.plot(sorted(list(datas[aoi][Azi].keys())),[datas[aoi][Azi][wvl].flatten()[i]/normfac for wvl in sorted(datas[aoi][Azi].keys())],c=colors[i],label='{}'.format(self.muellerNames[i]))
					else:
						for i in elements:
							plt.plot(sorted(list(datas[aoi][Azi].keys())),[datas[aoi][Azi][wvl].flatten()[i]/normfac for wvl in sorted(datas[aoi][Azi].keys())],c=colors[i],label='{}'.format(self.muellerNames[i]))
						
							
				plt.title('%s elements slice at AOI-%s , Azi-%s'%(what,AOI,Azi),fontsize=24)			
				plt.xlim([wvlmin,wvlmax])
				#plt.ylim([-0.1,0.1])
				plt.legend()
				plt.xlabel('Wavelength (nm)',fontsize=18)
				#plt.show()
		
		else:
			
			print('Building the plot')
			if what =='data':
				datas=self.data
			elif what =='Lm':
				datas=self.Lm
			elif what =='Lu':
				datas=self.Lu
			
			mpl.rcParams['axes.labelsize']	= 6
			mpl.rcParams['xtick.labelsize'] = 15
			mpl.rcParams['ytick.labelsize'] = 15
			mpl.rcParams['legend.fontsize'] = 15
			fig=plt.figure()
			colors=['k','r','brown','darkred','orange','darkorange','y','darksalmon','g','darkseagreen','c','darkcyan','b','darkslateblue','violet','darkviolet']
			
			for T in Temp:
				for Azi in AZI:
					for aoi in sorted(AOI):
						if elements==[]:
							for i in range(0,16):
								plt.plot(sorted(list(datas[T][aoi][Azi].keys())),[datas[T][aoi][Azi][wvl].flatten()[i]/normfac for wvl in sorted(datas[T][aoi][Azi].keys())],c=colors[i],label='{} at {}C'.format(self.muellerNames[i],T))
						else:
							for i in elements:
								WVL=sorted(list(datas[T][aoi][Azi].keys()))
								data=[datas[T][aoi][Azi][wvl].flatten()[i]/normfac for wvl in WVL]
								plt.plot(WVL,data,label='{} at {}C'.format(self.muellerNames[i],T))
						
							
				plt.title('%s elements slice at AOI-%s , Azi-%s'%(what,AOI,Azi),fontsize=24)			
				plt.xlim([wvlmin,wvlmax])
				#plt.ylim([-0.1,0.1])
				plt.legend()
				plt.xlabel('Wavelength (nm)',fontsize=18)
				#plt.show()
		


#*****************************************************************************************************************#	
	def depIndex(self,AOI,AZI):
		'''
		Calculate depolarization index from a MM at AOI and AZI.
		Returns a 2D plot of the depolarization index versus WVL.
		'''
		DI=[]
		for wvl in self.data[AOI][AZI].keys():				
			DI.append(np.sqrt((np.trace(np.dot(self.data[AOI][AZI][wvl].T,self.data[AOI][AZI][wvl]))-(self.data[AOI][AZI][wvl][0][0])**2)/3))
		plt. plot(sorted(self.data[AOI][AZI].keys()),DI)
		plt.title('Depolarization index at AOI-%s , Azi-%s'%(AOI,AZI),fontsize=24)			  
		plt.xlim([min(self.data[AOI][AZI].keys()),max(self.data[AOI][AZI].keys())])
		#plt.ylim([-0.1,0.1])
		plt.legend()
		plt.xlabel('Wavelength (nm)',fontsize=18)
		plt.show()
#*****************************************************************************************************************#	
	def polConvRead(self,directory, fileName,minWvl):
		'''
		Read polarization conversion data Psi_ps and Psi_sp from the dat file.
		Define directory and filename of the file. Also the minimal wavelength to read.
		Save data into self.ps and self.sp
		'''
		print('Reading polarization conversion data')
		Aps={}
		Asp={}
		with open(directory+"/"+fileName+'.dat','r') as f:
			for line in f:
				line=f.readline()
				line=re.split(r'\t+',line)
				
		
			
				if line[0]=='Aps' and float(line[1])>=minWvl:
					
					AOI=float(line[2])
					if AOI not in sorted(Aps.keys()):
						Aps[AOI]={}
				
					AZI=float(line[7])
					if AZI not in sorted(Aps[AOI].keys()):
						Aps[AOI][AZI]={}
					
					wvl=float(line[1])
					if wvl not in sorted(Aps[AOI][AZI].keys()):
						Aps[AOI][AZI][wvl]=[]
				
					Aps[AOI][AZI][wvl].append(float(line[3]))
				
			
				if line[0]=='Asp' and float(line[1])>=minWvl:
					
					AOI=float(line[2])
					if AOI not in sorted(Asp.keys()):
						Asp[AOI]={}
				
					AZI=float(line[7])
					if AZI not in sorted(Asp[AOI].keys()):
						Asp[AOI][AZI]={}
					
					wvl=float(line[1])
					if wvl not in sorted(Asp[AOI][AZI].keys()):
						Asp[AOI][AZI][wvl]=[]
				
					Asp[AOI][AZI][wvl].append(float(line[3]))
					
				
		self.Aps=Aps
		self.Asp=Asp
#*****************************************************************************************************************#	


	def plotPolConversion(self,type):
		'''
		Plots polarization conversion data from self.ps or self.sp
		'''

		
		if type=='sp':
			# Plot Psi_sp
			for AOI in sorted(self.Asp.keys()):
				data=[]
				for AZI in sorted(self.Asp[AOI].keys()):
					for wvl in sorted(self.Asp[AOI][AZI].keys()):
						data.append(self.Asp[AOI][AZI][wvl])
						
						
				
				AZIm=np.asarray(list(sorted(self.Asp[AOI].keys())))
				
				wvls=np.asarray(list(sorted(self.Asp[AOI][AZI].keys())))
				
				data=np.array(data)
				
				data=data.reshape(len(data)/len(wvls),len(wvls))
			
				
				plt.figure()
				cp = plt.contourf(wvls, AZIm, data,200,cmap='jet')
				plt.colorbar(cp,cmap='jet')
				axes = plt.gca()
				axes.set_ylim([0,360])
				#axes.set_xlim([min(self.data['0.000'][0].keys()),max(self.data['0.000'][0].keys())])
				plt.gca()
				plt.title(r'$\mathrm{\Psi_{sp}}$ at AOI-%s'%AOI,size=30)
				plt.xlabel('Wavelength (nm)',size=20)
				plt.ylabel('Azimuthal (deg)',size=20)
				plt.show()
			
			
		#Plot Psi_ps
		if type=='ps':
			# Plot Psi_ps
			for AOI in sorted(self.Aps.keys()):
				data=[]
				for AZI in sorted(self.Aps[AOI].keys()):
					for wvl in sorted(self.Aps[AOI][AZI].keys()):
						data.append(self.Aps[AOI][AZI][wvl])
						
						
				
				AZIm=np.asarray(list(sorted(self.Aps[AOI].keys())))
				
				wvls=np.asarray(list(sorted(self.Aps[AOI][AZI].keys())))
				
				data=np.array(data)
				
				data=data.reshape(len(data)/len(wvls),len(wvls))
			
				
				plt.figure()
				cp = plt.contourf(wvls, AZIm, data,200,cmap='jet')
				plt.colorbar(cp,cmap='jet')
				axes = plt.gca()
				axes.set_ylim([0,360])
				#axes.set_xlim([min(self.data['0.000'][0].keys()),max(self.data['0.000'][0].keys())])
				plt.gca()
				plt.title(r'$\mathrm{\Psi_{sp}}$ at AOI-%s'%AOI,size=30)
				plt.xlabel('Wavelength (nm)',size=20)
				plt.ylabel('Azimuthal (deg)',size=20)
				plt.show()
#******************************************************************************************************************************************#


	def multFileRead(self,angles,directory,filename,wvlcutofmin,wvlcutofmax):
		'''
		Read Mueller matrix data from the collection of files stored in the same folder.
		Each file must represent a separate azimuthal angle and the name of the file ends with the azimuthal angle.
		Saves Mueller matrix data into self.data
		'''
		for i in range(int((angles[1]-angles[0])/angles[2]+1)):		#Azimuthal angles are determined from file names
			print ("processing Azi %s"%(angles[0]+angles[2]*i)) 
			with open(directory+'/'+filename+str(angles[0]+angles[2]*i)+'.dat','r') as f:	 # Azimuthal angles are determined from file names
				
				for lines in f:
				
					if lines[0]=='m':									# skip first non mm-values
								
						lines = re.split(r'\t+', lines)					 
					
						AOI=lines[2]
						if AOI not in self.data.keys():
							self.data[AOI]=dict()
							
							
						
						AZI=float(angles[0]+angles[2]*i)
						if AZI not in self.data[AOI].keys():
							self.data[AOI][AZI]=dict()
							
							
						wvl=float(lines[1])									 
						if wvl>=wvlcutofmin and wvl<=wvlcutofmax and wvl not in self.data[AOI][AZI].keys():
							self.data[AOI][AZI][wvl]=list()
							
							
				
				
						
						if lines[0] in self.muellerNames and wvl<=wvlcutofmax and wvl>=wvlcutofmin:
							#if wvl==wvlcutofmax:
								#muellerMatrices[AOI][AZI][wvl].append(0.5)
							#elif wvl==wvlcutofmin:
								# muellerMatrices[AOI][AZI][wvl].append(-0.5) 
							#else:	  
							self.data[AOI][AZI][wvl].append(float(lines[3])) 
					
					
						
					
		#for AOI in sorted(muellerMatrices.keys()):
			#for azi in sorted(muellerMatrices[AOI].keys()):
				#for wl in sorted(muellerMatrices[AOI][azi].keys()):
					#mm = np.reshape(np.asarray(muellerMatrices[AOI][azi][wl]), (4,4)) 
					#if checkDepolarizingMuellerMatrix(mm):
						#muellerMatrices[AOI][azi][wl] = np.copy(mm)
					#else:
						#muellerMatrices[AOI][azi][wl] = np.copy(cloudeDecomposition(mm))
	
	
	
		for AOI in sorted(self.data.keys()):
			for azi in sorted(self.data[AOI].keys()):
				for wl in sorted(self.data[AOI][azi].keys()):
					mm = np.reshape(np.asarray(self.data[AOI][azi][wl]), (4,4))
					self.data[AOI][azi][wl] = np.copy(mm)
#******************************************************************************************************************************************#

	def MFR(self,directory,wvlcutofmin,wvlcutofmax):
		'''
		Read Mueller matrix data from the collection of files stored in the same folder.
		Each file must represent a separate azimuthal angle and the name of the file ends with the azimuthal angle.
		Saves Mueller matrix data into self.data
		'''
		os.chdir(directory)
		files=os.listdir()
		if 'collection.p' in files:
			with open('collection.p','rb') as f:
				self.data=pickle.load(f)
				print('Loaded from pickle')
		else:
			for file in files:
				if file.endswith('.dat'):
			
					print ("processing Azi {}".format(re.findall("-(\d+).dat", file)[0])) 
					with open(directory+'/'+file,'r') as f:	   # Azimuthal angles are determined from file names
					
						for lines in f:
					
							if lines[0]=='m':									# skip first non mm-values
								
								lines = re.split(r'\t+', lines)					 
					
								AOI=lines[2]
								if AOI not in self.data.keys():
									self.data[AOI]=dict()
							
							
								
								AZI=float(re.findall("-(\d+).dat", file)[0])
								if AZI not in self.data[AOI].keys():
									self.data[AOI][AZI]=dict()
							
							
								wvl=float(lines[1])									 
								if wvl>=wvlcutofmin and wvl<=wvlcutofmax and wvl not in self.data[AOI][AZI].keys():
									self.data[AOI][AZI][wvl]=list()
							
							
					
					
						
								if lines[0] in self.muellerNames and wvl<=wvlcutofmax and wvl>=wvlcutofmin:
								#if wvl==wvlcutofmax:
									#muellerMatrices[AOI][AZI][wvl].append(0.5)
								#elif wvl==wvlcutofmin:
									# muellerMatrices[AOI][AZI][wvl].append(-0.5) 
								#else:	  
									self.data[AOI][AZI][wvl].append(float(lines[3])) 
						
						
			pickle.dump( self.data, open("collection.p", "wb" ))				
		



		
		#for AOI in sorted(muellerMatrices.keys()):
			#for azi in sorted(muellerMatrices[AOI].keys()):
				#for wl in sorted(muellerMatrices[AOI][azi].keys()):
					#mm = np.reshape(np.asarray(muellerMatrices[AOI][azi][wl]), (4,4)) 
					#if checkDepolarizingMuellerMatrix(mm):
						#muellerMatrices[AOI][azi][wl] = np.copy(mm)
					#else:
						#muellerMatrices[AOI][azi][wl] = np.copy(cloudeDecomposition(mm))
		
		
		
		for AOI in sorted(self.data.keys()):
			for azi in sorted(self.data[AOI].keys()):
				for wl in sorted(self.data[AOI][azi].keys()):
					mm = np.reshape(np.asarray(self.data[AOI][azi][wl]), (4,4))
					self.data[AOI][azi][wl] = np.copy(mm)
#******************************************************************************************************************************************#


	def Dispersion(self,Az,wvlcutofmin,wvlcutofmax,cMap,scalemax):
		'''
		Calculates P and S polarization intensities based on the Mueller matrix data.
		Returns a dispersion plot of P and S light.
		'''
		muellerMatrices=self.muellerNames
		plt.rcParams['axes.facecolor'] = 'black'
		plt.rcParams['font.size']=15
		
		SoutP=[]
		SoutS=[]
		muellerMatrices=self.data
		
		SinP=[1,1,0,0]
		SinS=[1,-1,0,0]
		
			
		#Unnormalize MM and reshape it
		for AOI in sorted(muellerMatrices.keys()):
			for azi in sorted(muellerMatrices[AOI].keys()):
				for wl in sorted(muellerMatrices[AOI][azi].keys()):
					mm = np.reshape(np.asarray(muellerMatrices[AOI][azi][wl]), (4,4))
					i=mm[0][0]
					mm=mm*i
					mm[0][0]=i	
					muellerMatrices[AOI][azi][wl] = np.copy(mm)
		
		
		
		
		#Calculate Sout			  
		for AOI in sorted(muellerMatrices.keys()):		
			for wvl in list(sorted(muellerMatrices[AOI][Az].keys())):
				SoutP.append(muellerMatrices[AOI][Az][wvl].dot(SinP))
				SoutS.append(muellerMatrices[AOI][Az][wvl].dot(SinS))
		
		#Create datasets				
		Pp=[x[0] for x in SoutP]
		Sp=[x[0] for x in SoutS]
		AOIs=list(map(float,list(sorted(muellerMatrices.keys()))))
		AOIs=np.asarray(AOIs)
		AOIsin=np.sin(np.radians(AOIs))
		k=[]
		wvls=np.asarray(list(sorted(muellerMatrices[AOI][Az].keys())))
		
		for el in AOIsin:
			for wvl in wvls:
				k.append((2*m.pi/wvl)*el)
		k=np.asarray(k)
		k=k.reshape(len(k)/len(wvls),len(wvls))
		
		
		dataP=np.array(Pp)
		dataP=dataP.reshape(len(AOIs),len(dataP)/len(AOIs))
		dataS=np.array(Sp)
		dataS=dataS.reshape(len(AOIs),len(dataS)/len(AOIs))
		
		
		wvls=np.asarray(list(sorted(muellerMatrices[AOI][Az].keys()))*len(AOIs))
		energy=np.asarray([3*10**8*6.62607004e-34/(wvl*10**(-9))*6.242e18 for wvl in wvls ])
		energy=energy.reshape([len(AOIs),len(energy)/len(AOIs)])
		wvls=np.reshape(wvls,[len(AOIs),len(wvls)/len(AOIs)])
		
		# Plot P polarization
		fig=plt.figure()
		plt.subplot(121)
		cp = plt.contourf(k, energy, dataP,500,vmin=0,vmax=scalemax,cmap=cMap)
		plt.xlabel(r'$K_{\parallel} (10^{9} m^{-1})$',size=20)
		plt.ylabel('Energy (eV)',size=20)
		plt.title('P-pol.')
		
		
		#Plot S polarization
		plt.subplot(122)
		cp = plt.contourf(k, energy, dataS,500,vmin=0,vmax=scalemax,cmap=cMap)
		plt.xlabel(r'$K_{\parallel} (10^{9} m^{-1})$',size=20)
		plt.ylabel('Energy (eV)',size=20)
		plt.title('S-pol.')
		
		
		
		#Define the colourmap and colorbar
		cax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
		cNorm = mpl.colors.Normalize(vmin=0, vmax=scalemax)
		cb1 = mpl.colorbar.ColorbarBase(cax, norm=cNorm,cmap=cMap)
		axes = plt.gca()
		
		plt.suptitle('Pp and Sp Intensity Dispersion Azi={}'.format(Az),fontsize=24)
		plt.show()	
		
		
		
#******************************************************************************************************************************************#
	def DispersionAOI(self,Az,wvlcutofmin,wvlcutofmax,cMap,scalemax):
		'''
		Calculates P and S polarization intensities based on the Mueller matrix data.
		Returns a dispersion plot of P and S light.
		'''
		muellerMatrices=self.muellerNames
		plt.rcParams['axes.facecolor'] = 'black'
		plt.rcParams['font.size']=15
		
		SoutP=[]
		SoutS=[]
		muellerMatrices=self.data
		
		SinP=[1,1,0,0]
		SinS=[1,-1,0,0]
		
			
		#Unnormalize MM and reshape it
		for AOI in sorted(muellerMatrices.keys()):
			for azi in sorted(muellerMatrices[AOI].keys()):
				for wl in sorted(muellerMatrices[AOI][azi].keys()):
					mm = np.reshape(np.asarray(muellerMatrices[AOI][azi][wl]), (4,4))
					i=mm[0][0]
					mm=mm*i
					mm[0][0]=i	
					muellerMatrices[AOI][azi][wl] = np.copy(mm)
		
		
		
		
		#Calculate Sout			  
		for AOI in sorted(muellerMatrices.keys()):		
			for wvl in list(sorted(muellerMatrices[AOI][Az].keys())):
				SoutP.append(muellerMatrices[AOI][Az][wvl].dot(SinP))
				SoutS.append(muellerMatrices[AOI][Az][wvl].dot(SinS))
		
		#Create datasets				
		Pp=[x[0] for x in SoutP]
		Sp=[x[0] for x in SoutS]
		AOIs=list(map(float,list(sorted(muellerMatrices.keys()))))
		AOIs=np.asarray(AOIs)
		AOIsin=np.sin(np.radians(AOIs))
		k=[]
		wvls=np.asarray(list(sorted(muellerMatrices[AOI][Az].keys())))
		X,Y= np.meshgrid(wvl,AOIs)
		
		
		dataP=np.array(Pp)
		dataP=dataP.reshape(len(AOIs),len(dataP)/len(AOIs))
		dataS=np.array(Sp)
		dataS=dataS.reshape(len(AOIs),len(dataS)/len(AOIs))
		
		
		
		
		
		# Plot P polarization
		fig=plt.figure()
		plt.subplot(121)
		cp = plt.contourf(X, Y, dataP,500,vmin=0,vmax=scalemax,cmap=cMap)
		plt.xlabel(r'$K_{\parallel} (10^{9} m^{-1})$',size=20)
		plt.ylabel('Energy (eV)',size=20)
		plt.title('P-pol.')
		
		
		#Plot S polarization
		plt.subplot(122)
		cp = plt.contourf(X, Y, dataS,500,vmin=0,vmax=scalemax,cmap=cMap)
		plt.xlabel(r'$K_{\parallel} (10^{9} m^{-1})$',size=20)
		plt.ylabel('Energy (eV)',size=20)
		plt.title('S-pol.')
		
		
		
		#Define the colourmap and colorbar
		cax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
		cNorm = mpl.colors.Normalize(vmin=0, vmax=scalemax)
		cb1 = mpl.colorbar.ColorbarBase(cax, norm=cNorm,cmap=cMap)
		axes = plt.gca()
		
		plt.suptitle('Pp and Sp Intensity Dispersion Azi={}'.format(Az),fontsize=24)
		plt.show()		
#******************************************************************************************************************************************#

	def TdepRead(self,folder,wvlcutofmin,wvlcutofmax):
		'''
		Reads temperature resolved data files which contained in the same folder 'folder'.
		The name of the file must represent the temperature at which measured.
		After saves data in pickle and in hdf5 format.
		Saves a dictionary type data into the self.data container.
		
		'''
		files=os.listdir(folder)
		if 'allT.p' in files:
			with open(folder+'\\allT.p','rb') as f:
				self.data=pickle.load(f)
			print('Loaded from pickle file')
		
		else:
			for file in files:
				if file.endswith('.dat'):
					
					Temp=float(file[0:-4])
					if Temp not in self.data.keys():
						self.data[Temp]={}
						print('T={}'.format(Temp))
					with open(folder+'\\'+file,'r') as f:
						for lines in f:				
								if lines[0]=='m':								   # skip first non mm-values				 
									lines = re.split(r'\t+', lines)					   #formating the lines
						
									AOI=lines[2]
									if AOI not in self.data[Temp].keys():
										self.data[Temp][AOI]=dict()
							
									try:		
										AZI=float(lines[7])
									except IndexError:
										AZI=0.0
									finally:
										if AZI not in self.data[Temp][AOI].keys():
											self.data[Temp][AOI][AZI]=dict() 
									
									
									wvl=float(lines[1])									 
									if wvl>=wvlcutofmin and wvl<=wvlcutofmax and wvl not in self.data[Temp][AOI][AZI].keys():
										self.data[Temp][AOI][AZI][wvl]=list()
							
							
							
									if lines[0] in self.muellerNames and wvl<=wvlcutofmax and wvl>=wvlcutofmin:
										if lines[3]=='Infinity':
											self.data[Temp][AOI][AZI][wvl].append(float(1))
										else :						  
											self.data[Temp][AOI][AZI][wvl].append(float(lines[3]))					 
					
					
			os.chdir(folder)
			pickle.dump( self.data, open("allT.p", "wb" ) )
		
		for Temp in sorted(self.data.keys()):		
			for AOI in sorted(self.data[Temp].keys()):
				for azi in sorted(self.data[Temp][AOI].keys()):
					for wl in sorted(self.data[Temp][AOI][azi].keys()):
						try:
							self.data[Temp][AOI][azi][wl] = np.reshape(np.asarray(self.data[Temp][AOI][azi][wl]), (4,4))	
						except( ValueError):
							print('oops')
		print('Done')
	
#******************************************************************************************************************************************#	
	def FindHermitian(self,AOI,AZI,Temp,WVL):
		s0=np.asarray([[1,0],[0,1]]).reshape(2,2)
		s1=np.asarray([[0,1],[1,0]]).reshape(2,2)
		s2=np.asarray([[0,-1j],[1j,0]]).reshape(2,2)
		s3=np.asarray([[1,0],[0,-1]]).reshape(2,2)
		S=[s0,s1,s2,s3]
		if Temp==None:
			M=self.data[AOI][AZI][WVL]
		else:
			M=self.data[Temp][AOI][AZI][WVL]
		M[0][0]=1	
		h=np.empty((4,4),dtype='complex128')
		for i in range(4):
			for j in range(4):
				h+=M[i][j]*np.kron(S[i],np.conj(S[j]))
		H=h/4
		return H
				
#******************************************************************************************************************************************#		
	def isPhysical(self,AOI,AZI,Temp,WVL):
		# we assume the matrix is valid and change the
		# validity if one or more of the 4 conditions
		# are not met.
		if Temp==None:
			mm=self.data[AOI][AZI][WVL]
		else:
			mm=self.data[Temp][AOI][AZI][WVL]
		mm[0][0] = 1.0
		
		matrixIsValid = True
		b = np.sqrt(mm[0][1]**2 + mm[0][2]**2 + mm[0][3]**2) 
		
		# 1st condition
		if not (np.trace(mm.dot(mm.T)) <= 4.0*mm[0][0]**2):
			matrixIsValid = False
			print('Condition 1')
			
		# 2nd condition
		for i in range(4):
			for j in range(4):
				if not (mm[0][0] >= np.abs(mm[i][j])):
					matrixIsValid = False
					print('Condition 2')
		
		# 3rd condition
		if not (mm[0][0]**2 >= b**2):
			matrixIsValid = False
			print('Condition 3')
			
		# 4th condition
		ss = np.zeros(3)	# second sum, see pdf
		fs = np.zeros(3)	# first sum, see pdf
		for j in range(1,4):
			for k in range(1,4):
				ss[k-1]=mm[j][k]*(mm[0][j]/b)
			fs[j-1] = mm[0][j] - np.sum(ss)
		sumRes = np.sum(fs)
		if not ((mm[0][0]-b)**2 >= sumRes):
			matrixIsValid = False
			print('Condition 4')
		
		return matrixIsValid
#******************************************************************************************************************************************#
	def PolarDecomposition(self,AOI,Azi,wvl):
		I=np.asarray(np.eye(3))
		self.M_d={}
		self.M_D={}
		self.M_r={}
		for T in sorted(self.data.keys()):
			print('decomposing data at T={}'.format(T))
			mm=self.data[T][AOI][Azi][wvl]
			mm[0][0]=1
		
			D=np.asarray(np.empty([1,3]))
			Dv=np.asarray([mm[0][1],mm[0][2],mm[0][3]])
			dvs=np.sum(np.square(Dv))
			D[:]=[x/(dvs) for x in Dv]
		
			ds=np.sum(np.square(D))
		
		
			try:
				m_d=(m.sqrt(1-ds))*I+((1-m.sqrt(1-ds))*(D@D.T)[0][0])
			except ValueError:
				m_d=np.eye(3)
			
			M_d=np.zeros([4,4])
			M_d[0][0]=1
			M_d[0][1:4]=Dv 
			M_d[1:4,0]=Dv
			M_d[1:4,1:4]=m_d
			self.M_d[T]=M_d
		
		
		
			M_p=mm@(np.linalg.inv(M_d))
		
			M_p[0][0]=1
			
			m_p=M_p[1:4,1:4]
		
			
			#Depolarization submatrix m_D
			M=m_p@m_p.T
			try:
				eigval=np.linalg.eig(M)[0]
			except np.linalg.linalg.LinAlgError:
				eigval=[0,0,0]
			
			
			m_D1=(M+(m.sqrt(eigval[0]*eigval[1])+m.sqrt(eigval[1]*eigval[2])+m.sqrt(eigval[2]*eigval[0]))*I)
			
			if np.linalg.det(M)<0:
				m_D1=-m_D1
			
			m_D2=(m.sqrt(eigval[0]) +m.sqrt(eigval[1])+m.sqrt(eigval[2]))*M+m.sqrt(eigval[0]*eigval[1]*eigval[2])*I
			m_D=np.linalg.inv(m_D1)@m_D2
			M_D=np.zeros([4,4])
			M_D[0][0]=1
			M_D[1:4,1:4]=m_D
			self.M_D[T]=M_D
			
			
			
			#Finally calculate the retardence matrix
			
			M_r=np.linalg.inv(M_D)@M_p
			self.M_r[T]=M_r
					
		
		return
#******************************************************************************************************************************************#
	def PlotMMSlice(self,what,AOI,Azi):
		import matplotlib
		matplotlib.rcParams['font.size']=15
		figure= plt.figure()
		if what=='Lm':
			for aoi in AOI:
				for azi in Azi:
					k=1
					for i in range(0,4):
						for j in range(0,4):
							if i==0 and j==0:
								plt.subplot(4,4,k)
								plt.plot(sorted(list(self.Lm[aoi][azi].keys())),[self.Lm[aoi][azi][wl][i][j] for wl in self.Lm[aoi][azi].keys()],label='Azi={}'.format(azi))
								plt.legend()
								k+=1
							else:
								plt.subplot(4,4,k)
								plt.plot(sorted(list(self.Lm[aoi][azi].keys())),[self.Lm[aoi][azi][wl][i][j] for wl in self.Lm[aoi][azi].keys()])
								#plt.legend()
								k+=1
		
		else:
			for aoi in AOI:
				for azi in Azi:
					k=1
					for i in range(0,4):
						for j in range(0,4):
							if i==0 and j==0:
								plt.subplot(4,4,k)
								plt.plot(sorted(list(self.data[aoi][azi].keys())),[self.data[aoi][azi][wl][i][j] for wl in self.data[aoi][azi].keys()],label='Azi={}'.format(azi))
								plt.legend()
								k+=1
							else:
								plt.subplot(4,4,k)
								plt.plot(sorted(list(self.data[aoi][azi].keys())),[self.data[aoi][azi][wl][i][j] for wl in self.data[aoi][azi].keys()])
								#plt.legend()
								k+=1
		figure.text(0.45,0.00,'Wavelength (nm)')
		return
					
			
					
			