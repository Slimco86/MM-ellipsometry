import numpy as np
from scipy.linalg import logm, expm
import math
import pandas as pd
import re
import gc
import matplotlib as mpl
import os
import matplotlib.pyplot as plt
import math as m
import cmath
from tkinter import *

global angles
angles=[0,360,5]
gratingPeriod=312.0
order = 1
n1=1

def getMuellerMatrixElementNames():
	return ['mm11', 'mm12', 'mm13', 'mm14',
            'mm21', 'mm22', 'mm23', 'mm24',
            'mm31', 'mm32', 'mm33', 'mm34',
            'mm41', 'mm42', 'mm43', 'mm44']


     
def checkDepolarizingMuellerMatrix(mm):
	# we assume the matrix is valid and change the
	# validity if one or more of the 4 conditions
	# are not met.
	
	mm[0][0] = 1.0
	
	matrixIsValid = True
	b = np.sqrt(mm[0][1]**2 + mm[0][2]**2 + mm[0][3]**2) 
	
	# 1st condition
	if not (np.trace(mm.dot(mm.T)) <= 4.0*mm[0][0]**2):
		matrixIsValid = False
		
	# 2nd condition
	for i in range(4):
		for j in range(4):
			if not (mm[0][0] >= np.abs(mm[i][j])):
				matrixIsValid = False
	
	# 3rd condition
	if not (mm[0][0]**2 >= b**2):
		matrixIsValid = False
		
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
	
	return matrixIsValid

def cloudeDecomposition(mm):
	M = np.copy(mm)
	# calculate coherency matrix H of Mueller Matrix M
	M[0][0] = 1.0
	
	H = np.reshape(np.zeros(16, dtype='complex64'), (4, 4))
	H[0][0] = M[0][0]    + M[1][1]    + M[2][2]    + M[3][3]
	H[0][1] = M[0][1]    + M[1][0]    - 1j*M[2][3] + 1j*M[3][2]
	H[0][2] = M[0][2]    + M[2][0]    + 1j*M[1][3] - 1j*M[3][1]  
	H[0][3] = M[0][3]    - 1j*M[1][2] + 1j*M[2][1] + M[3][0]  
	
	H[1][0] = M[0][1]    + M[1][0]    + 1j*M[2][3] - 1j*M[3][2]
	H[1][1] = M[0][0]    + M[1][1]    - M[2][2]    - M[3][3]
	H[1][2] = 1j*M[0][3] + M[1][2]    + M[2][1]    - 1j*M[3][0]
	H[1][3] = 1j*M[2][0] - 1j*M[0][2] + M[1][3]    + M[3][1]
	
	H[2][0] = M[0][2]    + M[2][0]    - 1j*M[1][3] + 1j*M[3][1]
	H[2][1] = M[1][2]    - 1j*M[0][3] + M[2][1]    + 1j*M[3][0]
	H[2][2] = M[0][0]    - M[1][1]    + M[2][2]    -  M[3][3]
	H[2][3] = 1j*M[0][1] - 1j*M[1][0] + M[2][3]    + M[3][2]
	
	H[3][0] = M[0][3]    + 1j*M[1][2] - 1j*M[2][1] + M[3][0]
	H[3][1] = 1j*M[0][2] - 1j*M[2][0] + M[1][3]    + M[3][1]
	H[3][2] = 1j*M[1][0] - 1j*M[0][1] + M[2][3]    + M[3][2]
	H[3][3] = M[0][0]    - M[1][1]    - M[2][2]    + M[3][3]
	
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
	
	JtoMtransformer = np.array([[1.0, 0.0, 0.0,  1.0],
							    [1.0, 0.0, 0.0, -1.0],
							    [0.0, 1.0, 1.0,  0.0],
							    [0.0, 1j , -1j,  0.0]])
	invJtoMtransformer = np.array([[ 0.5,  0.5,  0.0,  0.0   ],
								   [ 0.0,  0.0,  0.5, -0.5*1j],
								   [ 0.0,  0.0,  0.5,  0.5*1j],
								   [ 0.5, -0.5,  0.0,  0.0   ]])
	
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
	
	return MM1,MM2,MM3,MM4,eigValH





def differentialDecomposition(mm):
    print("Differential decomposition")
    g=[1,0,0,0,
       0,-1,0,0,
       0,0,-1,0,
       0,0,0,-1]
    G=np.reshape(g,(4,4))
    L=dict()
    L_u=dict()
    L_m=dict()
    
    
    for AOI in mm.keys():
        if AOI not in L.keys():
            L[AOI]=dict()
            L_m[AOI]=dict()
            L_u[AOI]=dict()
        for AZI in sorted(mm[AOI].keys()):
            if AZI not in L[AOI].keys():
                L[AOI][AZI]=dict()
                L_m[AOI][AZI]=dict()
                L_u[AOI][AZI]=dict()            
            for wvl in sorted(mm[AOI][AZI].keys()):
                if wvl not in L[AOI][AZI].keys():
                    L[AOI][AZI][wvl]=np.empty([4,4],dtype=complex)
                    L_m[AOI][AZI][wvl]=np.empty([4,4],dtype=complex)
                    L_u[AOI][AZI][wvl]=np.empty([4,4],dtype=complex)
                        
                L[AOI][AZI][wvl] = logm(mm[AOI][AZI][wvl])    # Right option 
                L_m[AOI][AZI][wvl]=1/2*(L[AOI][AZI][wvl]-(G@np.transpose(L[AOI][AZI][wvl]))@G)
                L_u[AOI][AZI][wvl]=1/2*(L[AOI][AZI][wvl]+(G@np.transpose(L[AOI][AZI][wvl]))@G)
                
    
    return(L,L_m,L_u)  
    

def PolarDecomposition(mm):
    I=np.asarray(np.eye(3))
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
    
    
    
    #Finally calculate the retardence matrix
    
    M_r=np.linalg.inv(M_D)@M_p
    
    #Filter very small values
    
    for index,element in np.ndenumerate(M_r):
        if abs(element)<0.001:
            element=0
    for index,element2 in np.ndenumerate(M_d):
        if abs(element2)<0.001:
            element2=0
            
    for index,element3 in np.ndenumerate(M_D):
        if abs(element3)<0.001:
            element3=0        
    
    return M_d,M_r,M_D



def createMuellerMatrixPlotM_D(mm):
    # generate dictionary for each AOI to plot
    muellerGraphs = dict()
    for AOI in mm.keys():
        muellerGraphs[AOI] = dict()

	# generate dictionary entry for each mm-element
    for AOI in muellerGraphs.keys():
        for element in getMuellerMatrixElementNames():
            muellerGraphs[AOI][element] = list()

    # now insert all mm-values of each element
    for AOI in sorted(mm.keys()):
        for AZI in sorted(mm[AOI].keys()):
            for wvl in sorted(mm[AOI][AZI].keys()):
                currMM = mm[AOI][AZI][wvl].reshape(-1)
                for element, i in zip(sorted(getMuellerMatrixElementNames()), range(16)):
                    muellerGraphs[AOI][element].append(currMM[i])

    AOIkeys = list(mm.keys())
    #print(AOIkeys)
    aziKeys = list(mm[AOIkeys[0]].keys())
    
    #print(aziKeys)
    
    azimuth = np.radians(sorted(np.asarray(list(mm[AOIkeys[0]].keys()))))
    zenith  = np.asarray(sorted(np.asarray(list(mm[AOIkeys[0]][aziKeys[0]].keys()))))
    
    #print(zenith)
    # adjust plot settings
    mpl.rcParams['axes.labelsize']  = 6
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 6
    
    
        
    print('Depolarization Diatennuation plot')
    for AOI in sorted(muellerGraphs.keys()):
        fig, axs = plt.subplots(4, 4, subplot_kw=dict(projection='polar'))	
        for element, ax in zip(sorted(muellerGraphs[AOI].keys()), iter(axs.reshape(-1))):
            data = np.reshape(muellerGraphs[AOI][element], (azimuth.size, zenith.size)).T
            cs = ax.contourf(azimuth, zenith, data, 100, cmap=plt.cm.jet)
            ax.grid(False)
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            thetaticks = np.arange(0,360,90)
            ax.set_thetagrids(thetaticks, frac=1.3)
            minmax = np.asarray([np.min(data), np.max(data)])
            fig.colorbar(cs, ax=ax, orientation='vertical', ticks=[minmax[0], np.mean(minmax), minmax[1]], pad=0.2)
            
            #plt.tight_layout()
            #print('Saving L_m %s plot at %s'%(element, AOI))
            #plt.savefig(saveDirectory + 'L_m_AOI_%03.1f.png'%float(AOI), dpi=600, bbox_inches='tight')
        plt.suptitle("Depolarization at AOI-%s from %s to %s wvl."%(AOI,min(zenith),max(zenith)))
        plt.show()
        plt.close()
    print('Finished!!!!!!!!!')
    
    
def createMuellerMatrixPlotM_r(mm):
    # generate dictionary for each AOI to plot
    muellerGraphs = dict()
    for AOI in mm.keys():
        muellerGraphs[AOI] = dict()

	# generate dictionary entry for each mm-element
    for AOI in muellerGraphs.keys():
        for element in getMuellerMatrixElementNames():
            muellerGraphs[AOI][element] = list()

    # now insert all mm-values of each element
    for AOI in sorted(mm.keys()):
        for AZI in sorted(mm[AOI].keys()):
            for wvl in sorted(mm[AOI][AZI].keys()):
                currMM = mm[AOI][AZI][wvl].reshape(-1)
                for element, i in zip(sorted(getMuellerMatrixElementNames()), range(16)):
                    muellerGraphs[AOI][element].append(currMM[i])

    AOIkeys = list(mm.keys())
    #print(AOIkeys)
    aziKeys = list(mm[AOIkeys[0]].keys())
    
    #print(aziKeys)
    
    azimuth = np.radians(sorted(np.asarray(list(mm[AOIkeys[0]].keys()))))
    zenith  = np.asarray(sorted(np.asarray(list(mm[AOIkeys[0]][aziKeys[0]].keys()))))
    
    #print(zenith)
    # adjust plot settings
    mpl.rcParams['axes.labelsize']  = 6
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 6
    
    
        
    print('Building Retardance plot')
    for AOI in sorted(muellerGraphs.keys()):
        fig, axs = plt.subplots(4, 4, subplot_kw=dict(projection='polar'))	
        for element, ax in zip(sorted(muellerGraphs[AOI].keys()), iter(axs.reshape(-1))):
            data = np.reshape(muellerGraphs[AOI][element], (azimuth.size, zenith.size)).T
            cs = ax.contourf(azimuth, zenith, data, 100, cmap=plt.cm.jet)
            ax.grid(False)
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            thetaticks = np.arange(0,360,90)
            ax.set_thetagrids(thetaticks, frac=1.3)
            minmax = np.asarray([np.min(data), np.max(data)])
            fig.colorbar(cs, ax=ax, orientation='vertical', ticks=[minmax[0], np.mean(minmax), minmax[1]], pad=0.2)
            
            #plt.tight_layout()
            #print('Saving L_m %s plot at %s'%(element, AOI))
            #plt.savefig(saveDirectory + 'L_m_AOI_%03.1f.png'%float(AOI), dpi=600, bbox_inches='tight')
        plt.suptitle("Retardance at AOI-%s from %s to %s wvl."%(AOI,min(zenith),max(zenith)))
        plt.show()
        plt.close()
    print('Finished!!!!!!!!!')
    

def createMuellerMatrixPlotM_d(mm):
    # generate dictionary for each AOI to plot
    muellerGraphs = dict()
    for AOI in mm.keys():
        muellerGraphs[AOI] = dict()

	# generate dictionary entry for each mm-element
    for AOI in muellerGraphs.keys():
        for element in getMuellerMatrixElementNames():
            muellerGraphs[AOI][element] = list()

    # now insert all mm-values of each element
    for AOI in sorted(mm.keys()):
        for AZI in sorted(mm[AOI].keys()):
            for wvl in sorted(mm[AOI][AZI].keys()):
                currMM = mm[AOI][AZI][wvl].reshape(-1)
                for element, i in zip(sorted(getMuellerMatrixElementNames()), range(16)):
                    muellerGraphs[AOI][element].append(currMM[i])

    AOIkeys = list(mm.keys())
    #print(AOIkeys)
    aziKeys = list(mm[AOIkeys[0]].keys())
    
    #print(aziKeys)
    
    azimuth = np.radians(sorted(np.asarray(list(mm[AOIkeys[0]].keys()))))
    zenith  = np.asarray(sorted(np.asarray(list(mm[AOIkeys[0]][aziKeys[0]].keys()))))
    
    #print(zenith)
    # adjust plot settings
    mpl.rcParams['axes.labelsize']  = 6
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 6
    
    
        
    print('Building Diatennuation plot')
    for AOI in sorted(muellerGraphs.keys()):
        fig, axs = plt.subplots(4, 4, subplot_kw=dict(projection='polar'))	
        for element, ax in zip(sorted(muellerGraphs[AOI].keys()), iter(axs.reshape(-1))):
            data = np.reshape(muellerGraphs[AOI][element], (azimuth.size, zenith.size)).T
            cs = ax.contourf(azimuth, zenith, data, 100, cmap=plt.cm.jet)
            ax.grid(False)
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            thetaticks = np.arange(0,360,90)
            ax.set_thetagrids(thetaticks, frac=1.3)
            minmax = np.asarray([np.min(data), np.max(data)])
            fig.colorbar(cs, ax=ax, orientation='vertical', ticks=[minmax[0], np.mean(minmax), minmax[1]], pad=0.2)
            
            #plt.tight_layout()
            #print('Saving L_m %s plot at %s'%(element, AOI))
            #plt.savefig(saveDirectory + 'L_m_AOI_%03.1f.png'%float(AOI), dpi=600, bbox_inches='tight')
        plt.suptitle("Diatennuation at AOI-%s from %s to %s wvl."%(AOI,min(zenith),max(zenith)))
        plt.show()
        plt.close()
    print('Finished!!!!!!!!!')
    
    
    
def createMuellerMatrixPlotL_u(mm):
    # generate dictionary for each AOI to plot
    muellerGraphs = dict()
    for AOI in mm.keys():
        muellerGraphs[AOI] = dict()

	# generate dictionary entry for each mm-element
    for AOI in muellerGraphs.keys():
        for element in getMuellerMatrixElementNames():
            muellerGraphs[AOI][element] = list()

    # now insert all mm-values of each element
    for AOI in sorted(mm.keys()):
        for AZI in sorted(mm[AOI].keys()):
            for wvl in sorted(mm[AOI][AZI].keys()):
                currMM = mm[AOI][AZI][wvl].reshape(-1)
                for element, i in zip(sorted(getMuellerMatrixElementNames()), range(16)):
                    muellerGraphs[AOI][element].append(currMM[i])

    AOIkeys = list(mm.keys())
    #print(AOIkeys)
    aziKeys = list(mm[AOIkeys[0]].keys())
    
    #print(aziKeys)
    
    azimuth = np.radians(sorted(np.asarray(list(mm[AOIkeys[0]].keys()))))
    zenith  = np.asarray(sorted(np.asarray(list(mm[AOIkeys[0]][aziKeys[0]].keys()))))
    
    #print(zenith)
    # adjust plot settings
    mpl.rcParams['axes.labelsize']  = 6
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 6
    
    saveDirectory = 'D:/Ievgen Ell/Measurements/NanoDiscs/MM_transmition/OD5/result_MM_plot/'
    if not os.path.isdir(saveDirectory):
        os.makedirs(saveDirectory)
    print('Building L_u plot')
    for AOI in sorted(muellerGraphs.keys()):
        fig, axs = plt.subplots(4, 4, subplot_kw=dict(projection='polar'))	
        for element, ax in zip(sorted(muellerGraphs[AOI].keys()), iter(axs.reshape(-1))):
            data = np.reshape(muellerGraphs[AOI][element], (azimuth.size, zenith.size)).T
            cs = ax.contourf(azimuth, zenith, data, 100, cmap=plt.cm.jet)
            ax.grid(False)
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            thetaticks = np.arange(0,360,90)
            ax.set_thetagrids(thetaticks, frac=1.3)
            minmax = np.asarray([np.min(data), np.max(data)])
            fig.colorbar(cs, ax=ax, orientation='vertical', ticks=[minmax[0], np.mean(minmax), minmax[1]], pad=0.2)
            
            #plt.tight_layout()
            #print('Saving L_u %s plot at %s'%(element, AOI))
            #plt.savefig(saveDirectory + 'L_u_AOI_%03.1f.png'%float(AOI), dpi=600, bbox_inches='tight')
        plt.suptitle("L_u plot at AOI-%s from %s to %s wvl."%(AOI,min(zenith),max(zenith)))
        plt.show()
        plt.close()
    print('Finished!!!!!!!!!')


def createMuellerMatrixPlotL_m(mm):
    # generate dictionary for each AOI to plot
    muellerGraphs = dict()
    for AOI in mm.keys():
        muellerGraphs[AOI] = dict()

	# generate dictionary entry for each mm-element
    for AOI in muellerGraphs.keys():
        for element in getMuellerMatrixElementNames():
            muellerGraphs[AOI][element] = list()

    # now insert all mm-values of each element
    for AOI in sorted(mm.keys()):
        for AZI in sorted(mm[AOI].keys()):
            for wvl in sorted(mm[AOI][AZI].keys()):
                currMM = mm[AOI][AZI][wvl].reshape(-1)
                for element, i in zip(sorted(getMuellerMatrixElementNames()), range(16)):
                    muellerGraphs[AOI][element].append(currMM[i])

    AOIkeys = list(mm.keys())
    #print(AOIkeys)
    aziKeys = list(mm[AOIkeys[0]].keys())
    
    #print(aziKeys)
    
    azimuth = np.radians(sorted(np.asarray(list(mm[AOIkeys[0]].keys()))))
    zenith  = np.asarray(sorted(np.asarray(list(mm[AOIkeys[0]][aziKeys[0]].keys()))))
    
    #print(zenith)
    # adjust plot settings
    mpl.rcParams['axes.labelsize']  = 6
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 6
    
    saveDirectory = 'D:/Ievgen Ell/Measurements/NanoDiscs/MM_transmition/OD5/result_MM_plot/'
    if not os.path.isdir(saveDirectory):
        os.makedirs(saveDirectory)
        
    print('Building L_m plot')
    for AOI in sorted(muellerGraphs.keys()):
        fig, axs = plt.subplots(4, 4, subplot_kw=dict(projection='polar'))	
        for element, ax in zip(sorted(muellerGraphs[AOI].keys()), iter(axs.reshape(-1))):
            data = np.reshape(muellerGraphs[AOI][element], (azimuth.size, zenith.size)).T
            cs = ax.contourf(azimuth, zenith, data, 100, cmap=plt.cm.jet)
            ax.grid(False)
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            thetaticks = np.arange(0,360,90)
            ax.set_thetagrids(thetaticks, frac=1.3)
            minmax = np.asarray([np.min(data), np.max(data)])
            fig.colorbar(cs, ax=ax, orientation='vertical', ticks=[minmax[0], np.mean(minmax), minmax[1]], pad=0.2)
            
            #plt.tight_layout()
            #print('Saving L_m %s plot at %s'%(element, AOI))
            #plt.savefig(saveDirectory + 'L_m_AOI_%03.1f.png'%float(AOI), dpi=600, bbox_inches='tight')
        plt.suptitle("L_m plot at AOI-%s from %s to %s wvl."%(AOI,min(zenith),max(zenith)))
        plt.show()
        plt.close()
    print('Finished!!!!!!!!!')    
    
    
    
    
    
def createMuellerMatrixPlot(mm):
    # generate dictionary for each AOI to plot
    muellerGraphs = dict()
    for AOI in mm.keys():
        muellerGraphs[AOI] = dict()

	# generate dictionary entry for each mm-element
    for AOI in muellerGraphs.keys():
        for element in getMuellerMatrixElementNames():
            muellerGraphs[AOI][element] = list()

    # now insert all mm-values of each element
    for AOI in sorted(mm.keys()):
        for AZI in sorted(mm[AOI].keys()):
            for wvl in sorted(mm[AOI][AZI].keys()):
                currMM = mm[AOI][AZI][wvl].reshape(-1)
                for element, i in zip(sorted(getMuellerMatrixElementNames()), range(16)):
                    muellerGraphs[AOI][element].append(currMM[i])

    AOIkeys = list(mm.keys())
    #print(AOIkeys)
    aziKeys = list(mm[AOIkeys[0]].keys())
    
    #print(aziKeys)
    
    azimuth = np.radians(sorted(np.asarray(list(mm[AOIkeys[0]].keys()))))
    zenith  = np.asarray(sorted(np.asarray(list(mm[AOIkeys[0]][aziKeys[0]].keys()))))
    
    #print(zenith)
    # adjust plot settings
    mpl.rcParams['axes.labelsize']  = 6
    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6
    mpl.rcParams['legend.fontsize'] = 6
    
    saveDirectory = 'D:/Ievgen Ell/Measurements/NanoDiscs/MM_transmition/OD5/result_MM_plot/'
    if not os.path.isdir(saveDirectory):
        os.makedirs(saveDirectory)
    print('Building Measured MM plot')
    for AOI in sorted(muellerGraphs.keys()):
        
        fig, axs = plt.subplots(4, 4, subplot_kw=dict(projection='polar'))
        #azi=[]
        #azim2=[]
        #azim3=[]
        #azim4=[]
        #azim5=[]
        #azim6=[]
        #azim7=[]
        #azim8=[]
        #for i in range(0,46,1):
            #azi.append(m.radians(i))
            #azim2.append(m.radians(i)+m.pi)
            #azim3.append(m.radians(i)+(m.pi)/2)
            #azim4.append(m.radians(i)+(m.pi)*3/2)
            #azim5.append(m.radians(i)-(m.pi)/4)
            #azim6.append(m.radians(i)-(m.pi)*3/4)
            #azim7.append(m.radians(i)-(m.pi)*5/4)
            #azim8.append(m.radians(i)+(m.pi)/4)
            
            #bz=ReighleyOnBZ(AOI,439,0.0,2,1.4)
            #bz2=ReighleyOnBZ(AOI,439,0.0,4,1.0)
            #bz3=ReighleyOnBZ(AOI,621,45.0,1,1.4)
            #bz4=ReighleyOnBZ(AOI,439,0.0,2,1.4)
            #bz5=ReighleyOnBZ(AOI,621,45.0,2,1.4)
        
        for element, ax in zip(sorted(muellerGraphs[AOI].keys()), iter(axs.reshape(-1))):
            data = np.reshape(muellerGraphs[AOI][element], (azimuth.size, zenith.size)).T
            cs = ax.contourf(azimuth, zenith, data, 200, cmap=plt.cm.jet)
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            thetaticks = np.arange(0,360,90)
            ax.set_thetagrids(thetaticks, frac=1.3)
            ax.grid(False)
            minmax = np.asarray([np.min(data), np.max(data)])
            fig.colorbar(cs, ax=ax, orientation='vertical', ticks=[minmax[0], minmax[1]], pad=0.2)
            
          
            #ax.plot(azi,bz,linewidth=2,color='k')
            #ax.plot(azim2,bz,linewidth=2,color='k')
            #ax.plot(azim3,bz,linewidth=2,color='k')
            #ax.plot(azim4,bz,linewidth=2,color='k')
            #ax.plot(azim5,sorted(bz, reverse=True),linewidth=2,color='k')
            #ax.plot(azim6,sorted(bz, reverse=True),linewidth=2,color='k')
            #ax.plot(azim7,sorted(bz, reverse=True),linewidth=2,color='k')
            #ax.plot(azim8,sorted(bz, reverse=True),linewidth=2,color='k')
            
            #ax.plot(azi,bz2,linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x+m.pi/2,azi)),bz2,linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x+3*m.pi/2,azi)),bz2,linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x+m.pi,azi)),bz2,linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x-m.pi/4,azi)),sorted(bz2, reverse=True),linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x-3*m.pi/4,azi)),sorted(bz2, reverse=True),linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x-5*m.pi/4,azi)),sorted(bz2, reverse=True),linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x+m.pi/4,azi)),sorted(bz2, reverse=True),linewidth=2,color='w')
            
            #ax.plot(azi,bz3,linewidth=2,color='k')
            #ax.plot(list(map(lambda x:x+m.pi/2,azi)),bz3,linewidth=2,color='k')
            #ax.plot(list(map(lambda x:x+3*m.pi/2,azi)),bz3,linewidth=2,color='k')
            #ax.plot(list(map(lambda x:x+m.pi,azi)),bz3,linewidth=2,color='k')
            #ax.plot(list(map(lambda x:x-m.pi/4,azi)),sorted(bz3, reverse=False),linewidth=2,color='k')
            #ax.plot(list(map(lambda x:x-3*m.pi/4,azi)),sorted(bz3, reverse=False),linewidth=2,color='k')
            #ax.plot(list(map(lambda x:x-5*m.pi/4,azi)),sorted(bz3, reverse=False),linewidth=2,color='k')
            #ax.plot(list(map(lambda x:x+m.pi/4,azi)),sorted(bz3, reverse=False),linewidth=2,color='k')
            
            #ax.plot(azi,bz4,linewidth=2,color='k')
            #ax.plot(azim2,bz4,linewidth=2,color='k')
            #ax.plot(azim3,bz4,linewidth=2,color='k')
            #ax.plot(azim4,bz4,linewidth=2,color='k')
            #ax.plot(azim5,sorted(bz4, reverse=True),linewidth=2,color='k')
            #ax.plot(azim6,sorted(bz4, reverse=True),linewidth=2,color='k')
            #ax.plot(azim7,sorted(bz4, reverse=True),linewidth=2,color='k')
            #ax.plot(azim8,sorted(bz4, reverse=True),linewidth=2,color='k')
           
            #ax.plot(azi,bz5,linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x+m.pi/2,azi)),bz5,linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x+3*m.pi/2,azi)),bz5,linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x+m.pi,azi)),bz5,linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x-m.pi/4,azi)),sorted(bz5, reverse=False),linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x-3*m.pi/4,azi)),sorted(bz5, reverse=False),linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x-5*m.pi/4,azi)),sorted(bz5, reverse=False),linewidth=2,color='w')
            #ax.plot(list(map(lambda x:x+m.pi/4,azi)),sorted(bz5, reverse=False),linewidth=2,color='w')
            print('Ploting Measured %s'%(element))
            #plt.tight_layout()
            #print('Saving %s plot at %s'%(element, AOI))
            #plt.savefig(saveDirectory + 'AOI_%03.1f.png'%float(AOI), dpi=600, bbox_inches='tight')
        plt.suptitle("Mueller matrix plot at AOI-%s from %s to %s wvl."%(AOI,min(zenith),max(zenith)))
        plt.show()
        plt.close()
    print('Finished!!!!!!!!!')
    
    
    
    
def oneFile(directory,filename,wvlcutofmin,wvlcutofmax,plotMM,diffDecomp,LuChipman,Clode):
    
    muellerNames = getMuellerMatrixElementNames() 
    muellerMatrices = dict()
    if Clode:
        CDM1={}
        CDM2={}
        CDM3={}
        CDM4={}
        eigValH={}
    if LuChipman:
        M_d={}
        M_r={}
        M_D={}
    
            
    print ("Processing File")
    with open(directory+'/'+filename+ '.dat') as f:
        for lines in f:             
            if lines[0]=='m':                                  # skip first non mm-values                
                lines = re.split(r'\t+', lines)                    #formating the lines
                 
                AOI=lines[2]
                if AOI not in muellerMatrices.keys():
                    muellerMatrices[AOI]=dict()
                    if Clode:
                        CDM1[AOI]={}
                        CDM2[AOI]={}
                        CDM3[AOI]={}
                        CDM4[AOI]={}
                        eigValH[AOI]={}
                    if LuChipman:
                        M_d[AOI]={}
                        M_r[AOI]={}
                        M_D[AOI]={}
                        
                        
                AZI=float(lines[7])
                if AZI not in muellerMatrices[AOI].keys():
                    muellerMatrices[AOI][AZI]=dict() 
                    if Clode:
                        CDM1[AOI][AZI]={}
                        CDM2[AOI][AZI]={}
                        CDM3[AOI][AZI]={}
                        CDM4[AOI][AZI]={}
                        eigValH[AOI][AZI]={}
                    
                    if LuChipman:
                        M_d[AOI][AZI]={}
                        M_r[AOI][AZI]={}
                        M_D[AOI][AZI]={}
                        
                wvl=float(lines[1])                                  
                if wvl>=wvlcutofmin and wvl<=wvlcutofmax and wvl not in muellerMatrices[AOI][AZI].keys():
                    muellerMatrices[AOI][AZI][wvl]=list()
                    if Clode:
                        CDM1[AOI][AZI][wvl]=[]
                        CDM2[AOI][AZI][wvl]=[]
                        CDM3[AOI][AZI][wvl]=[]
                        CDM4[AOI][AZI][wvl]=[]
                        eigValH[AOI][AZI][wvl]=[]
                        
                    if LuChipman:
                        M_d[AOI][AZI][wvl]=[]
                        M_r[AOI][AZI][wvl]=[]
                        M_D[AOI][AZI][wvl]=[]
               
             
                    
                if lines[0] in muellerNames and wvl<=wvlcutofmax and wvl>=wvlcutofmin:
                    muellerMatrices[AOI][AZI][wvl].append(float(lines[3]))
                
                    				
                
                
                    
                   
    #for AOI in sorted(muellerMatrices.keys()):
        #for azi in sorted(muellerMatrices[AOI].keys()):
            #for wl in sorted(muellerMatrices[AOI][azi].keys()):
                #mm = np.reshape(np.asarray(muellerMatrices[AOI][azi][wl]), (4,4)) 
                #if checkDepolarizingMuellerMatrix(mm):
                    #muellerMatrices[AOI][azi][wl] = np.copy(mm)
                #else:
                    #muellerMatrices[AOI][azi][wl] = np.copy(cloudeDecomposition(mm))
   


    for AOI in sorted(muellerMatrices.keys()):
        for azi in sorted(muellerMatrices[AOI].keys()):
            for wl in sorted(muellerMatrices[AOI][azi].keys()):
                mm = np.reshape(np.asarray(muellerMatrices[AOI][azi][wl]), (4,4))
				
                if Clode:
                    MM1,MM2,MM3,MM4,eigVal=cloudeDecomposition(mm)
                    eigValH[AOI][azi][wl]=list(eigVal)
                    CDM1[AOI][azi][wl]=(np.reshape(np.asarray(MM1), (4,4)))
                    CDM2[AOI][azi][wl]=(np.reshape(np.asarray(MM2), (4,4)))
                    CDM3[AOI][azi][wl]=(np.reshape(np.asarray(MM3), (4,4)))
                    CDM4[AOI][azi][wl]=(np.reshape(np.asarray(MM4), (4,4)))
                
                if LuChipman:
                    M_d[AOI][azi][wl],M_r[AOI][azi][wl],M_D[AOI][azi][wl]=PolarDecomposition(mm)
                  
                muellerMatrices[AOI][azi][wl] = np.copy(mm)
                    
                
    if LuChipman:
        createMuellerMatrixPlotM_d(M_d)
        createMuellerMatrixPlotM_r(M_r)
        createMuellerMatrixPlotM_D(M_D)
        
    if Clode:
        plt.figure(1)
        lambda1=[]
        lambda2=[]
        lambda3=[]
        lambda4=[]
        AOI=list(eigValH.keys())[0]
        azi=sorted(list(eigValH[AOI].keys()))[9]
        wvls=list(eigValH[AOI][azi].keys())
        for wvl in sorted(list( eigValH[AOI][azi].keys())):
            lambda1.append(eigValH[AOI][azi][wvl][0])
            lambda2.append(eigValH[AOI][azi][wvl][1])
            lambda3.append(eigValH[AOI][azi][wvl][2])
            lambda4.append(eigValH[AOI][azi][wvl][3])
            
        plt.plot(wvls,lambda1,label="lambda1", color='r')
        plt.plot(wvls,lambda2,label="lambda2", color='g')
        plt.plot(wvls,lambda3,label="lambda3", color='b')
        plt.plot(wvls,lambda4,label="lambda4", color='k')
        plt.legend()
        plt.title('Cloude eigenvalues at AOI-%s, Azi-%s'%(AOI,azi))
        plt.show()
        #createMuellerMatrixPlot(CDM1)
        #createMuellerMatrixPlot(CDM2)
        #createMuellerMatrixPlot(CDM3)
        #createMuellerMatrixPlot(CDM4)
        #print( eigValH[0],eigValH[1],eigValH[2],eigValH[3])
    
    if plotMM:
        createMuellerMatrixPlot(muellerMatrices)
    
    if diffDecomp:
        L,L_m,L_u=differentialDecomposition(muellerMatrices)
        createMuellerMatrixPlotL_m(L_m)
        createMuellerMatrixPlotL_u(L_u)
        
   
    del muellerMatrices
    if Clode:
        del CDM1
        del CDM2
        del CDM3
        del CDM4
    if LuChipman:
        del M_d,M_r,M_D
    
    return     
 



   
#give the azimuthal angels in format e.g.: [0,180,step]
#directory= D:/Ievgen Ell/Measurements/NanoDiscs/MM_transmission/OD5
#mf is merged file
#plotMM, boolean, If you want to get polar plots of the measured MM elements put True
# Diffdecomp, boolean, if you want differential decomposition put True
#LuChipman, boolean, if you want polar decomposition yielding diatennuation,retardance and dipolarization put True
#Clode, boolean, to retrieve Clode matrix eigen values put True


def mergeFiles(angles,directory,filename,wvlcutofmin,wvlcutofmax,plotMM,diffDecomp,LuChipman,Clode):
    header=('AOI'+'\t'+'AZI'+'\t'+'Wvl'+'\t'+'MMe'+'\t'+'Value'+'\n'+'\n')
    muellerNames = getMuellerMatrixElementNames()
    #mf=open(directory+'/'+filename+'_merged.txt','w')
    #mf.write(header)
    if LuChipman:
        M_d={}
        M_r={}
        M_D={}
        #MMr={}
    if Clode:
        CDM1={}
        CDM2={}
        CDM3={}
        CDM4={}
        eigValH={}
    
    muellerMatrices = dict()
    for i in range(int((angles[1]-angles[0])/angles[2]+1)):     #Azimuthal angles are determined from file names
        print ("processing Azi %s"%(angles[0]+angles[2]*i)) 
        with open(directory+'/'+filename+str(angles[0]+angles[2]*i)+'.dat','r') as f:    # Azimuthal angles are determined from file names
            
            for lines in f:
            
                if lines[0]=='m':                                   # skip first non mm-values
                              
                    lines = re.split(r'\t+', lines)                    #formating the lines
                
                    #mf.write(str(lines[2]) +'\t'+ str((angles[0]+angles[2]*i))+'\t'+ str((lines[1]))+'\t'+str(lines[0])+'\t'+ str((lines[3])+'\n'))
                
                    AOI=lines[2]
                    if AOI not in muellerMatrices.keys():
                        muellerMatrices[AOI]=dict()
                        if LuChipman:
                             M_d[AOI]={}
                             M_r[AOI]={}
                             M_D[AOI]={}
                             #MMr[AOI]={}
                        if Clode:
                            CDM1[AOI]={}
                            CDM2[AOI]={}
                            CDM3[AOI]={}
                            CDM4[AOI]={}
                            eigValH[AOI]={}
                        
                    
                    AZI=float(angles[0]+angles[2]*i)
                    if AZI not in muellerMatrices[AOI].keys():
                        muellerMatrices[AOI][AZI]=dict()
                        if LuChipman:
                             M_d[AOI][AZI]={}
                             M_r[AOI][AZI]={}
                             M_D[AOI][AZI]={}
                             #MMr[AOI][AZI]={}
                        if Clode:
                            CDM1[AOI][AZI]={}
                            CDM2[AOI][AZI]={}
                            CDM3[AOI][AZI]={}
                            CDM4[AOI][AZI]={}
                            eigValH[AOI][AZI]={}
                        
                    wvl=float(lines[1])                                  
                    if wvl>=wvlcutofmin and wvl<=wvlcutofmax and wvl not in muellerMatrices[AOI][AZI].keys():
                        muellerMatrices[AOI][AZI][wvl]=list()
                        if LuChipman:
                             M_d[AOI][AZI][wvl]=[]
                             M_r[AOI][AZI][wvl]=[]
                             M_D[AOI][AZI][wvl]=[]
                             #MMr[AOI][AZI][wvl]=[]
                        if Clode:
                            CDM1[AOI][AZI][wvl]=[]
                            CDM2[AOI][AZI][wvl]=[]
                            CDM3[AOI][AZI][wvl]=[]
                            CDM4[AOI][AZI][wvl]=[]
                            eigValH[AOI][AZI][wvl]=[]
                        
               
             
                    
                    if lines[0] in muellerNames and wvl<=wvlcutofmax and wvl>=wvlcutofmin:
                        muellerMatrices[AOI][AZI][wvl].append(float(lines[3])) 
                
                
                    
                   
    #for AOI in sorted(muellerMatrices.keys()):
        #for azi in sorted(muellerMatrices[AOI].keys()):
            #for wl in sorted(muellerMatrices[AOI][azi].keys()):
                #mm = np.reshape(np.asarray(muellerMatrices[AOI][azi][wl]), (4,4)) 
                #if checkDepolarizingMuellerMatrix(mm):
                    #muellerMatrices[AOI][azi][wl] = np.copy(mm)
                #else:
                    #muellerMatrices[AOI][azi][wl] = np.copy(cloudeDecomposition(mm))
   


    for AOI in sorted(muellerMatrices.keys()):
        for azi in sorted(muellerMatrices[AOI].keys()):
            for wl in sorted(muellerMatrices[AOI][azi].keys()):
                mm = np.reshape(np.asarray(muellerMatrices[AOI][azi][wl]), (4,4))
                muellerMatrices[AOI][azi][wl] = np.copy(mm)
                if LuChipman:
                    M_d[AOI][azi][wl],M_r[AOI][azi][wl],M_D[AOI][azi][wl]=PolarDecomposition(mm)
                    #MMr[AOI][azi][wl]=M_r[AOI][azi][wl]@M_d[AOI][azi][wl]@M_D[AOI][azi][wl]
                if Clode:
                    MM1,MM2,MM3,MM4,eigVal=cloudeDecomposition(mm)
                    eigValH[AOI][azi][wl]=list(eigVal)
                    CDM1[AOI][azi][wl]=(np.reshape(np.asarray(MM1), (4,4)))
                    CDM2[AOI][azi][wl]=(np.reshape(np.asarray(MM2), (4,4)))
                    CDM3[AOI][azi][wl]=(np.reshape(np.asarray(MM3), (4,4)))
                    CDM4[AOI][azi][wl]=(np.reshape(np.asarray(MM4), (4,4)))
                
    if plotMM:
        createMuellerMatrixPlot(muellerMatrices)
    
    if diffDecomp:
        L,L_m,L_u=differentialDecomposition(muellerMatrices)
        createMuellerMatrixPlotL_m(L_m)
        createMuellerMatrixPlotL_u(L_u)
        
        
    
    if LuChipman:
        createMuellerMatrixPlotM_d(M_d)
        createMuellerMatrixPlotM_r(M_r)
        createMuellerMatrixPlotM_D(M_D)
    
    
    if Clode:
        plt.figure(1)
        lambda1=[]
        lambda2=[]
        lambda3=[]
        lambda4=[]
        AOI=list(eigValH.keys())[0]
        azi=sorted(list(eigValH[AOI].keys()))[9]
        wvls=list(eigValH[AOI][azi].keys())
        for wvl in sorted(list( eigValH[AOI][azi].keys())):
            lambda1.append(eigValH[AOI][azi][wvl][0])
            lambda2.append(eigValH[AOI][azi][wvl][1])
            lambda3.append(eigValH[AOI][azi][wvl][2])
            lambda4.append(eigValH[AOI][azi][wvl][3])
            
        plt.plot(wvls,lambda1,label="lambda1", color='r')
        plt.plot(wvls,lambda2,label="lambda2", color='g')
        plt.plot(wvls,lambda3,label="lambda3", color='b')
        plt.plot(wvls,lambda4,label="lambda4", color='k')
        plt.legend()
        plt.title('Cloude eigenvalues at AOI-%s, Azi-%s'%(AOI,azi))
        plt.show()
        #createMuellerMatrixPlot(CDM1)
        #createMuellerMatrixPlot(CDM2)
        #createMuellerMatrixPlot(CDM3)
        #createMuellerMatrixPlot(CDM4)
        #print( eigValH[0],eigValH[1],eigValH[2],eigValH[3])
    
        
        
    
   
    #mf.close()
    
    return 


#Program GUI

top = Tk()
top.title('Welcome to the Mueller Matrix analysis tool ver.0.1')

directory=StringVar()
fp = Label(top, text="File path")
fp.pack( side = TOP)
fpe = Entry(top,textvariable=directory,width=75, bd =5)
fpe.pack(side = TOP)





filename=StringVar()
fn = Label(top, text="File name")
fn.pack( side = TOP)
fne = Entry(top,textvariable=filename,width=75, bd =5)
fne.pack(side = TOP)




wvlcutofmin=DoubleVar()
mwl = Label(top, text="Minimal Wavelength")
mwl.pack( side = TOP)
mwle = Entry(top,textvariable=wvlcutofmin,width=6, bd =5)
mwle.pack(side = TOP)




wvlcutofmax=DoubleVar()
mawl = Label(top, text="Maximal Wavelength")
mawl.pack( side = TOP)
mawle = Entry(top,textvariable=wvlcutofmax,width=6, bd =5)
mawle.pack(side = TOP)








pMMe = BooleanVar()
Checkbutton(top, text="Plot MM", variable=pMMe).pack(side=TOP)
Label(top, text="Decomposition methods::").pack(side=TOP)
Cloudee = BooleanVar()
Checkbutton(top, text="Cloude decomposition", variable=Cloudee).pack(side=LEFT)
diffDecompe = BooleanVar()
Checkbutton(top, text="Differential decomposition", variable=diffDecompe).pack(side=LEFT)
LuChipmane = BooleanVar()
Checkbutton(top, text="Lu-Chipman decomposition", variable=LuChipmane).pack(side=LEFT)
mfe=BooleanVar()
Checkbutton(top, text="Several files", variable=mfe).pack(side=LEFT)

    

def close_window(): 
    top.destroy()

Start = Button(top, text="START", command=close_window)
Start.pack(side=BOTTOM)

top.mainloop()

directory=directory.get()
filename=filename.get()
wvlcutofmin=wvlcutofmin.get()
wvlcutofmax=wvlcutofmax.get()
plotMM=pMMe.get()
diffDecomp=diffDecompe.get()
LuChipman=LuChipmane.get()
Clode=Cloudee.get()
mf=mfe.get()

if mf:
     mergeFiles(angles,directory,filename,wvlcutofmin,wvlcutofmax,plotMM,diffDecomp,LuChipman,Clode)
else:
    oneFile(directory,filename,wvlcutofmin,wvlcutofmax,plotMM,diffDecomp,LuChipman,Clode)
    






