###Program with custom functions for 31P analyser
import _tkinter
import tkinter as tk
import tkinter.filedialog as filedialog
import re
import os
import sys
import math
import cmath
import numpy as np
from scipy.optimize import root
import matplotlib.pyplot as plt



#Importing all files
def import_files():
	root = tk.Tk()
	root.withdraw()
	file_paths = filedialog.askopenfilenames(parent=root,title='Select all files:')
	return file_paths

###SINGLE LINE OF DATA
#Opening a file and writing its contents into a list
def open_file_sing(path):
	filelist=[]
	with open(str(path)) as file:
		for line in file:
			filelist.append(line.strip())
	file.close()
	return filelist

###MULTIPLE LINES OF DATA
#Opening a single file and writing its contents into a list of lists
def open_file_mult(path):
	filelol=[]
	with open(str(path)) as file:
		for line in file:
			filelol.append(line.strip().split('	'))
	file.close()
	return filelol

	
###THE LEVENBERG-MARQUARDT ALGORITHM
#The LM algorithm is used to solve the set of differential equations
#NOTE: x in this case is going to become a vector in which x[0]^=sigma x[1]^=T_is
def levenberg4(x,A_is,t_m,A_ii_0):
	#This is our original function
	f=[
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[0]))*math.exp(-t_m[0]/float(x[1]))-A_is[0],
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[1]))*math.exp(-t_m[1]/float(x[1]))-A_is[1],
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[2]))*math.exp(-t_m[2]/float(x[1]))-A_is[2],
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[3]))*math.exp(-t_m[3]/float(x[1]))-A_is[3]
	]
	
	#This is our derivation for each parameter
	df=np.array([
	[(2*t_m[0])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[0])*math.exp(-t_m[0]/float(x[1])),
	((A_ii_0*t_m[0])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[0]))*math.exp(-t_m[0]/float(x[1]))],
	
	[(2*t_m[1])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[1])*math.exp(-t_m[1]/float(x[1])),
	((A_ii_0*t_m[1])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[1]))*math.exp(-t_m[1]/float(x[1]))],
	
	[(2*t_m[2])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[2])*math.exp(-t_m[2]/float(x[1])),
	((A_ii_0*t_m[2])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[2]))*math.exp(-t_m[2]/float(x[1]))],
	
	[(2*t_m[3])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[3])*math.exp(-t_m[3]/float(x[1])),
	((A_ii_0*t_m[3])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[3]))*math.exp(-t_m[3]/float(x[1]))]
	])
	
	return f, df

def levenberg5(x,A_is,t_m,A_ii_0):
	#This is our original function
	f=[
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[0]))*math.exp(-t_m[0]/float(x[1]))-A_is[0],
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[1]))*math.exp(-t_m[1]/float(x[1]))-A_is[1],
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[2]))*math.exp(-t_m[2]/float(x[1]))-A_is[2],
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[3]))*math.exp(-t_m[3]/float(x[1]))-A_is[3],
	(A_ii_0/float(2))*(1-math.exp(-2*x[0]*t_m[4]))*math.exp(-t_m[4]/float(x[1]))-A_is[4]
	]
	
	#This is our derivation for each parameter
	df=np.array([
	[(2*t_m[0])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[0])*math.exp(-t_m[0]/float(x[1])),
	((A_ii_0*t_m[0])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[0]))*math.exp(-t_m[0]/float(x[1]))],
	
	[(2*t_m[1])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[1])*math.exp(-t_m[1]/float(x[1])),
	((A_ii_0*t_m[1])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[1]))*math.exp(-t_m[1]/float(x[1]))],
	
	[(2*t_m[2])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[2])*math.exp(-t_m[2]/float(x[1])),
	((A_ii_0*t_m[2])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[2]))*math.exp(-t_m[2]/float(x[1]))],
	
	[(2*t_m[3])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[3])*math.exp(-t_m[3]/float(x[1])),
	((A_ii_0*t_m[3])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[3]))*math.exp(-t_m[3]/float(x[1]))],
	
	[(2*t_m[4])*(A_ii_0/float(2))*math.exp(-2*x[0]*t_m[4])*math.exp(-t_m[4]/float(x[1])),
	((A_ii_0*t_m[4])/float(2*x[1]**2))*(1-math.exp(-2*x[0]*t_m[4]))*math.exp(-t_m[4]/float(x[1]))]
	])
	
	return f, df

###THE PLOTTING OF THE FITTED CURVES WITH THE LM ALGORITHM
def plot_fit_curves(sigma,sigma_error,T_is,A_is,t_m,A_ii_0,peak_chemical,peak_lipid,chemical_data,lipid_data,dir):
	#Plotting the Integral points
	for i in range(len(A_is)):
		x=plt.plot(t_m[i],A_is[i],marker="s",color="#000000",linestyle="None")
	#Plotting the fit
	if len(t_m)>3:
		t_m_fine=np.linspace(0,0.5,num=1000)
		A_is_fine=[]
		for k in t_m_fine:
			a=(A_ii_0/float(2))*(1-math.exp(-2*sigma*k))*math.exp(-k/float(T_is))
			A_is_fine.append(a)
		y,=plt.plot(t_m_fine,A_is_fine,color='#ff5733',label="Fit")
		#Generating the legend
		plt.legend(handles=[y],loc=1)
	#Plotting the paramters
	plt.text(0.3, 0.9*max(A_is),'$\sigma$='+str(round(sigma,4)),verticalalignment='top', horizontalalignment='left', fontsize=14)
	plt.text(0.3, 0.85*max(A_is),'$\sigma_{err}$='+str(round(sigma_error,4)),verticalalignment='top', horizontalalignment='left', fontsize=14)
	plt.text(0.3, 0.80*max(A_is),'$T_{is}$='+str(round(T_is,4)),verticalalignment='top', horizontalalignment='left', fontsize=14)
	#Title of the plot
	plt.title('Fit of the '+str(chemical_data[peak_chemical+1])+'at '+str(lipid_data[peak_lipid+1]))
	#Labeling the axes
	plt.xlabel('Time [s]')
	plt.ylabel(str(lipid_data[peak_lipid+1]))
	#Saving the plot
	plt.savefig(dir+"/"+chemical_data[0]+"/"+str(chemical_data[peak_chemical+1])+"ppm"+"/"+str(chemical_data[peak_chemical+1])+"_"+str(lipid_data[peak_lipid+1])+".png")
	plt.clf()
	plt.close()


###THE PLOTTING OF THE CROSS RELAXATION RATES
def plot_cross_correlation(peak_chemical,chemical_data,lipid_sorted,sigma_new,sigma_error_new,dir,plot_data):
	#Plotting the bar graphs 
	x=range(len(lipid_sorted)-1)
	y=sigma_new[peak_chemical-1]
	y_err=sigma_error_new[peak_chemical-1]
	star=[]
	for s in y:
		if s==0:
			star.append(0.05*(np.max(sigma_new[peak_chemical-1])+max(y_err))*1.05)
		else:
			star.append(None)
	plt.scatter(x,star,marker="*",c="k")
	plt.bar(x, y,color=plot_data[1],linewidth="0.2",edgecolor=plot_data[3],yerr=y_err)
	#Adding the ticks
	plt.xticks(x, lipid_sorted[1::],rotation=45)
	#Scaling the y-Axis
	#i,j = np.where(sigma_new == np.max(sigma_new))
	#plt.ylim(0,(np.max(sigma_new[peak_chemical-1])+sigma_error_new[i[0]][j[0]])*1.05)
	plt.ylim(0,(np.max(sigma_new[peak_chemical-1])+max(y_err))*1.05)
	#Title of the plot
	plt.title('Cross relaxation $\sigma$ at '+str(chemical_data[peak_chemical])+'ppm')
	#Labeling the axes
	plt.ylabel('Cross relaxation $\sigma$')
	plt.xlabel(" ")
	###SIZES OF PLOTS
	fig = plt.gcf()
	fig.set_size_inches(6, 4.5)
	fig.set_dpi(300)
	#Adjust the bottom so that the ticks don't get cut off
	plt.gcf().subplots_adjust(bottom=0.15)
	#Saving the plot
	plt.savefig(dir+"/"+chemical_data[0]+"/"+str(chemical_data[0])+"_sigma_"+str(chemical_data[peak_chemical])+"ppm.png")
	plt.clf()
	plt.close()

'''
###MAKING A .txt FILE
def make_dat():
	parlis=["LB1","LB2","LB3","csa","axi","hex","iso","c","a","drift","g"]
	dat=open('parameters'+'.txt',"w")
	for par in parlis:
		dat.write(str(par))
		dat.write("	")
	dat.write("\n")
	dat.close()

###WRITING DATA INTO .txt FILE
def write_to_dat(parameter_list,g_sd):
	dat=open('parameters'+'.txt',"a")
	for par in parameter_list:
		dat.write(str(par))
		dat.write("	")
	dat.write(str(g_sd))
	dat.write("\n")
	dat.close()	
'''