###Program to calculate the fits and everything for a full NOESY run

#import section
import os
import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from modules.NOESY2_custom import *


###READING the initial files. In those files you should note your lipid peaks and  chemical ppm positions
###READ THE LIPID FILE
lipid_data=open_file_sing('ini_data/lipid.txt')
###READ THE CHEMICAL FILE
chemical_data=open_file_sing('ini_data/chemical.txt')

###MAKING the directories to store the data in
dir="zresults"
if not os.path.exists(dir):
	os.mkdir(dir)
result_path=dir+"/"+chemical_data[0]
if not os.path.exists(result_path):
	os.mkdir(result_path)
for element in range(1,len(chemical_data)):
	if not os.path.exists(result_path+"/"+str(chemical_data[element])+"ppm"):
		os.mkdir(result_path+"/"+str(chemical_data[element])+"ppm")


###IMPORTANT DATA
#This is the location in the list of your integral data, for example:
#Integral 50	45684168	otherdata	otherdata
#here your position would be 1
integral_location=1

lipid_len=len(lipid_data)-1
chemical_len=len(chemical_data)-1

#Define the mixing times
t_m=[0,0.1,0.2,0.3,0.5]


###SELECT THE CORRESPONDING FILES, select your files with integral data from topspin
#---------------------------------------------------------------------------------------------------------------
print("Select 0,1ms")
filepaths=['01ms.txt','100ms.txt','200ms.txt','300ms.txt','500ms.txt']
#filepaths=import_files()
N01ms=open_file_mult(filepaths[0])
del N01ms[0]

print("Select 100ms")
#filepaths=import_files()
N100ms=open_file_mult(filepaths[1])
del N100ms[0]

print("Select 200ms")
#filepaths=import_files()
N200ms=open_file_mult(filepaths[2])
del N200ms[0]

print("Select 300ms")
#filepaths=import_files()
N300ms=open_file_mult(filepaths[3])
del N300ms[0]

print("Select 500ms")
#filepaths=import_files()
N500ms=open_file_mult(filepaths[4])
del N500ms[0]
#---------------------------------------------------------------------------------------------------------------



###MAKING THE FINALISED FILES
###In this section the data is rearranged and put into easily imported files for Origin
#Final Format:
#Time(s)	Peak1	Peak2	etc.
#0			Int1
#0.1		Int2
#and so on
#---------------------------------------------------------------------------------------------------------------

###0.1 mixing time, saving the crosspeaks
dat=open(result_path+"/"+'crosspeaks.txt',"w")
for peak_lipid in range(1,len(lipid_data)):
	dat.write(str(lipid_data[peak_lipid]))
	dat.write("	")
dat.write("\n")

for peak in range(len(N01ms)):
	dat.write(N01ms[peak][integral_location])
	dat.write("	")
dat.write("\n")
dat.close()

###100-500ms mxixing times
for peak_chemical in range(len(chemical_data)-1):
	#Open the corresponding file
	dat=open(result_path+"/"+chemical_data[peak_chemical+1]+"ppm"+"/"+chemical_data[peak_chemical+1]+'ppm.txt',"w")
	
	#Writing the First Line
	dat.write('Zeit(s)')
	
	###Writing the Header line for Origin
	#Writing in the LIPID peaks
	for peak_lipid in range(1,len(lipid_data)):
		dat.write("	")
		dat.write(lipid_data[peak_lipid])
	dat.write("\n")
	
	#Getting the zero datapoint
	dat.write("0")
	for i in range(len(lipid_data)-1):
		dat.write("	")
		dat.write("0")
	dat.write("\n")
	
	###Writing the 100-500ms Signals into the files, the files can the be imported into ORIGIN
	#Writing in the INTEGRALS
	#100ms
	dat.write("0.1")
	for peak_lipid in range(len(lipid_data)-1):
		dat.write("	")
		dat.write(N100ms[peak_chemical*10+peak_lipid][integral_location])
	dat.write("\n")
	
	#200ms
	dat.write("0.2")
	for peak_lipid in range(len(lipid_data)-1):
		dat.write("	")
		dat.write(N200ms[peak_chemical*10+peak_lipid][integral_location])
	dat.write("\n")
	
	#300ms
	dat.write("0.3")
	for peak_lipid in range(len(lipid_data)-1):
		dat.write("	")
		dat.write(N300ms[peak_chemical*10+peak_lipid][integral_location])
	dat.write("\n")
	
	#500ms
	dat.write("0.5")
	for peak_lipid in range(len(lipid_data)-1):
		dat.write("	")
		dat.write(N500ms[peak_chemical*10+peak_lipid][integral_location])
	dat.write("\n")
	
	#Closing
	dat.close()
#---------------------------------------------------------------------------------------------------------------



###CALULATING THE MIXING TIME AND THE CROSSRELAXATION
#---------------------------------------------------------------------------------------------------------------
#Here the mixing time and crossrelaxation are calculated
#as in Origin the used algorithm is the LEVENBERG-MARQUARDT
#we won't reopen the files, because we still have all our data saved in our arrays

##REMINDER
#lists that are in use and important:
#N01ms,N100ms,N200ms,N300ms,N500ms
#chemical_data,lipid_data

###FUNCTION THAT the data will be fitted to:
#A_is=(A_ii_0/2)*(1-exp(-2*sigma*t_m))*exp(-t_m/T_is)
#known parameters:
#A_ii_0 - the integral of our diagonal crosspeaks at t=0.1, this we get from N01ms
#t_m is our mixing time, which will be either 100,200,300 or 500 ms
#sigma and T_is are to be calculated

#Defining the two array sigma_array and T_is_array where the values that are calculated are stored into for later use
#Format is x-Lipids y-chemical
sigma_array=[]
T_is_array=[]
sigma_error_array=[]

for peak_chemical in range(len(chemical_data)-1):
	print("CACLULATING CHEMICAL PEAK "+str(chemical_data[peak_chemical+1]))
	#Setting up empty lists to save the sigma and T_is into
	sigma_list=[]
	T_is_list=[]
	#Setting up the error lists
	sigma_error_list=[]
	
	##Checking if there is an elimination file and if not then creating it
	#-------------------------------
	if not os.path.isfile(dir+"/"+chemical_data[0]+"/"+str(chemical_data[peak_chemical+1])+"ppm"+"/"+str(chemical_data[peak_chemical+1])+"_eliminator.txt"):
		dat=open(dir+"/"+chemical_data[0]+"/"+str(chemical_data[peak_chemical+1])+"ppm"+"/"+str(chemical_data[peak_chemical+1])+"_eliminator.txt","w")
		dat.write("Zeit(s)")
		for peak_lipid in range(len(lipid_data)-1):
			dat.write("	")
			dat.write(str(lipid_data[peak_lipid+1]))
		dat.write("\n")
	
		for t in range(len(t_m)):
			dat.write(str(t_m[t]))
			for peak_lipid in range(len(lipid_data)-1):
				dat.write("	")
				dat.write(str(1))
			dat.write("\n")
		dat.close()
	#-------------------------------
	
	##Applying the elimination file
	#Open the eliminator file and cutting off the string part
	eliminator=open_file_mult(dir+"/"+chemical_data[0]+"/"+str(chemical_data[peak_chemical+1])+"ppm"+"/"+str(chemical_data[peak_chemical+1])+"_eliminator.txt")
	del eliminator[0]
	eliminator=list(np.array(eliminator).transpose())
	del eliminator[0]
	eliminator=np.array(eliminator).transpose()
	#Destringing
	eliminator = [[float(j) for j in i] for i in eliminator]
	eliminator=np.array(eliminator)

	
	for peak_lipid in range(len(lipid_data)-1):
		print(str(peak_lipid+1)+"/"+str(len(lipid_data)-1))
		#Make the A_is list for the function solver
		A_is=[]
		A_is.append(0)
		A_is.append(float(N100ms[peak_chemical*10+peak_lipid][integral_location]))
		A_is.append(float(N200ms[peak_chemical*10+peak_lipid][integral_location]))
		A_is.append(float(N300ms[peak_chemical*10+peak_lipid][integral_location]))
		A_is.append(float(N500ms[peak_chemical*10+peak_lipid][integral_location]))
		
		#Define the Integral at t=0
		A_ii_0=float(N01ms[peak_lipid][integral_location])
		
		#IN THIS PART THE ELIMINATION TAKES PLACE
		#Only values that are assigned a 1 in the elimination file are assigned to the lists whcih get calculated
		#For less than four values no caluclationas are made
		t_m_calc=[]
		A_is_calc=[]
		for t in range(len(t_m)):
			if eliminator[t][peak_lipid]==1:
				t_m_calc.append(t_m[t])
				A_is_calc.append(A_is[t])
		
		if len(t_m_calc)>3:
			if (len(t_m_calc)==4):
				#SOLVING 4 EQUATIIONS
				start_values=[[0.01,0.4],[1,1],[0.1,0.4],[0.1,7],[5,0.5],[0.5,5],[0.2,0.9],[0.27,0.42]]
				#Initial Solving of the fit, LEVENBERG MARQUARDT
				sol=root(levenberg4,[0.001,0.4],args=(A_is_calc,t_m_calc,A_ii_0),jac=True,method='lm')
				ini_error=np.sum((sol.fun-np.array(A_is_calc))**2)
				#Solving everything and checking for a better solution, LEVENBERG MARQUARDT
				for val in start_values:
					fsol=root(levenberg4,val,args=(A_is_calc,t_m_calc,A_ii_0),jac=True,method='lm')
					err=np.sum((fsol.fun-np.array(A_is_calc))**2)
					if err<ini_error:
						sol=fsol
						ini_error=err
				#The n for the RSS correction later on
				nrss=4 
				
			elif (len(t_m_calc)==5):
				#SOLVING 5 EQUATIONS
				start_values=[[0.01,0.4],[1,1],[0.1,0.4],[0.1,7],[5,0.5],[0.5,5],[0.2,0.9],[0.27,0.42]]
				#Initial Solving of the fit, LEVENBERG MARQUARDT
				sol=root(levenberg5,[0.001,0.4],args=(A_is_calc,t_m_calc,A_ii_0),jac=True,method='lm')
				ini_error=np.sum((sol.fun-np.array(A_is_calc))**2)
				#Solving everything and checking for a better solution, LEVENBERG MARQUARDT
				for val in start_values:
					fsol=root(levenberg5,val,args=(A_is_calc,t_m_calc,A_ii_0),jac=True,method='lm')
					err=np.sum((fsol.fun-np.array(A_is_calc))**2)
					if err<ini_error:
						sol=fsol
						ini_error=err
				#The n for the RSS correction later on
				nrss=5
				
			#WITHDRAW THE ARGUMENTS
			print(sol)
		
			#DRAW THE ARGUMENTS FROM THE SOLUTION
			sigma=sol.x[0]
			T_is=sol.x[1]
		
			#DRAW THE ERRORS FROM THE SOLUTION
			sigma_error=np.sqrt(np.diag(sol.cov_x))[0]
			print(str(sigma_error))
			#CORRECT THE ERROR TO GET FROM STANDARD DEVIATION TO STANDARD ERROR AS ORIGIN DOES
			#Short explanation here for later or for whoever has to modify this code(my sincerest condolences to you)
			#The covariance matrix apparently has the RSS-the residual sum of squares multiplied into it
			#that means we have to edit it out so our error is not too small
			#you can also look that up at https://www.originlab.com/doc/Origin-Help/NLFit-Theory#Parameter_Standard_Errors
			#there they explain it; so here we calculate our RSS
			RSS=0
			for i in range(len(t_m_calc)):
				RSS+=((A_ii_0/float(2))*(1-math.exp(-2*sigma*t_m_calc[i]))*math.exp(-t_m_calc[i]/float(T_is))-A_is_calc[i])**2
			#Then our RSS gets divided by n=Number of data points that are used for the fit, and
			#p=number of parameters that the fit has to calculate
			#so: 1 unknown->p=1	2 unknown->p=2
			#then the RSS is multiplied with the sqrt of the diagonal element in the covariance matrix and we get our STANDARD ERROR
			RSS=RSS/float(nrss-2)
			print(RSS)
			sigma_error=sigma_error*float(np.sqrt(RSS))
			print(str(sigma_error))
		
		elif (len(t_m_calc)<=3):
			#This part is true when so many datapoints have been eliminated that no analysis is fruitful anymore
			#The sigmas are here eliminated
			sigma=0
			sigma_error=0
			T_is=0
		
		
		#Save the arguments into lists to cleanly append to array
		sigma_list.append(sigma)
		T_is_list.append(T_is)
		sigma_error_list.append(sigma_error)

		
		
		#PLOT THE POINTS OF DATA
		plot_fit_curves(sigma,sigma_error,T_is,A_is_calc,t_m_calc,A_ii_0,peak_chemical,peak_lipid,chemical_data,lipid_data,dir)


	#SAVING THE CALULATED VALUES INTO AN ARRAY TO BE USED LATER
	sigma_array.append(sigma_list)
	T_is_array.append(T_is_list)
	sigma_error_array.append(sigma_error_list)
		


#---------------------------------------------------------------------------------------------------------------




###SAVING THE CALCULATED VALUES INTO A .txt FILE
#---------------------------------------------------------------------------------------------------------------
dat=open(result_path+"/"+"sigma"+'.txt',"w")
dat.write(str(chemical_data[0])+"/"+str(lipid_data[0]))
for peak_lipid in range(1,len(lipid_data)):
	dat.write("	")
	dat.write(str(lipid_data[peak_lipid]))
dat.write("\n")

for peak_chemical in range(1,len(chemical_data)):
	dat.write(str(chemical_data[peak_chemical]))
	for peak_lipid in range(1,len(lipid_data)):
		dat.write("	")
		dat.write(str(sigma_array[peak_chemical-1][peak_lipid-1]))
	dat.write("\n")
dat.close()


dat=open(result_path+"/"+"sigma_error"+'.txt',"w")
dat.write(str(chemical_data[0])+"/"+str(lipid_data[0]))
for peak_lipid in range(1,len(lipid_data)):
	dat.write("	")
	dat.write(str(lipid_data[peak_lipid]))
dat.write("\n")

for peak_chemical in range(1,len(chemical_data)):
	dat.write(str(chemical_data[peak_chemical]))
	for peak_lipid in range(1,len(lipid_data)):
		dat.write("	")
		dat.write(str(sigma_error_array[peak_chemical-1][peak_lipid-1]))
	dat.write("\n")
dat.close()


dat=open(result_path+"/"+"T_is"+'.txt',"w")
dat.write(str(chemical_data[0])+"/"+str(lipid_data[0]))
for peak_lipid in range(1,len(lipid_data)):
	dat.write("	")
	dat.write(str(lipid_data[peak_lipid]))
dat.write("\n")

for peak_chemical in range(1,len(chemical_data)):
	dat.write(str(chemical_data[peak_chemical]))
	for peak_lipid in range(1,len(lipid_data)):
		dat.write("	")
		dat.write(str(T_is_array[peak_chemical-1][peak_lipid-1]))
	dat.write("\n")
dat.close()
#---------------------------------------------------------------------------------------------------------------





###REARRANGING THE LIPID PEAKS FOR THE DESIRED OUTPUT
#---------------------------------------------------------------------------------------------------------------
#Importing the desired output file
lipid_sorted=open_file_sing('plot_data/lipid_sorted.txt')

#Creating an empty sigma array to store the rearranged data in
sigma_new = [[] for i in range(len(lipid_data)-1)]
sigma_error_new = [[] for i in range(len(lipid_data)-1)]

sigma_old_transposed=np.array(sigma_array).transpose()
sigma_error_old_transposed=np.array(sigma_error_array).transpose()

#Rearranging the sigma-array
for peak_lipid in range(1,len(lipid_data)):
	ni=lipid_sorted.index(lipid_data[peak_lipid])
	sigma_new[ni-1]=sigma_old_transposed[peak_lipid-1]
	sigma_error_new[ni-1]=sigma_error_old_transposed[peak_lipid-1]

sigma_new=list(np.array(sigma_new).transpose())
sigma_error_new=list(np.array(sigma_error_new).transpose())
#---------------------------------------------------------------------------------------------------------------


###PLOTTING THE CROSSCORRELATION SIGMA TO THE RELATIVE PEAKS
#---------------------------------------------------------------------------------------------------------------
#Importing the plot data
plot_data=open_file_sing('plot_data/plot_data.txt')

###PEAK ELIMINATION PART
#Here we check if there is a file "eliminating" certain peaks and if not we generate one
#Of course if there is such a file then we eliminate the required peaks
if not os.path.isfile(result_path+"/barplot_elimination.txt"):
	dat=open(result_path+"/"+"barplot_elimination"+'.txt',"w")
	
	dat.write(str(chemical_data[0])+"/"+str(lipid_data[0]))
	for peak_lipid in range(1,len(lipid_sorted)):
		dat.write("	")
		dat.write(str(lipid_sorted[peak_lipid]))
	dat.write("\n")

	for peak_chemical in range(1,len(chemical_data)):
		dat.write(str(chemical_data[peak_chemical]))
		for peak_lipid in range(1,len(lipid_data)):
			dat.write("	")
			dat.write(str(1))
		dat.write("\n")
	dat.close()

##Eliminating those unwanted peaks
#Open the eliminator file and cutting off the string part
eliminator_bar=open_file_mult(result_path+"/barplot_elimination.txt")
del eliminator_bar[0]
eliminator_bar=list(np.array(eliminator_bar).transpose())
del eliminator_bar[0]
eliminator_bar=np.array(eliminator_bar).transpose()
#Destringing
eliminator_bar = [[float(j) for j in i] for i in eliminator_bar]
eliminator_bar=np.array(eliminator_bar)

#Now we are ready for ELIMINATION
sigma_new=sigma_new*eliminator_bar
sigma_error_new=sigma_error_new*eliminator_bar


###PLOTTING
for peak_chemical in range(1,len(chemical_data)):
	print("Plot sigma="+str(chemical_data[peak_chemical]))
	plot_cross_correlation(peak_chemical,chemical_data,lipid_sorted,sigma_new,sigma_error_new,dir,plot_data)

#---------------------------------------------------------------------------------------------------------------












