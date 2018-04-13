#! /usr/bin/env python2
# -*- coding: utf-8 -*-
__author__="Roxane Lassis"
__version__="0.1"
__email__="lassis@etud.insa-toulouse.fr"
'''
    File name: GA_Main.py
    Author: Roxane Lassis 
    Description : Complete Test with Initialization and evolution of a population, Display and Saving of the results 
		- Subarray selection using Genetic Algorithm
                  with parameters given in a file '.txt'
    Context : ALMA internship
    Date created: 06/2017
    Date last modified: 08/2017
    Python Version: 2.7
    
    History:
    
    2018.04.10:
        - update for automation
        
    2018.04.13:
        - updating the output-
        
'''
import numpy as np
from pylab import *
import pandas as pd
import csv
import time
import sys
import os
import glob
import shutil
DOSSIER_COURRANT = os.path.dirname(os.path.abspath(__file__))
DOSSIER_PARENT = os.path.dirname(DOSSIER_COURRANT)
sys.path.append(os.path.join(DOSSIER_PARENT,"Codes"))
import classes as c
import evolution as ev
import Useful_Functions_for_GA_Main as uf
import compute_constraints_subarray as ccs

import pickle


#==================================================================================================================================
#==================================================================================================================================
# Import the inputs for the optimization : 
# - initialize a population thanks to a configuration file
# - evolve a population doing a genetic algorithm
# - save the subarrays solution at "casa" format
#!!! TO FULFILL : GA_Subarray_Selection/GA_Inputs.txt,  GA_Subarray_Selection/GA_Parameters.txt  !!!
#==================================================================================================================================
#==================================================================================================================================
File_Inputs='GA_Subarray_Selection/GA_Inputs_Automatic_pop.txt'
inputs= pd.read_csv(
        File_Inputs, comment='#', names=['NAMES','VALUES'], sep='\s+'+':'+'\s+')
Names=inputs[['NAMES']].values
Values=inputs[['VALUES']].values
for i in range(0,len(Names)):
	if Names[i]=='File_Parameters':
		File_Parameters=Values[i][0]
	elif Names[i]=='Display_Screen_Results':
		if Values[i][0]=="False":
    			Display_Screen_Results=False
		elif Values[i][0]=="True":
    			Display_Screen_Results=True  
		else:
    			raise Exception('Invalid input for Display_Screen_Results : should be "True" or "False"')
	elif Names[i]=='Save_Folder_Results':
		if Values[i][0]=="False":
    			Save_Folder_Results=False
		elif Values[i][0]=="True":
   			 Save_Folder_Results=True  
		else:
  			  raise Exception('Invalid input for Save_Folder_Results : should be "True" or "False"')
	elif Names[i]=='Launch_CASA_Simulation':
		if Values[i][0]=="False":
    			Launch_CASA_Simulation=False
		elif Values[i][0]=="True":
   			 Launch_CASA_Simulation=True
    			 if Save_Folder_Results==False :
        			raise Exception('Can NOT launch CASA Simulation because Save_Folder_Results should be "True"')
		else:
   			 raise Exception('Invalid input for Launch_CASA_Simulation : should be "True" or "False"')

	elif Names[i]=='Display_Figure_Score_Evolution':
		if Values[i][0]=="False":
    			Display_Figure_Score_Evolution=False
		elif Values[i][0]=="True":
    			Display_Figure_Score_Evolution=True  
		else:
  			  raise Exception('Invalid input for Display_Figure_Score_Evolution : should be "True" or "False"')
	else:
		print 'Invalid file '+File_Inputs+'. Check if there is all the inputs and if the names of the inputs are correct. See the example files in the folder GA_Subarray_Selection/Input_Files_To_Fulfill/Input_Files_Examples.'
	
#==================================================================================================================================
#==================================================================================================================================
# Import the inputs for the GA parameters : 
#!!! TO FULFILL : GA_Subarray_Selection/GA_Parameters.txt !!!
#==================================================================================================================================
#==================================================================================================================================


#==================================================================================================================================
# Configuration : 
#==================================================================================================================================
param = pd.read_csv(
        File_Parameters, comment='#', names=['NAMES', 'VALUES'], sep='\s+'+':'+'\s+')


cm=c.Configuration_Manager()
configuration_Manager=cm

Values=param[['VALUES']].values
Names=param[['NAMES']].values

list_num_pads_subarrays=[] 
list_subarrays_weights=[]
list_Number_of_Iterations=[]
list_Population_Size=[]
list_Termination_Condition=[]
list_Threshold=[]
list_Mutation_Rate=[]
list_Tournament_Size=[]
list_Number_for_Elitism=[]
for j in range(0,len(Names)):
	if Names[j]=='Configuration_Input_File':
		Configuration_Input_File=Values[j][0]
		if Configuration_Input_File=='Actual_Config':
			os.system("python GA_Subarray_Selection/CASA_Sim/arrayConfiguration_actual_config.py --aos")
			my_file=open("GA_Subarray_Selection/Cfg_Input_Files/Actual_Config.txt","r")
			Configuration_Input_File=my_file.read()
		num_pads=sum(1 for _ in open(str(Configuration_Input_File)))-3
	elif Names[j]=='Source_Declination':
		Source_Declination=float(Values[j][0])
	elif Names[j]=='Source_Hour_Angle':
		Source_Hour_Angle=float(Values[j][0])
	elif Names[j]=='Number_of_Subarrays':
		Number_of_Subarrays=int(Values[j][0])
		now=time.time()
		list_objective_constraints_res=[now]*Number_of_Subarrays
		list_objective_constraints_mrs=[now]*Number_of_Subarrays
		list_objective_constraints_e=[now]*Number_of_Subarrays
		list_objective_constraints_s=[now]*Number_of_Subarrays
              
		list_constraints_weights_res=[now]*Number_of_Subarrays
		list_constraints_weights_mrs=[now]*Number_of_Subarrays
		list_constraints_weights_e=[now]*Number_of_Subarrays
		list_constraints_weights_s=[now]*Number_of_Subarrays

		test_sum=[-0]*Number_of_Subarrays
	elif Names[j]=='Number_of_Pads_per_Subarray':
		param_list_num_pads_subarrays=Values[j][0]     
		param_list_num_pads_subarrays.split
		for i in param_list_num_pads_subarrays.split(","):
   			 list_num_pads_subarrays.append(int(i))
	elif Names[j]=='Weights_for_Subarrays':
		param_list_subarrays_weights=Values[j][0]
		list_int=[]
		s=0
		for i in param_list_subarrays_weights.split(","):
    			list_int.append(float(i))

		if len(list_int)!=Number_of_Subarrays:
    			raise Exception('Invalid number of subarrays_weights')
		for i in range(0,Number_of_Subarrays):
    			if list_int[i]<0:
 				raise Exception('Negative weight for subarray')
    			list_subarrays_weights.append(list_int[i])
    			s+=list_int[i]

		if s!=0 and s!=1:
  			for i in range(0,len(list_subarrays_weights)):
    				list_subarrays_weights[i]=round(list_subarrays_weights[i]/s,2)

#==================================================================================================================================
# Parameters for the genetic algorithm : 
#==================================================================================================================================
	
	elif Names[j]=='Number_of_Iterations':
		param_list_Number_of_Iterations=Values[j][0]     
		param_list_Number_of_Iterations.split
		for i in param_list_Number_of_Iterations.split(","):
   			 list_Number_of_Iterations.append(int(i))
	elif Names[j]=='Population_Size':
		param_list_Population_Size=Values[j][0]     
		param_list_Population_Size.split
		for i in param_list_Population_Size.split(","):
   			 list_Population_Size.append(int(i))
	elif Names[j]=='Termination_Condition':
		param_list_Termination_Condition=Values[j][0]     
		param_list_Termination_Condition.split
		for i in param_list_Termination_Condition.split(","):
			if i=="False":
    				list_Termination_Condition.append(i)
			elif i=="True":
    				list_Termination_Condition.append(i)
			else:
  			  	raise Exception('Invalid termination condition '+i+' : should be "True" or "False"')
	elif Names[j]=='Threshold':
		param_list_Threshold=Values[j][0]     
		param_list_Threshold.split
		for i in param_list_Threshold.split(","):
   			 list_Threshold.append(i)
	elif Names[j]=='Mutation_Rate':
		param_list_Mutation_Rate=Values[j][0]     
		param_list_Mutation_Rate.split
		for i in param_list_Mutation_Rate.split(","):
   			 list_Mutation_Rate.append(float(i))
	elif Names[j]=='Tournament_Size':
		param_list_Tournament_Size=Values[j][0]     
		param_list_Tournament_Size.split
		for i in param_list_Tournament_Size.split(","):
   			 list_Tournament_Size.append(int(i))
	elif Names[j]=='Number_for_Elitism':	
		param_list_Number_for_Elitism=Values[j][0]     
		param_list_Number_for_Elitism.split
		for i in param_list_Number_for_Elitism.split(","):
   			 list_Number_for_Elitism.append(int(i))
#==================================================================================================================================
# Constraints and weights : 
#==================================================================================================================================

	else:
		if int(Names[j][0][-1])>Number_of_Subarrays-1:
			raise Exception('Invalid number of subarray for the input .._Sub_'+str(int(Names[j][0][-1]))+'. Should be in [0,'+str(Number_of_Subarrays-1)+'].')
		elif Names[j][0][:-1]=='Spatial_Resolution_Sub_':
			list_objective_constraints_res[int(Names[j][0][-1])]=float(Values[j][0])
		elif Names[j][0][:-1]=='Maximum_Recoverable_Scale_Sub_':
			list_objective_constraints_mrs[int(Names[j][0][-1])]=float(Values[j][0])
		elif Names[j][0][:-1]=='Elongation_Sub_':
			list_objective_constraints_e[int(Names[j][0][-1])]=float(Values[j][0])
		elif Names[j][0][:-1]=='Percentage_of_Sidelobes_Sub_':
			list_objective_constraints_s[int(Names[j][0][-1])]=float(Values[j][0])
		elif Names[j][0][:-1]=='Weight_for_Spatial_Resolution_Sub_':
			list_constraints_weights_res[int(Names[j][0][-1])]=float(Values[j][0])
			test_sum[int(Names[j][0][-1])]+=list_constraints_weights_res[int(Names[j][0][-1])]
		elif Names[j][0][:-1]=='Weight_for_Maximum_Recoverable_Scale_Sub_':
			list_constraints_weights_mrs[int(Names[j][0][-1])]=float(Values[j][0])
			test_sum[int(Names[j][0][-1])]+=list_constraints_weights_mrs[int(Names[j][0][-1])]	
		elif Names[j][0][:-1]=='Weight_for_Elongation_Sub_':
			list_constraints_weights_e[int(Names[j][0][-1])]=float(Values[j][0])
			test_sum[int(Names[j][0][-1])]+=list_constraints_weights_e[int(Names[j][0][-1])]
		elif Names[j][0][:-1]=='Weight_for_Percentage_of_Sidelobes_Sub_':
			list_constraints_weights_s[int(Names[j][0][-1])]=float(Values[j][0])
			test_sum[int(Names[j][0][-1])]+=list_constraints_weights_s[int(Names[j][0][-1])]

			
#==================================================================================================================================
# Exceptions : 
#==================================================================================================================================
if len(Values)!=13+Number_of_Subarrays*8:
	raise Exception('Invalid number of inputs in the file '+ File_Parameters+'. Check if there is the inputs for all the subarrays. See the example files in the folder GA_Subarray_Selection/Input_Files_To_Fulfill/Input_Files_Examples')
	
if sum(list_num_pads_subarrays)!=num_pads:
        raise Exception('Invalid numbers of pads per subarray : the sum should be equal to the total number of pads = '+str(num_pads)+'.'+'\n'+'But '+str(sum(list_num_pads_subarrays))+' pads were selected.')
elif len(list_num_pads_subarrays)!=Number_of_Subarrays:
        raise Exception('Invalid numbers of pads per subarray : should have '+str(Number_of_Subarrays)+' numbers of pads')
    
test=sum(len(list_Number_of_Iterations)+len(list_Population_Size)+len(list_Termination_Condition)+len(list_Threshold)+len(list_Mutation_Rate)+len(list_Tournament_Size)+len(list_Number_for_Elitism))
#if test/7!=len(list_Number_of_Iterations):
	#raise Exception('Not good number of inputs in the lists for GA Parameters')
    
for i in range(0,Number_of_Subarrays):    
    if test_sum[i]!=0.0 and test_sum[i]!=1.0:
	if list_constraints_weights_res[i]<0 or list_constraints_weights_mrs[i]<0 or list_constraints_weights_e[i]<0 or list_constraints_weights_s[i]<0:
		raise Exception('Negative constraint weight')
	list_constraints_weights_res[i]=round(list_constraints_weights_res[i]/test_sum[i],2)
    	list_constraints_weights_mrs[i]=round(list_constraints_weights_mrs[i]/test_sum[i],2)
    	list_constraints_weights_e[i]=round(list_constraints_weights_e[i]/test_sum[i],2)
    	list_constraints_weights_s[i]=round(list_constraints_weights_s[i]/test_sum[i],2)
     

list_objective_constraints=[]
list_objective_constraints.append(list_objective_constraints_res)
list_objective_constraints.append(list_objective_constraints_mrs)
list_objective_constraints.append(list_objective_constraints_e)
list_objective_constraints.append(list_objective_constraints_s)

if len(list_objective_constraints_res)!=Number_of_Subarrays:
	raise Exception('Invalid number of objective constraints for Spatial Resolution') 
elif len(list_objective_constraints_mrs)!=Number_of_Subarrays:
	raise Exception('Invalid number of objective constraints for Maximum Recoverable Scale ') 
elif len(list_objective_constraints_e)!=Number_of_Subarrays:
	raise Exception('Invalid number of objective constraints for Elongation')
elif len(list_objective_constraints_s)!=Number_of_Subarrays:
	raise Exception('Invalid number of objective constraints for Sidelobe Percentage')
elif now in list_objective_constraints_res:
	raise Exception('Parameter Weight_for_Spatial_Resolution_Sub_'+str		       (list_objective_constraints_res.index(now))+' not read')
elif now in list_objective_constraints_mrs:
	raise Exception('Parameter Weight_for_Maximum_Recoverable_Scale_Sub_'+str(list_objective_constraints_mrs.index(now))+' not read')
elif now in list_objective_constraints_e:
	raise Exception('Parameter Weight_for_Elongation_Sub_'+str(list_objective_constraints_e.index(now))+' not read')
elif now in list_objective_constraints_s:
	raise Exception('Parameter Weight_for_Percentage_of_Sidelobes_Sub_'+str(list_objective_constraints_s.index(now))+' not read')

list_constraints_weights=[]
list_constraints_weights.append(list_constraints_weights_res)
list_constraints_weights.append(list_constraints_weights_mrs)
list_constraints_weights.append(list_constraints_weights_e)
list_constraints_weights.append(list_constraints_weights_s)

if len(list_constraints_weights_res)!=Number_of_Subarrays:
	raise Exception('Invalid number of objective constraints weights for Spatial Resolution') 
elif len(list_constraints_weights_mrs)!=Number_of_Subarrays:
	raise Exception('Invalid number of objective constraints weights for Maximum Recoverable Scale ') 
elif len(list_constraints_weights_e)!=Number_of_Subarrays:
	raise Exception('Invalid number of objective constraints weights for Elongation')
elif len(list_constraints_weights_s)!=Number_of_Subarrays:
	raise Exception('Invalid number of objective constraints weights for Sidelobe Percentage')
elif now in list_constraints_weights_res:
	raise Exception('Parameter Spatial_Resolution_Sub_'+str		       (list_constraints_weights_res.index(now))+' not read')
elif now in list_constraints_weights_mrs:
	raise Exception('Parameter Maximum_Recoverable_Scale_Sub_'+str(list_constraints_weights_mrs.index(now))+' not read')
elif now in list_constraints_weights_e:
	raise Exception('Parameter Elongation_Sub_'+str(list_constraints_weights_e.index(now))+' not read')
elif now in list_constraints_weights_s:
	raise Exception('Parameter Percentage_of_Sidelobes_Sub_'+str(list_constraints_weights_s.index(now))+' not read')

#==================================================================================================================================
#==================================================================================================================================
# Parameters summary : 
#==================================================================================================================================
#==================================================================================================================================
Score_Storage=[]
Iter_Storage=[]


print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARAMETERS FOR THE TEST !!!!!!!!!!!!!!!!!!!!!"
print list_Number_of_Iterations
print list_Population_Size
print list_Termination_Condition
print list_Threshold
print list_Mutation_Rate
print list_Tournament_Size
print list_Number_for_Elitism
print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

for w in range(0,len(list_Number_of_Iterations)):
	Number_of_Iterations=list_Number_of_Iterations[w]
	Iter_Storage.append(Number_of_Iterations)
	Population_Size=list_Population_Size[w]
	Termination_Condition=list_Termination_Condition[w]
	Threshold=list_Threshold[w]
	Mutation_Rate=list_Mutation_Rate[w]
	Tournament_Size=list_Tournament_Size[w]
	Number_for_Elitism=list_Number_for_Elitism[w]

	if not(Mutation_Rate>=0 and Mutation_Rate<=1):
  	  raise Exception('Invalid mutation rate '+str(w)+' : should be a float in [0,1]')
	elif (Tournament_Size>Population_Size or Tournament_Size<0):
	    raise Exception('Invalid tournament size '+str(w)+' : should be an integer in [0,'+str(Population_Size)+'] (according to the Population size)')
	elif not(Number_for_Elitism>=0 and Number_for_Elitism<=Population_Size):
 	   raise Exception('Invalid number for elitism '+str(w)+' : should be an integer in [0,'+str(Population_Size)+'] (according to the Population size)')

	print '\n'
	print "________________________________________________________________"
	print "_____________________ START OF GA (test : "+str(w)+")_________________"
	print "________________________________________________________________"
	print '\n'

	print "Version : ",__version__

	if True:

	    print '\n'
	    print "----------------------------------------------------"
	    print 'PARAMETERS'
	    print "Configuration file : ",Configuration_Input_File
	    print "Source declination : ",Source_Declination
	    print "Source hour angle : ",Source_Hour_Angle
	    print "Total number of pads : ",num_pads
	    print "Number of subarrays : ",Number_of_Subarrays
	    print "Number of pads per subarray : ",list_num_pads_subarrays
	    print "Subarrays weights : ",list_subarrays_weights
	    print "Number of iterations : ", Number_of_Iterations
	    print "Population size : ",Population_Size
	    print "Termination condition : ",Termination_Condition
 	    if Termination_Condition==True:
		print "Threshold : ",Threshold
	    print "Mutation rate : ",Mutation_Rate
	    print "Tournament size : ",Tournament_Size
	    print "Number for elitism : ",Number_for_Elitism
	    print "Objective constraints : ",list_objective_constraints
 	    print "Constraints weights : ",list_constraints_weights
	    print "----------------------------------------------------"


#==================================================================================================================================
#==================================================================================================================================
# Test : 
# - initialize a population thanks to a configuration file
# - evolve a population doing a genetic algorithm
#==================================================================================================================================
#==================================================================================================================================

	print '\n'
	print "----------------------------------------------------"
	print 'SCORE EVOLUTION '         
	pop=uf.Init_Pop(Configuration_Input_File,configuration_Manager,Source_Declination,Source_Hour_Angle, Number_of_Subarrays, 	list_num_pads_subarrays,Population_Size,list_objective_constraints,list_constraints_weights,list_subarrays_weights)
	t0=time.clock()
	t2=time.time()
	new_pop,Scores, Indexes=uf.Evolve_Pop(pop,Number_of_Iterations,Termination_Condition,Threshold,configuration_Manager,Mutation_Rate,Tournament_Size,Number_for_Elitism)
	t1=time.clock()
	t3=time.time()
	t=t1-t0
	t4=t3-t2
	print "----------------------------------------------------"

#==================================================================================================================================
#==================================================================================================================================
# Displays :
#==================================================================================================================================
#==================================================================================================================================


#==================================================================================================================================
# Display of the initial configuration : 
#==================================================================================================================================
	if Display_Screen_Results==True:

	   print '\n'
 	   print "----------------------------------------------------"
 	   print 'INITIAL RESULTS ' 
 	   best_init,pos=pop.get_Fittest()
 	   print "Best score = ", best_init.get_Score()
 	   print "Best array = "
 	   best_init.display_Array()
 	   print "----------------------------------------------------"
 	   if True:
	        print '\n'
	        print "----------------------------------------------------"
	        print 'INITIAL CONSTRAINTS RESULTS '
	        Constraints=[]
	        for i in range(0,best_init.configuration_Manager.get_Number_Subarrays()):
	            baselines = ccs.calc_baselines(best_init,i)
	            res,mrs,elongation,sidelobe_percentage=ccs.calc_constraints(-1., -50., baselines)
	            Constraints.append([])
	            Constraints[i].append(res)
	            Constraints[i].append(mrs)
	            Constraints[i].append(elongation)
	            Constraints[i].append(sidelobe_percentage)
	            print '------ res, mrs, elongation, sidelobe_percentage for subarray ',i,' ------'
	            print Constraints[i]
	        print "----------------------------------------------------"

    

#==================================================================================================================================
# Display of final result : 
#==================================================================================================================================

 	   print '\n'
	   print "----------------------------------------------------"
	   print 'FINAL RESULTS'
	   print 'Result with ', Number_of_Iterations, 'generations in ',t,' CPU time(s) and ',t4,' real time(s) : '  
	best,pos=new_pop.get_Fittest()
	Score_Storage.append(best.get_Score())
	if Display_Screen_Results==True:

	    print "Best score = ", best.get_Score()
	    print "Best array = "
	    best.display_Array()
	    print "----------------------------------------------------"
	    print '\n'
	    print "----------------------------------------------------"
	    print 'FINAL CONSTRAINTS RESULTS '
	Constraints=[[],[],[],[]]
	for i in range(0,best.configuration_Manager.get_Number_Subarrays()):
	    baselines = ccs.calc_baselines(best,i)
	    res,mrs,elongation,sidelobe_percentage=ccs.calc_constraints(-1., -50., baselines)
	    Constraints.append([])
	    Constraints[i].append(res)
	    Constraints[i].append(mrs)
	    Constraints[i].append(elongation)
	    Constraints[i].append(sidelobe_percentage)

	    if Display_Screen_Results==True:

	      print '------ Spatial Resolution, Maximum Recoverable Scale, Elongation, Sidelobe Percentage for Subarray ',i,' ------'
	      print Constraints[i]
	      print "----------------------------------------------------"


#==================================================================================================================================
#==================================================================================================================================
# Saving of the test : 
#==================================================================================================================================
#==================================================================================================================================


	if Save_Folder_Results==True:

	    print '\n'
	    print "----------------------------------------------------"
	    print 'SAVING OF THE RESULTS'
	    stockage_file=uf.Subarrays_to_cfg(best,Configuration_Input_File,Save_Folder_Results,Constraints,best.get_Score(),t)
	    print 'Subarrays files at CASA format in folder : GA_Subarray_Selection/Results/'+stockage_file[0]+"/Subarrays_Storage"
	    print 'Final constraints results in folder : GA_Subarray_Selection/Results/'+stockage_file[0]
	    copy_file=os.path.basename(File_Parameters)
	    shutil.copy(File_Parameters,'GA_Subarray_Selection/Results/'+stockage_file[0]+'/'+copy_file)
	    print 'Parameters file in folder : GA_Subarray_Selection/Results/'+stockage_file[0]
	    print "----------------------------------------------------"


#==================================================================================================================================
# Display of the evolution of the score : 
#==================================================================================================================================
	if Display_Figure_Score_Evolution==True:
 	    plt.figure(1)
	    plt.plot(Indexes,Scores)
	    plt.title('Evolution of the score')
	    if Save_Folder_Results==True:
	       plt.savefig('GA_Subarray_Selection/Results/'+stockage_file[0]+'/Score_Evolution.png')
	       print 'Figure of score evolution in folder : GA_Subarray_Selection/Results/'+stockage_file[0]
	    plt.show()
	    print "----------------------------------------------------"
        
	print '\n'
	print "________________________________________________________________"
	print "_____________________ END OF GA (test : "+str(w)+")_____________________"
	print "________________________________________________________________"
	print '\n'

#==================================================================================================================================
#==================================================================================================================================
# CASA simulation on the subarrays from the array solution : 
#==================================================================================================================================
#==================================================================================================================================


	if Launch_CASA_Simulation==True and Save_Folder_Results==True:

	    print '\n'
	    print "________________________________________________________________"
	    print "_____________________ START OF CASA SIMULATION _________________"
	    print "________________________________________________________________"
	    print '\n'


	    for subarray_file in os.listdir('GA_Subarray_Selection/Results/'+stockage_file[0]+"/Subarrays_Storage"):
		print "............ CASA simulation of "+subarray_file+"............."
		time.sleep(3)
		os.system("python GA_Subarray_Selection/CASA_Sim/arrayConfiguration.py -i GA_Subarray_Selection/Results/"+stockage_file[0]+"/Subarrays_Storage/"+subarray_file+" -t casa -d "+str(Source_Declination))
		time.sleep(5)


	    print '\n'
	    print "________________________________________________________________"
	    print "_____________________ END OF CASA SIMULATION ___________________"
	    print "________________________________________________________________"
	    print '\n'


data = {}

data['Score'] = Score_Storage
data['Population'] = list_Population_Size
data['Mutation'] = list_Mutation_Rate
data['Tournament'] = list_Tournament_Size
data['Elitism'] = list_Number_for_Elitism
data['Iteration'] = list_Number_of_Iterations 

pkl_file = open('GA_Subarray_Selection/Results/score_population.pkl','wb')
pickle.dump(data,pkl_file)
pkl_file.close()

x = list_Population_Size
y = Score_Storage


t=time.strftime("%d-%m-%Y_%H:%M:%S")
plt.figure()
plt.plot(x,y)
plt.xlabel('Number for Population',fontsize=15)
plt.ylabel('Fitness value',fontsize=15)
plt.savefig('GA_Subarray_Selection/Results/Dispersion_Score_Population_'+str(t)+'.png')

