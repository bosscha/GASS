#! /usr/bin/env python2
# -*- coding: utf-8 -*-
__author__="Roxane Lassis"
__version__="0.2"
__email__="lassis@etud.insa-toulouse.fr"
'''
    File name: Useful_Functions_for_GA_Main.py
    Author: Roxane Lassis 
    Description : Useful functions to make optimization - Subarray selection using Genetic Algorithm
                  with parameters given in the file '.txt'
    Context : ALMA internship
    Date created: 06/2017
    Date last modified: 09/2017
    Python Version: 2.7
'''

import numpy as np
from pylab import *
import pandas as pd
import csv
import time
import sys
import os
DOSSIER_COURRANT = os.path.dirname(os.path.abspath(__file__))
DOSSIER_PARENT = os.path.dirname(DOSSIER_COURRANT)
sys.path.append(os.path.join(DOSSIER_PARENT,"Codes"))
import classes as c
import evolution as ev
import compute_constraints_subarray as ccs


#==================================================================================================================================
#==================================================================================================================================
# Useful functions in order to :
# - initialize a population thanks to a configuration file
# - evolve a population doing a genetic algorithm
# - save the subarrays solution to "casa" format and the final cosntaints results
#==================================================================================================================================
#==================================================================================================================================


def Init_Pop(cfg_file,configuration_Manager,source_declination,source_hour_angle, num_subarrays, list_num_pads_subarrays, num_arrays,list_objective_constraints,list_constraints_weights,list_subarrays_weights):
    """
    Initialize a population with :
    - cfg_file : a file containing the pads to consider
    - configuration_Manager : a class wich allows to configurate the population
    - num_subarrays : a number of subarrays in each array
    - list_num_pads_subarrays : a list of the number of pads per subarrays
    - num_arrays : a number of arrays for the population
    """
    cfg = pd.read_csv(
        cfg_file, comment='#', names=['E', 'N', 'U', 'diam', 'pad'], sep='\s+')

    cm=configuration_Manager
    
    xx, yy , zz =cfg[['E', 'N', 'U']].values.transpose()
    diam=cfg[['diam']].values.transpose()
    name =cfg[['pad']].values.transpose()
    s=0
    for i in range(0,num_subarrays):
        s+=list_num_pads_subarrays[i] 
    
    cm.clear()
    for i in range(0,len(xx)):
        pad=c.Pad(xx[i],yy[i],zz[i],diam[0][i],name[0][i])
        cm.add_Pad(pad)

    cm.set_Source_Declination(source_declination)
    cm.set_Source_Hour_Angle(source_hour_angle)

    cm.set_Number_Subarrays(num_subarrays)
    cm.set_Number_Pads_Per_Subarray(list_num_pads_subarrays)
      
    cm.set_Objective_Constraints(list_objective_constraints)
    cm.set_Constraints_Weights(list_constraints_weights)
    cm.set_Subarrays_Weights(list_subarrays_weights)
    
    pop=c.Population(cm,num_arrays,True)
    
    return pop
    
#==================================================================================================================================    

def Evolve_Pop(pop,num_generations,termination_condition,threshold,configuration_Manager,mutation_Rate,tournament_Size,elitism_Num,display=True):
    """
    Evolve a population for K=num_generations generations with parameters :
    - cfg_file : a file containing the pads to consider
    - num_generations : a number of generations
    - configuration_Manager : a class wich allows to configurate the population
    - mutation_Rate : a rate for the mutation
    - tournament_Size : a number of arrays that can participate to the tournament wich select the parents for crossover
    - elitism_Num : a number of survivors to keep for the next population
    - display : boolean to display or not the result every 10 generations
    """
    
    ga=ev.GA(configuration_Manager,mutation_Rate,tournament_Size,elitism_Num)
    counter=0
    Scores=[]
    Indexes=[]
    condition=False
    while counter<=num_generations and not condition:
        pop=ga.evolve_Population(pop)
	best,pos=pop.get_Fittest()
	Scores.append([best.get_Score()])
        Indexes.append(counter)
        if counter%10==0 and display==True : 
            print "Best score _ Generation n°", counter," = ", best.get_Score()
	if termination_condition==True:
		condition=Scores[counter][0]>threshold
		if condition==True:
			print "Threshold reached _  Best score _ Generation n°",counter," = ",Scores[counter][0]	

        counter+=1
    return pop, Scores, Indexes
    
#==================================================================================================================================    
    
def Subarrays_to_cfg(array,cfg_file,saving_constraints_results,Constraints,score,cpu):
    """
    Save the subarrays of the array solution 
    and the final constraints results of each subarray of the array solution
    in the folder "GA_Subarray_Selection/Results/"
    """
    ##Saving of the subarrays at CASA format
    Table=[]
    for i in range(0,array.configuration_Manager.get_Number_Subarrays()):
        Table.append([])
        for j in range(0,array.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):
            Table[i].append([])
            pad=array.get_Pad(i,j)
            Table[i][j].append(pad.get_E())
            Table[i][j].append(pad.get_N())
            Table[i][j].append(pad.get_U())
            Table[i][j].append(pad.get_Diam())
            Table[i][j].append(pad.get_Name())

    stockage_fichier=[]
    t=time.strftime("%d-%m-%Y_%H:%M:%S") #int(time.strftime("%Y/%m/%d/%H:%M:%S"))
    os.mkdir("GA_Subarray_Selection/Results/Results_"+str(t))	
    os.mkdir("GA_Subarray_Selection/Results/Results_"+str(t)+"/Subarrays_Storage")
    stockage_fichier.append("Results_"+str(t))
    for i in range(0,array.configuration_Manager.get_Number_Subarrays()):
        file_0="Subarray_"+str(i)+".cfg"
        with open("GA_Subarray_Selection/Results/Results_"+str(t)+"/Subarrays_Storage/"+file_0,"wb") as fw:
            fw.write("# observatory=ALMA \n")
            fw.write("# coordsys=LOC (local tangent plane) \n")
            fw.write("# x y z diam pad# \n")
    
            co = csv.writer(fw, delimiter=' ', quotechar='"') #, quoting=csv.QUOTE_NONNUMERIC)
            for ligne in Table[i]:
                co.writerow(ligne)
                
    ##Saving of the final constraints results of the array solution of the GA  
    file_1="Final_Constraints_Results.txt"
    with open("GA_Subarray_Selection/Results/Results_"+str(t)+"/"+file_1,"wb") as fw: 
        fw.write('------------------------------------------------------------\n') 
        fw.write("Final Constraints Results for the Array solution\n")
        fw.write("Version : "+str(__version__)+"\n")
        fw.write("Best score : "+str(score)+"\n")
	fw.write("CPU Time : "+str(cpu)+" sec"+"\n")
        fw.write('------------------------------------------------------------\n')
        fw.write("\n")   
        for i in range(0,array.configuration_Manager.get_Number_Subarrays()):
            fw.write('------ Spatial Resolution, Maximum Recoverable Scale, Elongation, Sidelobe Percentage for Subarray '+str(i)+' ------ \n')
            co = csv.writer(fw, delimiter=',', quotechar='"')
            co.writerow(Constraints[i])
            fw.write("\n")  
                
                
    return stockage_fichier  
