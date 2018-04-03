#! /usr/bin/env python2
# -*- coding: utf-8 -*-
__author__="Roxane Lassis"
__version__="0.1"
__email__="lassis@etud.insa-toulouse.fr"
'''
    File name: Codes.py
    Author: Roxane Lassis 
    Description : Classes to define the elements of a population - Subarray selection using Genetic Algorithm
    Context : ALMA internship
    Date created: 06/2017
    Date last modified: 08/2017
    Python Version: 2.7
'''

import numpy as np
import random
import pandas as pd

import compute_constraints_subarray as ccs

#==================================================================================================================================
# Configuration_Manager   
#-> list_pads : list of all the pads to consider
#-> configuration : list of 2 lists : 1) number of subarrays ; 2) number of pads per subarray
#-> objective_constraints : list of objectives for the constraints per subarray : 1) resolution ; 2) maximum recoverable scale ; 3) elongation ; 4) percetage of sidelobes
#-> constraints_weights : list of constraints weights 
#-> subarrays_weights : list of subarrays weights 
#==================================================================================================================================
     
class Configuration_Manager:  
    list_pads=[]
    configuration=[[],[],[]]
    objective_constraints=[[],[],[],[]]
    constraints_weights=[[],[],[],[]]
    subarrays_weights=[]
    
    def add_Pad(self,pad):
        self.list_pads.append(pad)
        
    def get_List_Pads(self):
        return self.list_pads
        
    def get_Pad(self,index):
        return self.list_pads[index]  
    
    def get_Number_Pads(self):
        return len(self.list_pads)
    
    def set_Number_Subarrays(self,num):
        self.configuration[0]=num
        
    def get_Number_Subarrays(self):
        return self.configuration[0]
        
    def set_Number_Pads_Per_Subarray(self,list_num):
        self.configuration[1][:]=[]
        for i in range(0,self.get_Number_Subarrays()):
            self.configuration[1].append(list_num[i])
    
    def set_Source_Declination(self,dec):
        self.configuration.append(dec)
        
    def get_Source_Declination(self):
        return self.configuration[2] 

    def set_Source_Hour_Angle(self,ha):
        self.configuration.append(ha)
        
    def get_Source_Hour_Angle(self):
        return self.configuration[3]      
        
    def get_Number_Pads_Per_Subarray(self):
        list_num=[]
        for i in range(0,len(self.configuration[1])):
            list_num.append(self.configuration[1][i])
        return list_num
                    
    def set_Objective_Constraints(self,list_objective_constraints):
        self.objective_constraints=list_objective_constraints
                    
    def get_Objective_Constraints(self):
        return self.objective_constraints
                
    def set_Constraints_Weights(self,list_constraints_weights):
        self.constraints_weights=list_constraints_weights[:]
        
    def get_Constraints_Weights(self):
        return self.constraints_weights
        
    def set_Subarrays_Weights(self,list_subarrays_weights):
        self.subarrays_weights=list_subarrays_weights
        
    def get_Subarrays_Weights(self):
        return self.subarrays_weights
        
    def clear(self):
        self.list_pads=[]
        self.configuration=[[],[]]
        self.objective_constraints=[[],[],[],[]]
        self.constraints_weights=[]
        self.subarrays_weights=[]
        

#==================================================================================================================================
# Class of pads   
#-> E : position of the pad on east axis
#-> N : position of the pad on north axis
#-> U : position of the pad on vertical axis
#-> diam : diameter of the antenna on this pad
#-> name 
#==================================================================================================================================

class Pad:
    def __init__(self,E, N, U,diam, name):
        self.E=E
        self.N=N
        self.U=U
        self.diam=diam
        self.name=name
    
    def get_E(self):
        return self.E
        
    def get_N(self):
        return self.N
        
    def get_U(self):
        return self.U
        
    def get_Diam(self):
        return self.diam
        
    def get_Name(self):
        return self.name
        

#==================================================================================================================================
# Class of arrays  
#-> configuration_Manager
#-> array : a list of subarrays (wich are a list of pads)
#-> score : score according to how it fits the constraints
#==================================================================================================================================
                         
class Array:
    def __init__(self, configuration_Manager, array=None):
        """
        Initialize an array according to the instructions
        (number of subarrays and of pads per subarrays) 
        given to the configuration_Manager
        """
        self.configuration_Manager=configuration_Manager
        self.array=[]
        self.score=0.0 
        if array is not None:
            self.array=array
        else:        
            #creation of the subarrays of Ni pads
            for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
                self.array.append([])
                #creation of the Ni pads in each subarray
                for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):
                    self.array[i].append(None)  
   
        
    def generate_Array(self):
        """
        Generate an array randomly from the list of pads 
        given to the configuration_Manager
        """
        """
        k=index of the pads in the list of pads 
        i=index of the subarrays
        j=index of the pads in each subarray
        """
        random_list_pads=self.configuration_Manager.get_List_Pads()
        random.shuffle(random_list_pads)
        
        k=0
        for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
            for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):
                pad=random_list_pads[k]
                self.set_Pad(i,j,pad)
                k+=1
        
    def get_Pad(self,num_subarray,index):
        return self.array[num_subarray][index]
        
    def set_Pad(self,num_subarray,index,value):
        self.array[num_subarray][index]=value  
       
    def get_Array(self):
       return self.array
       
    def get_Subarray(self,index):
       return self.array[index]
    
    def access_Score(self):
        return self.score

    def get_Score(self):
        """
        Compute the score of an array 
        thanks to the score of each subarray weighted according to the subarrays_weights
        GOAL : to have the smaller score 
        (since it corresponds to the difference to the objective_constrains)
        """
        list_subarrays_weights=self.configuration_Manager.get_Subarrays_Weights()
        if self.score==0.0:
            for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
                if list_subarrays_weights[i]!=0:
                    self.score+=self.get_Score_Subarray(i)*list_subarrays_weights[i]
        return -self.score*(self.score!=0)+self.score*(self.score==0)
    
    def get_Score_Subarray(self,num_subarray):
        """
        Compute the score of a subarray 
        according to the difference between the constraints computed (res,mrs,elongation,sidelobe_percentage)
        and the 4 objective_constraints 
        and weighted according to the constraints_weights
        """
        list_objective_constraints=self.configuration_Manager.get_Objective_Constraints()
        list_constraints_weights=self.configuration_Manager.get_Constraints_Weights()
        baselines = ccs.calc_baselines(self,num_subarray)
	dec=self.configuration_Manager.get_Source_Declination()
	ha=self.configuration_Manager.get_Source_Hour_Angle()
        res,mrs,elongation,sidelobe_percentage=ccs.calc_constraints(ha,dec, baselines)
        score_subarray=0.0
        score_subarray=list_constraints_weights[0][num_subarray]*abs(res-list_objective_constraints[0][num_subarray])/list_objective_constraints[0][num_subarray]+list_constraints_weights[1][num_subarray]*(list_objective_constraints[1][num_subarray]-mrs)/list_objective_constraints[1][num_subarray]+list_constraints_weights[2][num_subarray]*(elongation-list_objective_constraints[2][num_subarray])/list_objective_constraints[2][num_subarray]+list_constraints_weights[3][num_subarray]*(sidelobe_percentage-list_objective_constraints[3][num_subarray])/list_objective_constraints[3][num_subarray]
        return score_subarray    
        
    def size_Array(self):
        n=0
        for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
            for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):
                n+=1
        return n
    

    def subarray_contain_pad(self,pad,index):
        return pad in self.get_Subarray(index)
        
    def contain_pad(self,pad):
        res=[]
        for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
            res.append(self.subarray_contain_pad(pad,i))
        return True in res            
        
    def contain_None_Pad(self):
        res=False
        for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
            for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):
                pad_int=self.get_Pad(i,j)
                if pad_int==None:
                    res=True
        return res

    def display_Array(self):
        for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
            print '['
            for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):
                pad=self.get_Pad(i,j)
                print pad.get_Name()
            print ']'

                        
#==================================================================================================================================
# Class of populations   
#-> population : list of arrays
#==================================================================================================================================

class Population:
    

    def __init__(self,configuration_Manager, pop_size, init):
        """
        Initialize a population of arrays with size N=pop_size 
        If init=True, it initializes with generate_Array, ie. N arrays randomly
        constructed with the list of pads given to the construction_Manager
        """
        self.configuration_Manager=configuration_Manager
        self.pop_size=pop_size
        self.population=[]
        for i in range(0,self.pop_size):
            self.population.append(None)

        if init:
            for i in range(0,self.pop_size):
                new_array=Array(configuration_Manager)
                new_array.generate_Array()
                self.save_Array(i,new_array)
                
    def set_Array(self,index, array_value):
        self.population[index]=array_value
        
    def get_Array(self, index):
        return self.population[index]
        
        
    def get_Fittest(self):
        """
        Compute the best array of the population
        return the array and its position in the population
        """
        fittest=self.population[0]
        pos=0
        for i in range(0, self.population_Size()):
            if fittest.get_Score()<=self.get_Array(i).get_Score():
                fittest=self.get_Array(i)
                pos=i
        return fittest,pos
        
    def population_Size(self):
        return len(self.population)
        
    def save_Array(self,index,array):
        self.population[index]=array
    
    def del_Array(self,index):
        del self.population[index]
        
    def copy_Population(self):
        copy_pop=Population(self.configuration_Manager, self.pop_size, False)
        copy_pop.population=self.population[:]
        return copy_pop
