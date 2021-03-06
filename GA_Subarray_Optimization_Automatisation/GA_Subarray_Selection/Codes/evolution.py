#! /usr/bin/env python2
# -*- coding: utf-8 -*-
__author__="Roxane Lassis"
__version__="0.1"
__email__="lassis@etud.insa-toulouse.fr"
'''
    File name: evolution.py
    Author: Roxane Lassis 
    Description : Class GA to evolve a population - Subarray selection using Genetic Algorithm
    Context : ALMA internship
    Date created: 06/2017
    Date last modified: 08/2017
    Python Version: 2.7
'''

import numpy as np
import random
import pandas as pd

import classes as c


#==================================================================================================================================
# Class of the genetic algorithm with operators :   
#-> elitism 
#-> selection_parent 
#-> crossover 
#-> mutation 
#==================================================================================================================================
     

class GA:
    def __init__(self, configuration_Manager,mutation_Rate,tournament_Size,elitism_Num):
        """
       Initialize the genetic algorithm with parameters :
       - configuration_Manager : a class wich allows to configurate the population
       - mutation_Rate : a rate for the mutation
       - tournament_Size : a number of arrays that can participate to the tournament wich select the parents for crossover
       - elitism_Num : a number of survivors to keep for the next population
        """
            
        self.configuration_Manager=configuration_Manager
        self.mutation_Rate=mutation_Rate
        self.tournament_Size=tournament_Size
        self.elitism_Num=elitism_Num
        
#==================================================================================================================================

    def evolve_Population(self,pop):
        """
       Make a new population from the actual population doing :
       - elitism : choice of arrays to keep
       - selection_parent : choice of parents to cross in order to get knew arrays
       - crossover : generate new arrays
       - mutation : introduce mutation in the new population
        """
        new_pop=c.Population(self.configuration_Manager,pop.population_Size(),False)
        self.elitism(pop,new_pop)
        
        for i in range(self.elitism_Num,new_pop.population_Size()):
                parent_1=self.selection_Parent(pop)
                parent_2=self.selection_Parent(pop)
                child=self.crossover(parent_1,parent_2)
                new_pop.save_Array(i, child)
        for i in range(self.elitism_Num,new_pop.population_Size()):
                 self.mutation(new_pop.get_Array(i))
            
        return new_pop
        
#==================================================================================================================================
       
    def elitism(self,pop,new_pop):
            """
            Fitness Based Selection of Survivors :
            save K=elitism_Num arrays (those with the best score)
            in the new population (at the begin of the list of arrays)
            from the last generation
            and return their indexes
            """
            pop_int=pop.copy_Population()
            for i in range(0,self.elitism_Num):
                best_array,pos=pop_int.get_Fittest()
                new_pop.save_Array(i,best_array)
                pop_int.del_Array(pos)

#==================================================================================================================================
    def crossover(self,parent_1,parent_2):
        """
        One point crossover :
        a random crossover point is selected 
        and the tails of the two parents are swapped to get a child
        """
        child=c.Array(self.configuration_Manager)
        pos=random.randint(0,parent_1.size_Array()-1)
        counter=0
        pads_from_p1=[] 
        #copy of the pads of parent_1 situated before pos
        for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
            for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):   
                if counter<pos:
                    child.set_Pad(i,j,parent_1.get_Pad(i,j))
                    pads_from_p1.append(parent_1.get_Pad(i,j))
                counter+=1  
        #list of the pads of parents_2 not present in child
        pads_from_p2=[]       
        for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
            for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):
                if child.contain_pad(parent_2.get_Pad(i,j))==False:
                    pads_from_p2.append(parent_2.get_Pad(i,j))
        #copy those pads from parents_2 in child
        k=0
        for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
            for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):   
                if child.get_Pad(i,j)==None:
                    child.set_Pad(i,j,pads_from_p2[k])
                    k+=1   
        return child
                 
#==================================================================================================================================        

    def mutation(self,array):
            """
            Swap mutation :
            select 2 random pads into 2 different subarrays 
            and interchange them 
            with rate mutation_Rate
            """
            for i in range(0,self.configuration_Manager.get_Number_Subarrays()):
                for j in range(0,self.configuration_Manager.get_Number_Pads_Per_Subarray()[i]):
                    pos_pad_1=[i,j]
                    random_subarray=random.randint(0,self.configuration_Manager.get_Number_Subarrays()-1)
                    if random.random()<self.mutation_Rate and random_subarray!=i:
                        number_pads=self.configuration_Manager.get_Number_Pads_Per_Subarray()[random_subarray]
                        random_pad=random.randint(0,number_pads-1)
                        pos_pad_2=[random_subarray,random_pad]
                        pad_1=array.get_Pad(pos_pad_1[0],pos_pad_1[1])
                        pad_2=array.get_Pad(pos_pad_2[0],pos_pad_2[1])
                        
                        array.set_Pad(pos_pad_1[0],pos_pad_1[1],pad_2)
                        array.set_Pad(pos_pad_2[0],pos_pad_2[1],pad_1)
                        
#==================================================================================================================================                        
            
    def selection_Parent(self,pop):
            """
            Tournament selection :
            select a parent between K=tournament_Size arrays 
            and choose the best one
            """
            tournament=c.Population(self.configuration_Manager,self.tournament_Size,False)
            for i in range(0,self.tournament_Size):
                random_array_index=random.randint(0,pop.population_Size()-1)
                array=pop.get_Array(random_array_index)
                tournament.save_Array(i,array)
            fittest,pos=tournament.get_Fittest()
            return fittest
         
        
