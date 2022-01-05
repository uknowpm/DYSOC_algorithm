#----------------*****  Contact Information *****--------------------------
#   Primary Contact (Implementation Questions, Bug Reports, etc.):
#       Min:uknowpm@gmail.com     Haoyu Jia: leonjiahaoyu@gmail.com
#       Juliane Mueller: juliane.mueller2901@gmail.com
#   Secondary Contact:
#       Christine A. Shoemaker: cas12@cornell.edu
#----------------********************************--------------------------
from utility import *
import copy
import numpy as np
from SLHDstandard import SLHDstandard
from DYCORS_opt import DYCORS_opt
# import cPickle as p
import pickle as p
import os
import fileinput
import sys
import subprocess
import time

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def DYCORSrestart(data, maxeval, NumberNewSamples,data_file,FinalRun,start_program,no_trail):
    m = 2 * (data.dim + 1) # number of points in initial experimental design
    # initialize arrays for collecting results of current trial
    # number of restarts (when algorithm restarts within one trial after
    # encountering local optimum)
    numstart = 0 # collect all objective function values of the current trial here
    Y_all = None # collect all sample points of the current trial here
    S_all = None # best objective function value found so far in the current trial 
    value = np.inf # best objective function value found so far in the current trial 
    numevals = 0 # number of function evaluations done so far
    Fevaltime_all = None # collect all objective function evaluation times of the current trial here
    incre=1
    while numevals < maxeval: # do until max. number of allowed f-evals reached
        numstart = numstart + 1 # increment number of algorithm restarts

        #if numstart ==1:
        if numstart == 0:
           data.S=np.loadtxt(os.getcwd()+'/'+str(no_trail)+'/S_initial.txt') #os.getcwd()+'\\'+str(no_trail)+'\
        
        #print data.S
        # create initial experimental design by symmetric Latin hypercube sampling
        # for cubic and thin-plate spline RBF: rank_P must be Data.dim+1
        else:
            rank_P = 0
            # regenerate initial experimental design until matrix rank is dimension+1
            while rank_P != data.dim + 1:
                data.S = SLHDstandard(data.dim, m)
            # matrix augmented with vector of ones for computing RBF model parameters
                P = np.concatenate((np.ones((m, 1)), data.S), axis = 1)
                rank_P = np.linalg.matrix_rank(P)
            #save initial data.S value in file, in case being stopped when evaluating the initial values
                np.savetxt(os.getcwd()+'/'+str(no_trail)+'/S_initial.txt', data.S, fmt='%.18e', delimiter=' ', newline='\n')
        
        # for the current number of starts, run local optimization
        data,localminflag,flag_resubmit= DYCORS_opt(data, maxeval - numevals, NumberNewSamples,data_file,numevals,numstart,FinalRun,incre,no_trail)
        
        if numevals+data.NumberFevals+NumberNewSamples>maxeval:
            if localminflag==0:
                maxeval=1 # so that we can get out of the loop

        # update best solution found if current solution is better than best
        # point found so far
        if data.Fbest < value:
            solution = data.xbest # best point
            value = data.Fbest # best function value
        solution=np.array(solution)
        solution=solution.reshape(1,data.dim)
        if type(value) is not int:
            value=np.array(value)
            value= value.reshape(1)
        else:
            value=np.array([value])
        sol_v=np.concatenate((solution,[value]),axis=1)
        np.savetxt(os.getcwd()+'/'+str(no_trail)+'/solution_value.txt',sol_v)
        
        if type(Fevaltime_all) == type(None):
            Fevaltime_all = data.fevaltime
            Y_all = data.Y
            S_all = data.S
        else:
            Fevaltime_all = np.concatenate((Fevaltime_all, data.fevaltime), axis = 0)
            Y_all = np.concatenate((Y_all, data.Y), axis = 0)
            S_all = np.concatenate((S_all, data.S), axis = 0)
        
        
        numevals = numevals + data.NumberFevals
        
        
        #check number of evaluations done within this submitted job

    data.S = S_all
    data.Y = Y_all
    data.fevaltime = Fevaltime_all
    data.xbest = solution
    data.Fbest = value
    data.NumberFevals = numevals
    data.NumberOfRestarts = numstart
    
    
    return data #user_time,user_time_initial
