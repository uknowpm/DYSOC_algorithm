#----------------********************************--------------------------
# Copyright (C) 2014 Cornell University
# This file is part of the program DYCORS.py
#
#    DYCORS.py is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DYCORS.py is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with DYCORS.py.  If not, see <http://www.gnu.org/licenses/>.
#----------------********************************--------------------------



#----------------*****  Contact Information *****--------------------------
#   Primary Contact (Implementation Questions, Bug Reports, etc.):
#   Juliane Mueller: juliane.mueller2901@gmail.com
#       
#   Secondary Contact:
#       Christine A. Shoemaker: cas12@cornell.edu
#----------------********************************--------------------------
import sys
import os
import importlib
import numpy as np
#from pylab import *
import pickle as p
#import cPickle as p
import time
from utility import *
from DYCORSrestartManager import DYCORSrestartManager

def DYCORS(data_file, maxeval = None, Ntrials = None, NumberNewSamples = None):
    #os.chdir(os.getcwd()+'/DYCORS_blaine12_8_1')
    start_program=time.time()
    PlotResult = 0
    FinalRun,NeedRe = 1,0 # used for old HP system, deleted
    st_pt1, st_pt2 = 0,0 # used for old HP system, deleted
    try:
        ## Start input check
        data = read_check_data_file(data_file)
        maxeval, Ntrials, PlotResult, NumberNewSamples = \
            check_set_parameters(data, maxeval, Ntrials, PlotResult, NumberNewSamples)
        ## End input check

        ## Optimization
        solution = perform_optimization(data, maxeval, Ntrials, NumberNewSamples,data_file,NeedRe,FinalRun,st_pt1, st_pt2,start_program)
        ## End Optimization


        ## Plot Result
        if PlotResult:
            plot_results(solution, maxeval, Ntrials)
        ## End Plot Result
        return solution

    except myException as error:
        print (error.msg)

def perform_optimization(data, maxeval, Ntrials, NumberNewSamples,data_file,NeedRe,FinalRun,st_pt1, st_pt2,start_program):
    data.Ncand = min(100 * data.dim,5000) #revised by Min
    data.fNcand =  data.Ncand/3.0 # added by Min
    data.phifunction = 'cubic'
    data.polynomial = 'linear'
    solution = DYCORSrestartManager(data, maxeval, Ntrials, NumberNewSamples,data_file,NeedRe,FinalRun,st_pt1, st_pt2,start_program)
    return solution

def plot_results(solution, maxeval, Ntrials):
    Y_cur_best = np.zeros((maxeval, Ntrials))
    for ii in range(Ntrials): # go through all trials
        Y_cur = solution.FuncVal[:, ii] # unction values of current trial (trial ii)
        Y_cur_best[0, ii] = Y_cur[0] # first best function value is first function value computed
        for j in range(1, maxeval):
            if Y_cur[j] < Y_cur_best[j-1, ii]:
                Y_cur_best[j, ii] = Y_cur[j]
            else:
                Y_cur_best[j, ii] = Y_cur_best[j-1, ii]
    # compute means over matrix of current best values (Y_cur_best has dimension 
    # maxeval x Ntrials)
    Ymean = np.mean(Y_cur_best, axis = 1)

    Yplot = np.zeros((maxeval, 1)) # initialize vector for plotting results
    # sort results according to best point found till iteration
    # Seriously, why do we need that????

    X = np.arange(1, maxeval + 1)
    plot(X, Ymean)
    xlabel('Number Of Function Evaluations')
    ylabel('Average Best Objective Function Value In %d Trials' % Ntrials)
    draw()
    #show()
    savefig('DYCORS_Plot')

def read_check_data_file(data_file):
    if not isinstance(data_file, str):
        raise myException('''You have to supply a file name with your data. \
            \n\tSee example files and tutorial for information how to define problems.''')
    try:
        module = importlib.import_module(data_file)
        data = getattr(module, data_file)()
    except ImportError:
        raise myException('''The data file is not found in the current path\
            \n\tPlease place the data file in the path.''')
    except AttributeError:
        raise myException('''The function name must be the same with the data file name.\
            \n\tSee example files and tutorial for information how to define the function.''')

    data.validate()
    return data

def check_set_parameters(data, maxeval, Ntrials, PlotResult, NumberNewSamples):
    if maxeval == None:
        print ('''No maximal number of allowed function evaluations given.\
                \n\tI use default value maxeval = 20 * dimension.''')
        maxeval = 20 * data.dim
    if not isinstance(maxeval, int) or maxeval <= 0:
        raise myException('Maximal number of allowed function evaluations must be positive integer.\n')

    if Ntrials == None:
        print ('''No maximal number of trials given.\
                \n\tI use default value NumberOfTrials=1.''')
        Ntrials = 1
    if not isinstance(Ntrials, int) or Ntrials <= 0:
        raise myException('Maximal number of trials must be positive integer.\n')

    if PlotResult == None:
        print ('''No indication if result plot wanted.\
                \n\tI use default value PlotResult=1.''')
        PlotResult = 1
    elif abs(PlotResult) > 0:
        PlotResult = 1

    if NumberNewSamples == None:
        print ('''No number of desired new sample sites given.\
                \n\tI use default value NumberNewSamples=1.''')
        NumberNewSamples = 1
    if not isinstance(NumberNewSamples, int) or NumberNewSamples < 0:
        raise myException('Number of new sample sites must be positive integer.\n')

    return maxeval, Ntrials, PlotResult, NumberNewSamples


if __name__ == "__main__":
    print( 'This is start for DYCORS')
    solution = DYCORS(sys.argv[1], int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),\
                             int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),int(sys.argv[9]))
    print ('BestValues', solution.BestValues)
    print ('BestPoints', solution.BestPoints)
    print ('NumFuncEval', solution.NumFuncEval)
    print ('AvgFUncEvalTime', solution.AvgFuncEvalTime)
    print ('DMatrix', solution.DMatrix.shape)
    print ('NumberOfRestarts', solution.NumberOfRestarts)
