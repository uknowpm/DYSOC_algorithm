#----------------*****  Contact Information *****--------------------------
#   Primary Contact (Implementation Questions, Bug Reports, etc.):
#       Haoyu Jia: leonjiahaoyu@gmail.com
#       Juliane Mueller: juliane.mueller2901@gmail.com
#   Secondary Contact:
#       Christine A. Shoemaker: cas12@cornell.edu
#----------------********************************--------------------------
from utility import *
import copy,os
import numpy as np
from DYCORSrestart import DYCORSrestart

def DYCORSrestartManager(data, maxeval, Ntrials, NumberNewSamples,data_file,NeedRe,FinalRun,st_pt1, st_pt2,start_program):
    solution = Solution()
    solution.BestPoints = np.zeros((Ntrials, data.dim))
    solution.BestValues = np.zeros((Ntrials, 1))
    solution.NumFuncEval = np.zeros((Ntrials, 1))
    solution.AvgFuncEvalTime = np.zeros((Ntrials, 1))
    solution.FuncVal = np.asmatrix(np.zeros((maxeval, Ntrials)))
    solution.DMatrix = np.zeros((maxeval, data.dim, Ntrials))
    solution.NumberOfRestarts = np.zeros((Ntrials, 1))

    a = np.asmatrix(np.zeros((maxeval, Ntrials)))
    for i in range(Ntrials):
        path=os.getcwd()
#         print('path in DycorseManger ',path)
        
        if not os.path.exists(path+'/'+str(i)):
            os.makedirs(path+'/'+str(i))
    for j in range(Ntrials):
        
        # np.random.seed(j + 1)

        # Call the surrogate optimization function
        # Python pass parameter by reference, so we must copy the object
        data_temp = copy.copy(data)
        if NeedRe==0:
            data_temp = DYCORSrestart(data_temp, maxeval, NumberNewSamples,data_file,FinalRun,start_program,j)
        # elif NeedRe==1:
        #     data_temp = DYCORSrestart_re(data_temp, maxeval, NumberNewSamples,data_file,st_pt1,st_pt2,FinalRun,start_program,j)
        
        else:
            raise Exception( "NeedRe can only be 0 or 1")
        
        if FinalRun==2: # never record solution, use result.txt file  
            a[:, j] = np.copy(data_temp.Y)
            
            #data.numeval=0#added by min ignore Ntrial case here
            
            # Gather results in "solution" struct-variable
            solution.BestValues[j] = data_temp.Fbest
            solution.BestPoints[j, :] = data_temp.xbest
            solution.NumFuncEval[j] = data_temp.NumberFevals
            solution.AvgFuncEvalTime[j] = np.mean(data_temp.fevaltime)
            solution.FuncVal[:, j] = np.copy(data_temp.Y)
            solution.DMatrix[:, :, j] = np.copy(data_temp.S)
            solution.NumberOfRestarts[j] = data_temp.NumberOfRestarts
    return solution
