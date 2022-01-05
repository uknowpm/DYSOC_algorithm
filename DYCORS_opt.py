#----------------*****  Contact Information *****--------------------------
#   Primary Contact (Implementation Questions, Bug Reports, etc.):
#       Haoyu Jia: leonjiahaoyu@gmail.com
#       Juliane Mueller: juliane.mueller2901@gmail.com
#   Secondary Contact:
#       Christine A. Shoemaker: cas12@cornell.edu
#----------------********************************--------------------------
from utility import *
import numpy as np
import scipy.spatial as scp
import time
import math, copy
from InitialRBFMatrices import InitialRBFMatrices
from Minimize_Merit_Function_rbf import Minimize_Merit_Function_rbf
from phi import phi
import os, shutil
import copy
from ComputeRBF import ComputeRBF 

import subprocess
#Min create this
import multiprocessing
from multiprocessing import Pool, Process
import importlib


##### Min created this for blaine ##################
def run_folders_parallel(a):
    data_file = a[0]
    da = a[1]
    num = a[2]
    numevals = a[3]
    p = a[4]
    #num is number of processors using
    current = multiprocessing.current_process() #for each pool, only four continous four numbers could be
    # p='/data/ese-pangm/new_Scenario1_dsep'
    pp=p+'/proc/' # '/' using mac/linux(in windows should be '\\')

#     print('p in DYCORS_opt run_folder',pp)
    
    if current._identity[0]%num==0:
        x=num
    else:
        x=current._identity[0]%num
    print('p in DYCORS_opt run_folder',x)

    mypath=pp+str(x) 
    #change directory in to each folder and run the expensive function
    os.chdir(mypath)
#     print('current_path',mypath)
    module = importlib.import_module(data_file)
    data = getattr(module, data_file)()
    time_start = time.time()
#     print('da',type(da))
    da.append(x)
#     import pdb
#     pdb.set_trace()
    
#     print('da',da)
    ret_value = data.objfunction(da)
    obj_value = ret_value[0]
    #all constraint values
    cst_value = ret_value[1]
    ret_time = time.time() - time_start
    #change directory back in case of errors outside this function
    os.chdir(p)
    iterationnum=numevals+x
    total_Q=ret_value[2]
    area_subsidence5=ret_value[5]
    sum_sub=ret_value[4]
    status=ret_value[6]
    da = da[:-1]
    return obj_value ,  ret_time, da, iterationnum, cst_value,total_Q,sum_sub,area_subsidence5,status

##### End created  ##################


def DYCORS_opt(data, maxeval, NumberNewSamples,data_file,numevals,numrestat,FinalRun,incre,no_trail):
    '''LocalStochRBFstop is the local optimization routine. It iterates at most
    until totally maxeval points have been evaluated, or, if a local minimum
    has been found, the routine terminates in less than maxeval evaluations,
    but will restart from scratch to use up all remaining function evaluations
    points.

    Input: 
    Data: struct-variable with problem information (variable bounds,
           objective/simulation function handle, etc.)
    maxeval: maximum number of function evaluations for each trial
    NumberNewSamples: number of points where objective/simulation function
                       is evaluated in every iteration of the algorithm; if
                       NumberNewSamples > 1, evaluations are in parallel, in
                       which case Matlab Parallel Computing Toolbox is
                       required.

    Output:
    Data: updated struct-variable containing the results of the current run
           until stop at local minimum, or stop after maxeval evaluations
    '''
    #print data.xup
    dictory = os.getcwd()+'/'#'/data/ese-pangm/new_Scenario1_dsep/'
    penalty_coeff=10**10
    incre_resubmit=2000
    xrange = data.xup - data.xlow
    minxrange = np.amin(xrange) # smallest variable range
    m = data.S.shape[0] # number of points already sampled
    pdim=0
    phi0=1
    flag_resubmit=0
    #maxeval = 5
    # scale design points to actual dimensions
    data.S = np.multiply(np.tile(data.xup - data.xlow, (m, 1)), data.S) + np.tile(data.xlow, (m, 1))

    data.m = min(m, maxeval) # in case number of point in initial experimental design exceed max. number of allowed evaluations
    data.fevaltime = np.asmatrix(np.zeros((maxeval, 1))) # initialize vector with time for function evaluations
    data.Y = np.asmatrix(np.zeros((maxeval, 1))) # initialize array with function values
    data.S = data.S[:data.m, :] # in case Data.m>maxeval, throw additional points away
    if maxeval > data.m:
        # initialize array with sample points (first Data.m points are from initial experimental design)
        data.S = np.concatenate((data.S, np.zeros((maxeval - data.m, data.dim))), axis = 0) 
    
    # algorithm parameters
    sigma_stdev_default = 0.2 * minxrange
    sigma_stdev = sigma_stdev_default # current mutation rate 
    maxshrinkparam = 5 # maximal number of shrikage of standard deviation for normal distribution when generating the candidate points
    failtolerance = max(5,data.dim)
    succtolerance =3
    
    # initializations for record in the "flag.txt" file
    iterctr = 0 # number of iterations
    shrinkctr = 0 # number of times sigma_stdev was shrunk
    failctr = 0 # number of consecutive unsuccessful iterations
    localminflag = 0  # indicates whether or not xbest is at a local minimum
    succctr=0 # number of consecutive successful iterations
    
    #Min---------generate "Flag.txt"file to locate where it was before the cut by Yellowstone after 12hours
    flag_locate=dictory+str(no_trail)+'/Flag.txt'
    with open(flag_locate,'a+') as flt:
        flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0))
    
    # add by Min for constraint handling
    #data.Cst_no = 1# number of constraint exists, only one so far
    data.Cst = np.ones((maxeval,data.Cst_no))*np.inf
    #---------------------------------
    # For serial evaluation of points in initial starting desing:
    # --------------------------SERIAL------------------------------------------
    '''
    for ii in range(data.m): # go through all Data.m points
        time1 = time.time() # start timer for recording function evaluation time
        pp=os.getcwd()
        print 'pp_inLocalstochRBFstop',pp
        
        res = data.objfunction(np.array(data.S[ii, :])) #expensive simulation
        data.fevaltime[ii] = time.time() - time1 # record time for expensive evaluation
        data.Y[ii, 0] = res
        if ii == 0: # initialize best point found so far = first evaluated point
            data.xbest = data.S[0, :]
            data.Fbest = data.Y[ii]
        else:
            if data.Y[ii] < data.Fbest:
                data.Fbest = data.Y[ii]
                data.xbest = data.S[ii, :]
                
    #--------------------------END SERIAL----------------------------------------
    # for parallel evaluation of points in initial starting design delete
    # comments in the following and comment out the serial code above
    #--------------------------PARALLEL------------------------------------------
    '''
    m = data.m # m now takes maxeval in to consideration 
    Y = data.Y
    data.Y_final = copy.deepcopy(Y)
    Y_final = data.Y_final
    Cst = np.zeros([Y.shape[0],data.Cst_no])
    S = data.S
    Time = data.fevaltime
    
    #Min create this for parallel
    pdim=0
    nprocs = NumberNewSamples
    # # of process using (or we can use other number defined)
    sort=np.array([nprocs,data.dim+3])
    result_iter_filename = dictory+str(no_trail)+'/result.txt'
    
    result_restart_filename = dictory+str(no_trail)+'/result_restart_'+str(numrestat)+'.txt'
    # created folders of # of processors
    #print 'LocalStop',os.getcwd()
    for i in range(1,nprocs+1):
        pp=dictory # / is in mac need to change when using windows(should be '\\')
        mypath=pp+'proc/'+str(i)
        if os.path.exists(mypath):
            pass
        else:
            if i==1: # for each folder only one copy of content
                src=pp+'proc/'  #src=pp+'PYTHON'# yellowstone version
                shutil.copytree(dictory+'proc/0',dictory+'/proc/1')
            else:
                src=pp+'proc/0'
                shutil.copytree(src,dictory+'proc/'+str(i))
#             print('create_folder_source',src)
#             print('created_folder','/work/ese-pangm/Min/subsidence/new_Scenario1_sep/proc/'+str(i))
    #start_time_count_initial=resource.getrusage(resource.RUSAGE_SELF)
    
    ########start
    time_start_commui = time.time()
    if nprocs<m:
        #eval_initial is recording how many previous evaluation from the restart here is not necessarily starting as there maybe restart 
        eval_initial=0
        pool = Pool(processes=nprocs) 
#         import pdb 
#         pdb.set_trace()
        #print 'S.tolist()[eval_initial:eval_initial+nprocs]',S[eval_initial:eval_initial+nprocs].shape
        while (eval_initial+nprocs)<m: # to use whether it's appropriate to use all processors have
            #numevals used for writing txt file (calculateing iteration number)
#             import pdb
#             pdb.set_trace()
            result = pool.map(run_folders_parallel,((data_file,ij,nprocs,numevals,dictory) for ij in S.tolist()[eval_initial:eval_initial+nprocs]))
            #--Generate result.txt in case 12hr constraint meet in yellowstone ----------------------------
            #sort result in order make #iteration in order
            #timing on opening the pool and getting result in time_end_opti
            time_end_commui = time.time()-time_start_commui
            time_start_commui = time.time()
            with open(result_iter_filename,'a+') as fid:
                for ii in range(len(result)):
                    dummy_array=np.array(result[ii][4])
                    dummy_array[dummy_array<0]=0
                    final_result = result[ii][0]+penalty_coeff*sum(np.square(dummy_array))
                    fid.write('\n'+str('%-5d'%result[ii][3])+' '+' '.join(str('%12e'%num) for num in result[ii][2])+' '+\
                        str('%12e'%result[ii][0])+' '+' '.join(str('%12e'%num) for num in result[ii][4])+' '+str('%12e'%final_result )\
                        +' '+str('%12e'%result[ii][6])+' '+str('%12e'%result[ii][7])+' '+str('%12e'%result[ii][1]))
        
            with open(result_restart_filename,'a+') as fid1:
                for ii in range(len(result)):
                    dummy_array=np.array(result[ii][4])
                    dummy_array[dummy_array<0]=0
                    final_result = result[ii][0]+penalty_coeff*sum(np.square(dummy_array))
                    fid1.write('\n'+str('%-5d'%result[ii][3])+' '+' '.join(str('%12e'%num) for num in result[ii][2])+' '+\
                        str('%12e'%result[ii][0])+' '+' '.join(str('%12e'%num) for num in result[ii][4])+' '+str('%12e'%final_result )\
                        +' '+str('%12e'%result[ii][6])+' '+str('%12e'%result[ii][7])+' '+str('%12e'%result[ii][1]))
            for rst in range(len(result)):
                Y[eval_initial+rst, 0] = result[rst][0]
                Time[eval_initial+rst, 0] = result[rst][1]
                S[eval_initial+rst,:]=result[rst][2]
                Cst[eval_initial+rst,:]=result[rst][4]
                dummy_array=np.array(result[rst][4])
                dummy_array[dummy_array<0]=0
                final_result = result[rst][0]+penalty_coeff*sum(np.square(dummy_array))
                Y_final[eval_initial+rst, 0] = final_result
                    
            eval_initial=eval_initial+nprocs
            numevals=numevals+nprocs # numevals keep track of total number of evals, including previous restarts part
            with open(flag_locate,'a+') as flt:
                    flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% e,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,time_end_commui,'communi'))
            
            if (numevals/incre_resubmit) >= (incre*NumberNewSamples):
                flag_resubmit=1
                incre=incre+1
                with open(flag_locate,'a+') as flt:
                    flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% d,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,incre,'incre'))
                break
        pool.close()
        pool.join()       
        if flag_resubmit==0: #need to resubmit the job    
            #finish the rest evals in initialization part 
            pool_1= Pool(processes=m-eval_initial) # use up all the processors owned (same as nprocs)
            #result = pool.map(run_folders_parallel,((data_file,i,nprocs,numevals) for i in S.tolist()[0:m])) #numevals used for writing txt file (calculateing iteration number)
            result = pool_1.map(run_folders_parallel,((data_file,ijk,m-eval_initial,numevals,dictory) for ijk in S.tolist()[eval_initial:m]))
            pool_1.close()
            pool_1.join()
            time_end_commui= time.time()-time_start_commui
            with open(result_iter_filename,'a+') as fid:
                for ii in range(len(result)):
                    dummy_array=np.array(result[ii][4])
                    dummy_array[dummy_array<0]=0
                    final_result = result[ii][0]+penalty_coeff*sum(np.square(dummy_array))
                    fid.write('\n'+str('%-5d'%result[ii][3])+' '+' '.join(str('%12e'%num) for num in result[ii][2])+' '+\
                            str('%12e'%result[ii][0])+' '+' '.join(str('%12e'%num) for num in result[ii][4])+' '+str('%12e'%final_result)\
                        +' '+str('%12e'%result[ii][6])+' '+str('%12e'%result[ii][7])+' '+str('%12e'%result[ii][1]))
            with open(result_restart_filename,'a+') as fid1:
                for ii in range(len(result)):
                    dummy_array=np.array(result[ii][4])
                    dummy_array[dummy_array<0]=0
                    final_result = result[ii][0]+penalty_coeff*sum(np.square(dummy_array))
                    fid1.write('\n'+str('%-5d'%result[ii][3])+' '+' '.join(str('%12e'%num) for num in result[ii][2])+' '+\
                        str('%12e'%result[ii][0])+' '+' '.join(str('%12e'%num) for num in result[ii][4])+' '+str('%12e'%final_result)\
                        +' '+str('%12e'%result[ii][6])+' '+str('%12e'%result[ii][7])+' '+str('%12e'%result[ii][1]))
            
            for rst in range(len(result)):
                Y[eval_initial+rst, 0] = result[rst][0]
                Time[eval_initial+rst, 0] = result[rst][1]
                S[eval_initial+rst,:]=result[rst][2]
                Cst[eval_initial+rst,:]=result[rst][4]
                dummy_array=np.array(result[rst][4])
                dummy_array[dummy_array<0]=0
                final_result = result[rst][0]+penalty_coeff*sum(np.square(dummy_array))
                Y_final[eval_initial+rst, 0] = final_result
            numevals=numevals+m-eval_initial
            with open(flag_locate,'a+') as flt:
                        flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% e,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,time_end_commui,'communi'))
            
    else:
        pool = Pool(processes=m)   
        result = pool.map(run_folders_parallel,((data_file,iq,m,numevals,dictory) for iq in S.tolist()[0:m]))
        pool.close()
        pool.join()
        time_end_commui = time.time()-time_start_commui
        with open(result_iter_filename,'a+') as fid:
            for ii in range(len(result)):
                dummy_array=np.array(result[ii][4])
                dummy_array[dummy_array<0]=0
                final_result = result[ii][0]+penalty_coeff*sum(np.square(dummy_array))
                fid.write('\n'+str('%-5d'%result[ii][3])+' '+' '.join(str('%12e'%num) for num in result[ii][2])+' '+\
                        str('%12e'%result[ii][0])+' '+' '.join(str('%12e'%num) for num in result[ii][4])+' '+str('%12e'%final_result)\
                        +' '+str('%12e'%result[ii][6])+' '+str('%12e'%result[ii][7])+' '+str('%12e'%result[ii][1]))
        with open(result_restart_filename,'a+') as fid1:
            for ii in range(len(result)):
                dummy_array=np.array(result[ii][4])
                dummy_array[dummy_array<0]=0
                final_result = result[ii][0]+penalty_coeff*sum(np.square(dummy_array))
                fid1.write('\n'+str('%-5d'%result[ii][3])+' '+' '.join(str('%12e'%num) for num in result[ii][2])+' '+\
                    str('%12e'%result[ii][0])+' '+' '.join(str('%12e'%num) for num in result[ii][4])+' '+str('%12e'%final_result)\
                        +' '+str('%12e'%result[ii][6])+' '+str('%12e'%result[ii][7])+' '+str('%12e'%result[ii][1]))
        for rst in range(len(result)):
                Y[rst, 0] = result[rst][0]
                Time[rst, 0] = result[rst][1]
                S[rst,:]=result[rst][2]
                Cst[rst,:]=result[rst][4]
                dummy_array=np.array(result[rst][4])
                dummy_array[dummy_array<0]=0
                final_result = result[rst][0]+penalty_coeff*sum(np.square(dummy_array))
                Y_final[rst, 0] = final_result
        numevals=numevals+m
        with open(flag_locate,'a+') as flt:
                    flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% e,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,time_end_commui,'communi'))
        
        
    ########stop
    
    # time for updating RS
    time_start_optiupdate=time.time()
    
    data.Fbest = np.amin(Y_final[0:m])
    IDfbest = np.argmin(Y_final[0:m])
    data.xbest = S[IDfbest]
    data.Y = Y
    data.Y_final=Y_final
    data.fevaltime = Time
    data.Cst = Cst

    #--------------------------END PARALLEL--------------------------------------

    # determine pairwise distance between points
    PairwiseDistance = scp.distance.cdist(data.S[0:data.m, :], data.S[0:data.m, :], 'euclidean')
    # initial RBF matrices
    PHI, phi0, P, pdim = InitialRBFMatrices(maxeval, data, PairwiseDistance)
    np.savetxt(dictory+str(no_trail)+'/PHI.txt', PHI, fmt='%.10e', delimiter=' ', newline='\n')
    np.savetxt(dictory+str(no_trail)+'/P.txt', P, fmt='%.10e', delimiter=' ', newline='\n')
    #np.savetxt('phi0.txt', phi0, fmt='%.10e', delimiter=' ', newline='\n') phi0 is integer
    
    # tolerance parameters
    data.tolerance = 0.001 * minxrange * np.linalg.norm(np.ones((1, 3)))

    

    # initializations
    iterctr = 0 # number of iterations
    shrinkctr = 0 # number of times sigma_stdev was shrunk
    failctr = 0 # number of consecutive unsuccessful iterations
    localminflag = 0  # indicates whether or not xbest is at a local minimum
    succctr=0 # number of consecutive successful iterations

    p = data.Y[0, 0]
    # do until max number of f-evals reached or local min found
    weightpattern=np.array([0.3,0.5,0.8,0.95])
    
    #number of constraint
   
    
    time_end_optiupdate=time.time()-time_start_optiupdate
    #Min---------generate "Flag.txt"file to locate where it was before the cut by Yellowstone after 12hours
    flag_locate=dictory+str(no_trail)+'/Flag.txt'
    with open(flag_locate,'a+') as flt:
        flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% e,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,time_end_optiupdate,'opti_update_initial'))
    with open(flag_locate,'a+') as flt:
                    flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% d,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,incre,'incre'))
    if flag_resubmit==0:
        # do until max number of f-evals reached or local min found
        while data.m < maxeval and localminflag == 0:
            time_opti_start=time.time()
            mw=iterctr%len(weightpattern)
            w_r=weightpattern[mw]
            iterctr = iterctr + 1 # increment iteration counter
    
            print ('\n Iteration: %d \n' % iterctr)
            print ('\n fEvals: %d \n' % data.m)
            print ('\n data.Fbest %f \n' % data.Fbest)
        
            # number of new samples in an iteration
            NumberNewSamples = min(NumberNewSamples,maxeval - data.m)
            
           
            Ftransform=[]
            medianF=[]
            coeff=[]
            data.llambda=[]
            data.ctail=[]
            
            # Compute RBF parameters
            a_part1 = np.concatenate((PHI[0:data.m, 0:data.m], P[0:data.m, :]), axis = 1)
            a_part2 = np.concatenate((np.transpose(P[0:data.m, :]), np.zeros((pdim, pdim))), axis = 1)
            a = np.concatenate((a_part1, a_part2), axis = 0)
            eta = math.sqrt((1e-16) * np.linalg.norm(a, 1) * np.linalg.norm(a, np.inf))
    
            for no_cst in range(data.Cst_no):
                
                Cst_data=copy.deepcopy(data.Cst[0:data.m,no_cst])
                Cst_data=np.array(np.asmatrix(Cst_data).T)
                
                Ftransform.append(Cst_data)
               
                medianF.append(np.mean(Cst_data)) # used to be median
                Ftransform[no_cst][Ftransform[no_cst] >  medianF[no_cst]] = medianF[no_cst]
                
                #Ftransform[no_cst]=np.log(Ftransform[no_cst]+1) # log transform to avoid large value
                
                # fit the response surface
                coeff=np.linalg.solve((a + eta * np.eye(data.m + pdim)),\
                    np.concatenate((Ftransform[no_cst], np.zeros((pdim, 1))), axis = 0))
    
                # llambda is not a typo, lambda is a python keyword
                data.llambda.append( coeff[0:data.m])
                data.ctail.append(coeff[data.m: data.m + pdim])
                
    #-------------------------------------------------------------------------------------
            Ftransform1 = np.copy(np.asarray(data.Y_final)[0:data.m])
            medianF1 = np.median(np.asarray(data.Y_final)[0:data.m])
            Ftransform1[Ftransform1 > medianF1] = medianF1
            eta = math.sqrt((1e-16) * np.linalg.norm(a, 1) * np.linalg.norm(a, np.inf))
            coeff1 = np.linalg.solve((a + eta * np.eye(data.m + pdim)),\
                    np.concatenate((Ftransform1, np.zeros((pdim, 1))), axis = 0))
    
            # llambda is not a typo, lambda is a python keyword
            data.llambda1 = coeff1[0:data.m]
            data.ctail1 = coeff1[data.m: data.m + pdim]
            # select the next function evaluation point by DYCORS perturbation:
            # Perturbation probability
            # DDSprob1 = 50/38*min(50.0/data.dim,1)*(1-(math.log(data.m-2*(data.dim+1)+1)/math.log(maxeval-2*(data.dim+1))))
            # DDSprob2 = 50/12*min(50.0/data.dim,1)*(1-(math.log(data.m-2*(data.dim+1)+1)/math.log(maxeval-2*(data.dim+1))))
            # DDSprob= np.hstack((np.ones(38)*DDSprob1,np.ones(12)*DDSprob2))
            DDSprob = np.ones(data.dim)*min(50.0/data.dim,1)*(1-(math.log(data.m-2*(data.dim+1)+1)/math.log(maxeval-2*(data.dim+1))))
            #create candidate points
            CandPoint = np.asmatrix(np.kron(np.ones((data.Ncand, 1)), data.xbest))
            xlow=np.ravel(np.asarray(data.xlow))
            xup=np.ravel(np.asarray(data.xup))
            for ii in range(data.Ncand):
                r=np.random.rand(data.dim)
                ar=r<DDSprob
                if not(any(ar)):
                    r = np.random.permutation(data.dim)
                    ar[r[0]]=True
                for jj in range(data.dim):
                    if ar[jj]:
                        CandPoint[ii,jj] = CandPoint[ii,jj] +sigma_stdev*np.random.randn(1)
                        
                        if CandPoint[ii,jj] < xlow[jj]:
                            CandPoint[ii,jj] = xlow[jj]+ (xlow[jj]-CandPoint[ii,jj])
                            if CandPoint[ii,jj] >xup[jj]:
                                CandPoint[ii,jj] = xlow[jj]
                        elif CandPoint[ii,jj] > xup[jj]:
                            CandPoint[ii,jj] = xup[jj]- (CandPoint[ii,jj]-xup[jj])
                            if CandPoint[ii,jj] <xlow[jj]:
                                CandPoint[ii,jj] = xup[jj]
            #Predict constraint values using RBFs MinPang
            Predict_constraint, NormValue = ComputeRBF(CandPoint, data)

            #NormValue is the distance between points,
            #Predict_constraint is the predicted constriant by RBF
            
            #now need to find predicted feasible points as feasible Candidate points for select [CandPoint_f]
            #print 'Predict_constraint info', type(Predict_constraint), Predict_constraint.shape
           # meet_constraint=np.zeros([data.Ncand,data.Cst_no])
            CandPoint=np.array(CandPoint)
            
            Predict_constraint_C =copy.copy(Predict_constraint)
            Predict_constraint_C[Predict_constraint_C>0]=1
            Predict_constraint_C[Predict_constraint_C<=0]=0
            com_constraint=np.sum(Predict_constraint_C,axis=1)
            predi_feaible_pts=CandPoint[com_constraint==0]
            Pc_need_s = Predict_constraint[com_constraint==0]
            Predict_constraint[Predict_constraint<0]=0
            #print 'predi_feaible_pts',predi_feaible_pts
            #print 'com_constraint',com_constraint
            #print 'data.Cst_no',data.Cst_no.shape
            #print data.fNcand
            result_cand_filename = dictory+str(no_trail)+'/candidate_no.txt'
            with open(result_cand_filename,'a+') as flt:
                flt.write('\n%d,% d,% d,% d' %(iterctr,predi_feaible_pts.shape[0],data.fNcand,data.Ncand))
    
            for i in range(data.Cst_no):
                if predi_feaible_pts.shape[0]<data.fNcand:
                    #Predict_constraint_C[Predict_constraint_C>0]=1
                    #com_constraint=np.sum(Predict_constraint,axis=1) #count how many constraint violations exist
                    predi_feaible_pts_p=CandPoint[com_constraint==1+i]
                    Pc_need = Predict_constraint[com_constraint==1+i]
                    predi_feaible_pts = np.vstack((predi_feaible_pts,predi_feaible_pts_p))
                    Pc_need_s = np.vstack((Pc_need_s,Pc_need)) 
                    with open(result_cand_filename,'a+') as flt:
                        flt.write('\n%d,% d,% d,% d' %(iterctr,predi_feaible_pts.shape[0],data.fNcand,data.Ncand))
                else:
                    break
            
            CandPoint_f=predi_feaible_pts 
            Pc_need_s[Pc_need_s<0]=0
            Predict_constraint_after_filter=Pc_need_s
            # build constraint research to add into weighted merit functions
            data.Predict_constraint=penalty_coeff*np.transpose(np.atleast_2d(np.sum(Predict_constraint_after_filter**2,axis=1)))
            xselected, normval = Minimize_Merit_Function_rbf(data, CandPoint_f, NumberNewSamples,w_r)
            #print 'normval',normval
        
            xselected= np.asmatrix( xselected)    
            # time for optimization
            time_opti_end=time.time()-time_opti_start
           
            # Min created this for parallel######################
            # This is case is special as xselected would be equal to number of processors using as nprocs=NewSamplepoints
            time_commu_start=time.time()
            # goto is a flag to consider the case if it's not final run, then we should stop even maxeval doesn't reach this time,
            # but in next step, we'll exceed maxeval
            goto=0
            # eval_initial is to make sure Fselected and xselected only includes the eval in CAND points with no initialization points
            eval_initial=0
            if nprocs == NumberNewSamples:
                pool1=Pool(processes=nprocs)
                result1 = pool1.map(run_folders_parallel,((data_file,i,nprocs,numevals,dictory) for i in xselected.tolist()))
                pool1.close()
                pool1.join()
                
            else:
                if FinalRun==1: #if it's final run, run till maxeval
                    pool1=Pool(processes=NumberNewSamples)
                    result1 = pool1.map(run_folders_parallel,((data_file,i,NumberNewSamples,numevals,dictory) for i in xselected[0:NumberNewSamples,:].tolist()))
                    np.savetxt('xselected.txt', xselected, fmt='%.10e', delimiter=' ', newline='\n')
                    xselected= xselected[0:NumberNewSamples,:]
                    print("This is final run!")
                else:
                    goto=1 #nothing needs to write, need to get out of this file
                    localminflag = 1 # in order to get out of the while loop
                    maxeval=numevals
                    print('maximum evaluation number will be reached, start another iteration with st_pt1=0, st_pt2=2 please!')
            time_commu_end=time.time()-time_commu_start    
           
            if goto==0:
                Fselected = np.zeros((xselected.shape[0], 1))
                Fselected_final = np.zeros((xselected.shape[0], 1))
                Time = np.zeros((xselected.shape[0], 1))
                lens_xselected=len(xselected.tolist())
                Cst = np.zeros((xselected.shape[0], data.Cst_no)) # only one constraint Min 
                
                with open(result_iter_filename,'a+') as fid1:
                    for ii in range(len(result1)):
                        dummy_array=np.array(result1[ii][4])
                        dummy_array[dummy_array<0]=0
                        final_result = result1[ii][0]+penalty_coeff*sum(np.square(dummy_array))
                        fid1.write('\n'+str('%-5d'%result1[ii][3])+' '+' '.join(str('%12e'%num) for num in result1[ii][2])+' '+\
                        str('%12e'%result1[ii][0])+' '+' '.join(str('%12e'%num) for num in result1[ii][4])+' '+str('%12e'%final_result )\
                        +' '+str('%12e'%result1[ii][6])+' '+str('%12e'%result1[ii][7])+' '+str('%12e'%result1[ii][1]))
                           
                with open(result_restart_filename,'a+') as fid1:
                    for ii in range(len(result1)):
                        dummy_array=np.array(result1[ii][4])
                        dummy_array[dummy_array<0]=0
                        final_result = result1[ii][0]+penalty_coeff*sum(np.square(dummy_array))
                        fid1.write('\n'+str('%-5d'%result1[ii][3])+' '+' '.join(str('%12e'%num) for num in result1[ii][2])+' '+\
                        str('%12e'%result1[ii][0])+' '+' '.join(str('%12e'%num) for num in result1[ii][4])+' '+str('%12e'%final_result )\
                        +' '+str('%12e'%result1[ii][6])+' '+str('%12e'%result1[ii][7])+' '+str('%12e'%result1[ii][1]))
                        
                with open(flag_locate,'a+') as flt:
                    flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% e,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,time_commu_end,'communication')) #this communication time include each computation         
                for rst in range(len(result1)):
                    Fselected[eval_initial+rst, 0] = result1[rst][0]
                    Time[eval_initial+rst, 0] = result1[rst][1]
                    xselected[eval_initial+rst,:]=result1[rst][2]
                    Cst [eval_initial+rst, :] = result1[rst][4] # Min add constraint
                    dummy_array=np.array(result1[rst][4])
                    dummy_array[dummy_array<0]=0
                    final_result = result1[rst][0]+penalty_coeff*sum(np.square(dummy_array))
                    Fselected_final[eval_initial+rst, 0] = final_result
                numevals=numevals+len(result1)
                time_opti_start1=time.time()
                ########### end parallel ####################
                '''
                pool=Pool()
                pool_res = pool.map_async(wrapper_func, ((i, data.objfunction) for i in xselected.tolist()))
                
                print 'pool_res.get()[0]: ',pool_res.get()[0],'\n','pool_res.get()[1]: ',pool_res.get()[1]
                pool.close()
                pool.join()
                result = pool_res.get()
                '''
                data.fevaltime[data.m:data.m+xselected.shape[0], 0] = Time
                data.S[data.m:data.m+xselected.shape[0], :] = xselected
                data.Y[data.m:data.m+xselected.shape[0], 0] = Fselected
                data.Y_final[data.m:data.m+xselected.shape[0], 0] = Fselected_final
                data.Cst [data.m:data.m+xselected.shape[0], :] = Cst # Min add constraint
                data.m = data.m + xselected.shape[0]
    
                # determine best one of newly sampled points
                minSelected = np.amin(Fselected_final)
                IDminSelected = np.argmin(Fselected_final)
                xMinSelected = xselected[IDminSelected, :]
                if minSelected < data.Fbest:
                    if data.Fbest - minSelected > (1e-3)*math.fabs(data.Fbest):
                        # "significant" improvement
                        failctr = 0
                        # count the number of  great success in this iteration
                        diff= data.Fbest - Fselected_final
                        num_succc_iter =len(diff [diff > (1e-3)*math.fabs(data.Fbest)]) 
                        succctr = succctr + num_succc_iter 
                    else:
                        #failctr = failctr + 1 used to be changed by Min for parallel version
                        failctr = failctr + len(result1)/2 # maybe half is okay as all my numbers are even, maybe add all number 
                        succctr = 0
                    data.xbest = xMinSelected
                    data.Fbest = minSelected
                else:
                    failctr = failctr + len(result1)/2 # maybe half is okay as all my numbers are even, maybe add all number 
                    succctr = 0
    
                # check if algorithm is in a local minimum
                shrinkflag = 1
                if failctr >= failtolerance:
                    if shrinkctr >= maxshrinkparam:
                        shrinkflag = 0
                        print('Stopped reducing sigma because the maximum reduction has been reached.')
                    failctr = 0
    
                    if shrinkflag == 1:
                        shrinkctr = shrinkctr + 1
                        sigma_stdev = sigma_stdev / 2
                        print('Reducing sigma by a half!')
                    else:
                        localminflag = 1
                        print( 'Algorithm is probably in a local minimum! Restarting the algorithm from scratch.')
    
                if succctr >= succtolerance:
                    sigma_stdev = min(2 * sigma_stdev, sigma_stdev_default)
                    succctr = 0
                # update PHI matrix only if planning to do another iteration
                if data.m < maxeval and localminflag == 0:
                    n_old = data.m - xselected.shape[0]
                    for kk in range(xselected.shape[0]):
                        new_phi = phi(normval[kk], data.phifunction)
                        #print new_phi.shape,new_phi
                        PHI[n_old + kk, 0: n_old + kk] = new_phi
                        PHI[0:n_old+kk, n_old+kk] = (new_phi).T
                        PHI[n_old+kk, n_old+kk] = phi0
                        P[n_old+kk, 1:data.dim+1] = xselected[kk, :]
                
                #Min---------generate "Flag.txt"file to locate where it was before the cut by Yellowstone after 12hours
                time_opti_end1=time.time()-time_opti_start1+time_opti_end
                with open(flag_locate,'a+') as flt:
                    flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% e,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,time_opti_end1,'optimization'))
                
                #np.savetxt('PHI.txt', PHI, fmt='%.10e', delimiter=' ', newline='\n')
                #np.savetxt('P.txt', P, fmt='%.10e', delimiter=' ', newline='\n')
                #np.savetxt('phi0.txt', phi0, fmt='%.10e', delimiter=' ', newline='\n')
                if (numevals/incre_resubmit) >= (incre*NumberNewSamples):
                    flag_resubmit=1
                    incre=incre+1
                    with open(flag_locate,'a+') as flt:
                        flt.write('\n%d,% d,% d,% d,% d,% d,% e,% d,% e,% d,% s' %(numevals,iterctr,shrinkctr,failctr,localminflag,succctr,sigma_stdev,pdim,phi0,incre,'incre'))
                    break
      
        data.S = data.S[0:data.m, :]
        data.Y = data.Y[0:data.m, :]
        data.Y_final = data.Y_final[0:data.m, :]
        data.Cst = data.Cst[0:data.m, :]
        data.fevaltime = data.fevaltime[0:data.m, :]
        data.NumberFevals = data.m   
    return data,localminflag,flag_resubmit
