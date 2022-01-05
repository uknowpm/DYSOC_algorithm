import matplotlib.pyplot as plt
import numpy as np
#from utility import *
import sys
import os
import itertools
from os import listdir
import shutil
import fileinput
import sys

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)
    

def main():
    #"change the file in 5 places (mainly change the folder name) in order to reinitialize the code "
    mypath=os.getcwd()
    parepath=os.path.dirname(mypath)
    folder= mypath.replace(parepath,"")[1:]
    print folder
    if "Umatilla" in folder:
        name="Umatilla"
    elif "Blaine" in folder:
        name="Blaine"
    else:
        print "it's not a Umatilla or Blaine problem"
    #get processor number
    nproc=folder.split("_")[1]
    print nproc
    
    # change the StochasticRBF.py line 24 (directory)
    lines_SRBF=file(mypath+'/StochasticRBF.py','r').readlines()
    print lines_SRBF[24]
    k=lines_SRBF[24].split('\n')[0]
    replaceAll(mypath+'/StochasticRBF.py',k,'    os.chdir(os.getcwd()'+'+\'/'+folder+'\')')
    
    # change each BlaineTest or UmatillaTest.sh file
    if name=="Umatilla":
        test_sh='UmatillaTest.sh'
        print int(nproc)>16
        if int(nproc)>16:
            maxeval = 2000
        else:
            maxeval = 1500
        if int(nproc)==1:
            maxeval = 1000
    elif name=="Blaine":
        test_sh='BlaineTest.sh'
        if int(nproc)>16:
            maxeval = 600
        else:
            maxeval = 400
        if int(nproc)==1:
            maxeval = 400
    lines_SH=file(mypath+'/'+test_sh,'r').readlines()
    #set processor number
    proc=lines_SH[9].split('\n')[0]
    replaceAll(mypath+'/'+test_sh,proc,'set NumberNewSamples='+nproc)
    meval=lines_SH[6].split('\n')[0]
    replaceAll(mypath+'/'+test_sh,meval,'set maxeval='+str(maxeval))
    NeedR=lines_SH[10].split('\n')[0]
    replaceAll(mypath+'/'+test_sh,NeedR,'set NeedRe=0')
    k=lines_SH[12].split('\n')[0]
    replaceAll(mypath+'/'+test_sh,k,'set st_pt1=0')
    k=lines_SH[13].split('\n')[0]
    replaceAll(mypath+'/'+test_sh,k,'set st_pt2=0')
    direct1=lines_SH[15].split('\n')[0]
    replaceAll(mypath+'/'+test_sh,direct1,'set PYTHONPATH=/glade/u/home/minpang/'+folder)
    direct2=lines_SH[18].split('\n')[0]
    replaceAll(mypath+'/'+test_sh,direct2,\
        'python -m cProfile -o timing.txt /glade/p/work/minpang/'+folder+'/DYCORS.py ${data_file} $maxeval $Ntrials $PlotResult $NumberNewSamples $NeedRe $FinalRun $st_pt1 $st_pt2')
    
    #change local_restart and local_restart_re
    # local_restart
    lines_lre=file(mypath+'/LocalStochRBFrestart.py','r').readlines()
    a=lines_lre[121].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart.py',a,'            lines=file(\''+test_sh+'\',\'r\').readlines()')
    b=lines_lre[125].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart.py',b,'            replaceAll(\''+test_sh+'\',NeedRe_line,\'set NeedRe=\'+str(1))')
    c=lines_lre[126].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart.py',c,'            replaceAll(\''+test_sh+'\',st_pt1_line,\'set st_pt1=\'+str(st_pt1))')
    d=lines_lre[127].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart.py',d,'            replaceAll(\''+test_sh+'\',st_pt2_line,\'set st_pt2=\'+str(st_pt2))')
    e=lines_lre[138].split('\n')[0]
    if int(nproc) <16:
        replaceAll(mypath+'/LocalStochRBFrestart.py',e,'            pathSubmitJob = [\"bsub\", \"-P\", \"UCOR0001\", \"-W\",\"12:00\",\"-n\",\"'+str(nproc)+'\",\"-R\", \"\\\"span[ptile='+str(nproc)+']\\\"\",\\')
    else:
        replaceAll(mypath+'/LocalStochRBFrestart.py',e,'            pathSubmitJob = [\"bsub\", \"-P\", \"UCOR0001\", \"-W\",\"12:00\",\"-n\",\"'+str(nproc)+'\",\"-R\", \"\\\"span[ptile='+'16'+']\\\"\",\\')
    f1=lines_lre[139].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart.py',f1,'                             \"-J\",\"myjob\",\"-o\",\"/glade/p/work/minpang/'+folder+'/myjob.%J.out \",\"-e\",\"/glade/p/work/minpang/'+folder+'/myjob.%J.err\",\"-q\", \"regular\",\\')
    f=lines_lre[140].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart.py',f,'                             "/glade/p/work/minpang/'+folder+'/'+test_sh+'"]')
    #change local_restart_re
    lines_lrere=file(mypath+'/LocalStochRBFrestart_re.py','r').readlines()
    a_re=lines_lrere[158].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart_re.py',a_re,'            lines=file(\''+test_sh+'\',\'r\').readlines()')
    c_re=lines_lrere[163].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart_re.py',c_re,'            replaceAll(\''+test_sh+'\',st_pt1_line,\'set st_pt1=\'+str(st_pt1))')
    d_re=lines_lrere[164].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart_re.py',d_re,'            replaceAll(\''+test_sh+'\',st_pt2_line,\'set st_pt2=\'+str(st_pt2))')
    e_re=lines_lrere[175].split('\n')[0]
    if int(nproc) <16:
        replaceAll(mypath+'/LocalStochRBFrestart_re.py',e_re,'            pathSubmitJob = [\"bsub\", \"-P\", \"UCOR0001\", \"-W\",\"12:00\",\"-n\",\"'+str(nproc)+'\",\"-R\", \"\\\"span[ptile='+str(nproc)+']\\\"\",\\')
    else:
        replaceAll(mypath+'/LocalStochRBFrestart_re.py',e_re,'            pathSubmitJob = [\"bsub\", \"-P\", \"UCOR0001\", \"-W\",\"12:00\",\"-n\",\"'+str(nproc)+'\",\"-R\", \"\\\"span[ptile='+'16'+']\\\"\",\\')
    f1_re=lines_lrere[176].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart_re.py',f1_re,'                             \"-J\",\"myjob\",\"-o\",\"/glade/p/work/minpang/'+folder+'/myjob.%J.out \",\"-e\",\"/glade/p/work/minpang/'+folder+'/myjob.%J.err\",\"-q\", \"regular\",\\')
    f_re=lines_lrere[177].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFrestart_re.py',f_re,'                             "/glade/p/work/minpang/'+folder+'/'+test_sh+'"]')
    
    #change local_stop and local_stop_re
    if name=="Umatilla":
        incre_submit=201
    elif name=="Blaine":
        incre_submit=15
    lines_ls=file(mypath+'/LocalStochRBFstop.py','r').readlines()
    a=lines_ls[72].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFstop.py',a,'    incre_resubmit='+str(incre_submit))
    lines_ls=file(mypath+'/LocalStochRBFstop_re.py','r').readlines()
    a=lines_ls[75].split('\n')[0]
    replaceAll(mypath+'/LocalStochRBFstop_re.py',a,'    incre_resubmit='+str(incre_submit))
    
    
    
if __name__ == "__main__":
    main()
 