import numpy as np
import struct
'''
# read binary file from MODFLOW based on
The array data can be either in single-precision or double-precision format.
There is no data field that indicates which format is used.  Instead, you will need to try both and see which one works.
First read the following variables in order
KSTP: the time step number, an integer, 4 bytes.
KPER: the stress period number, an integer, 4 bytes.
PERTIM: the time in the current stress period, a real number, either 4 or 8 bytes.
TOTIM, the total elapsed time, a real number, either 4 or 8 bytes.
DESC, a description of the array, 16 ANSI characters, 16 bytes.
NCOL, the number of columns in the array, an integer, 4 bytes.
NROW, the number of rows in the array, an integer, 4 bytes.
ILAY, the layer number, an integer, 4 bytes.
Next come a list of NROW x NCOL real numbers that represent the values of the array.
The values are in row major order.  Each value in the array occupies either 4 or 8 bytes
depending on whether the values are in single- or double-precision.
After reading one set of values, start over with KSTP. Continue until reaching the end of the file.
The following is a list of possible values for DESC.  The list may be incomplete and is subject to change.
Check the values passed to the subroutine ULASAV in the MODFLOW source code for other possible values.
'            HEAD'
'        DRAWDOWN'
'      SUBSIDENCE'
'      COMPACTION'
'   CRITICAL HEAD'
'     HEAD IN HGU'
'NDSYS COMPACTION'
'  Z DISPLACEMENT'
' D CRITICAL HEAD'
'LAYER COMPACTION'
' DSYS COMPACTION'
'ND CRITICAL HEAD'
'LAYER COMPACTION'
'SYSTM COMPACTION'
'PRECONSOL STRESS'
'CHANGE IN PCSTRS'
'EFFECTIVE STRESS'
'CHANGE IN EFF-ST'
'      VOID RATIO'
'       THICKNESS'
'CENTER ELEVATION'
'GEOSTATIC STRESS'
'CHANGE IN G-STRS'
One way to determine whether the file has been saved with single- or double-precision,
is to read the file up through DESC using either single- or  double-precision numbers and
see if the value read for DESC matches one of the above values.
'''

def read_binary_file(f,nrow,ncol,nsp):
    ffile=open(f,'rb')
    result=np.zeros([nsp*3,nrow,ncol])
    for i in range(4):
        ff=struct.unpack('s',ffile.read(1))[0]
    for j in range(nsp*3):
        #print 'i',i
        first=struct.unpack('i',ffile.read(4))[0]
        #print 'dont know',first
        KSTP=struct.unpack('i',ffile.read(4))[0]
        #print kstp
        KPER=struct.unpack('i',ffile.read(4))[0]
        #print pertim
        TIME2=struct.unpack('f',ffile.read(4))[0]
        TEXT=''
        for i in range(16):
            TEXT += struct.unpack('s',ffile.read(1))[0]
        #print totaltime
#         ffile.seek(20,1) 
#         print(TEXT)
        ncol=struct.unpack('i',ffile.read(4))[0]
#         print ncol
        nrow=struct.unpack('i',ffile.read(4))[0]
#         print nrow
        nlay=struct.unpack('i',ffile.read(4))[0]
        xy=struct.unpack('i',ffile.read(4))[0]
        xyy = struct.unpack('i',ffile.read(4))[0]
#         print 'after lay',xyy
#         p=struct.unpack('i',ffile.read(4))[0]
#         print 'after lay',p
#         q=struct.unpack('i',ffile.read(4))[0]
#         print 'after lay',q
#         ffile.seek(12,1)
        ncell=ncol*nrow

        temp=np.zeros(ncell)
        for ij in range(ncell):
            temp[ij]=struct.unpack('f',ffile.read(4))[0]
        result_plan = temp.reshape(nrow,ncol)
        result[j,:,:]=result_plan
        
        if j!=nsp*3-1:
            for ijj in range(8):
                t=struct.unpack('s',ffile.read(1))[0]
    length=ffile.tell()
    
    #print length
    ffile.close()
    '''
    for i in range(nrow):
        for j in range(ncol):
            if result[nsp*3,]
    '''
    max_sub=np.max(result)
    i,j,k=np.unravel_index(result.argmax(),result.shape) # in order to find the place of the maximum subsidence
#     print max_sub
#     import pdb 
#     pdb.set_trace()
#     print i
#     print j
#     print k
    return result
        
        