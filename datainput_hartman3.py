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
from utility import *
from numpy import *
from numpy.matlib import *
from operator import mul

def datainput_hartman3():
    data = Data()
    data.xlow = zeros(3) #variable lower bounds
    data.xup = ones(3) #variable upper bounds
    data.objfunction = myfun #objective function -must return scalar
    data.dim = 3 #problem dimension
    data.Cst_no=2
    return data

def myfun(x):
    x=np.array(x)
    folder = x[-1]
    x=x[:-1]
#     folder = x[-1]
#     x=x[:-1]
    #define your objective function here
    c = array([1, 1.2, 3, 3.2])
    A = array([[3, 10, 30], [0.1, 10, 35], [3, 10, 30], [0.1, 10, 35]])
    P = array([[0.3689, 0.1170, 0.2673], 
        [0.4699, 0.4387, 0.747], 
        [0.1091, 0.8732, 0.5547],
        [0.0382, 0.5743, 0.8828]])
    y = -sum(c * exp(-sum(A * (repmat(x, 4, 1) - P) ** 2, axis = 1)))
    cts = sum(x)-0.5
    cts1 = 0.9-sum(x)
    return y,[cts,cts1],cts,cts1,c[0],c[1],c[0]

if __name__ == '__main__':
    print( myfun(array([0.5, 0.9, 0.3,1])))
