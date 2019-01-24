# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 22:18:23 2019

@author: Semih

Functions in this module:
    
-LU decomposition
-Gauss-Siedel

"""
from matrices import Matrix,FMatrix,Identity

class linearEq:
    """
    coefs: Coefficients' matrix
    consts: Constants' matrix
    """
    def __init__(self,coefs,consts):
        try:
            assert isinstance(coefs,Matrix) and isinstance(consts,Matrix)
        except:
            return None
        else:
            self.coefs=coefs
            self.consts=consts
            
    @property
    def eq(self):
        """
        Print out the linear system of equations
        """
        temp=[]
        for r in self.coefs.matrix:
            for els in r:
                temp.append(len(str(els)))
                
        tab_size=max(temp)+2
        
        for i in range(self.coefs.dim[1]):
            print(" "*(tab_size-2),"x{0}".format(i),end="|")
            
        print(" "*(tab_size//2 - 3),"CONSTANTS")
        print("*"*((tab_size+4)*self.coefs.dim[1]))
        
        for rows in range(self.coefs.dim[0]):
            
            for unknowns in range(self.coefs.dim[1]):
                num=self.coefs[rows][unknowns]
                print(" "*(tab_size-len(str(num))),num,end="|")   
                
            print(" "*(tab_size-len(str(self.consts[rows][0]))),self.consts[rows][0])
        
e1=linearEq(Matrix(8,ranged=[-5,5]),Matrix([8,1],ranged=[-25,50]))
e1.eq

def LU(coef_matrix,const_matrix):
    pass

def Gauss_Siebel(coef_matrix,const_matrix,init_pred=None,error_upper_bound=1e-2):
    """
    coef_matrix: Matrix that holds coefficients of the system
    const_matrix: Matrix that holds the constants/right side of the equation
    init_pred: Initial predictions for the unknowns, leave empty if not desired
    error_upper_bound: Relative error rate's upper bound to stop the iteration
    """
    pass
    

