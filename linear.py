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


def LU(coef_matrix,const_matrix):
    try:
        assert coef_matrix.dim[0]==coef_matrix.dim[1] and coef_matrix.dim[1]==const_matrix.dim[0] and const_matrix.dim[1]==1
        L=coef_matrix.lowtri.copy
        U=coef_matrix.uptri.copy
        #UX=Y & LY=B
        X=FMatrix(dim=[L.dim[0],1],randomFill=0,decimal=8)
        Y=FMatrix(dim=[U.dim[0],1],randomFill=0,decimal=8)
        B=const_matrix.copy
    except AttributeError:
        print("System has infinite solutions")
        return "Inf"
    except AssertionError:
        print("Matrices should be in nxn and nx1 dimensions")
    else:
        print("")
        #Solve LY=B using forward substitution
        for ys in range(B.dim[0]):
            Y[ys][0]=B[ys][0]
            for rest in range(ys):
                Y[ys][0]-=L[ys][ys-rest-1]*Y[ys-rest-1][0]
        #Solve UX=Y using backward substitution
        for xs in range(1,Y.dim[0]+1):
            X[-xs][0]=Y[-xs][0]/U[-xs][-xs]
            for rest in range(1,xs):
                X[-xs][0]-=U[-xs][-rest]*X[-rest][0]/U[-xs][-xs]
        for i in range(X.dim[0]):
            print("x{} =".format(i),X[i][0])
        return X

print("Data table:")
e1=linearEq(Matrix(5,ranged=[-5,5]),Matrix([5,1],ranged=[-15,10]))
e1.eq    
answer=LU(e1.coefs,e1.consts)

def Gauss_Siebel(coef_matrix,const_matrix,init_pred=None,error_upper_bound=1e-2):
    """
    coef_matrix: Matrix that holds coefficients of the system
    const_matrix: Matrix that holds the constants/right side of the equation
    init_pred: Initial predictions for the unknowns, leave empty if not desired
    error_upper_bound: Relative error rate's upper bound to stop the iteration
    """
    pass
    

