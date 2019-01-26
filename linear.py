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

def Gauss_Siebel(coef_matrix,const_matrix,init_pred=None,iterations=2):
    """
    ***DOESN'T WORK***
    coef_matrix: Matrix that holds coefficients of the system
    const_matrix: Matrix that holds the constants/right side of the equation
    init_pred: Initial predictions for the unknowns, leave empty if not desired
    iterations: How many iterations to evaluate the roots for
    """
    def diagonalCheck(matr1,matr2):
        try:
            assert isinstance(matr1,Matrix) or isinstance(matr1,FMatrix)
        except:
            print("err")
            return False
        else:
            i=0
            while i<matr1.dim[0]:
                if matr1[i][i]==0:
                    print("0 on diagonal, swapping rows")
                    i2=i
                    while matr1[i2][i]==0:
                        i2+=1
                    temp=matr1[i][:]
                    matr1[i]=matr1[i2][:]
                    matr1[i2]=temp
                    
                    temp=matr2[i][0]
                    matr2[i][0]=matr2[i2][0]
                    matr2[i2][0]=temp
                    
                    return diagonalCheck(matr1,matr2)
                i+=1
            return True
        
    try:
        assert coef_matrix.dim[0]==coef_matrix.dim[1] and coef_matrix.dim[1]==const_matrix.dim[0] and const_matrix.dim[1]==1
        assert diagonalCheck(coef_matrix,const_matrix)==True

        if init_pred!=None:
            assert isinstance(init_pred,list)
            assert len(init_pred)==coef_matrix.dim[1]
        else:
            init_pred=[1]*coef_matrix.dim[1]
        iters=0

        while iters<iterations:

            for i in range(coef_matrix.dim[0]):
                num=0
                for j in range(coef_matrix.dim[0]):
                    if j==i:
                        continue
                    else:
                        num+=coef_matrix[i][j]*init_pred[j]
                init_pred[i]=(const_matrix[i][0]-num)/coef_matrix[i][i]
            iters+=1

    except Exception as err:
        print(err)
    else:
        for i in range(const_matrix.dim[0]):
            print("x{} =".format(i),init_pred[i])
        return FMatrix([const_matrix.dim[0],1],listed=[init_pred])

print("Data table:")
e1=linearEq(Matrix(5,ranged=[-5,5]),Matrix([5,1],ranged=[-15,10]))
e1.eq    

answerLU=LU(e1.coefs,e1.consts)

#answerGS=Gauss_Siebel(e1.coefs,e1.consts)