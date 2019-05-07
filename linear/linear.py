# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 22:18:23 2019

@author: Semih

Functions in this module:
    
-LU decomposition
-Gauss-Siedel

"""
from MatricesM.matrix import Matrix

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
                num=self.coefs.matrix[rows][unknowns]
                print(" "*(tab_size-len(str(num))+len(str(unknowns))-1),num,end="|")   
                
            print(" "*(tab_size-len(str(self.consts.matrix[rows][0]))),self.consts.matrix[rows][0])
        print("")
def LU(coef_matrix,const_matrix):
    """
    coef_matrix: Coefficient matrix of the linear system
    const_matrix: Right-hand side of the equations in a column matrix
    """
    try:
        assert coef_matrix.dim[0]==coef_matrix.dim[1] and coef_matrix.dim[1]==const_matrix.dim[0] and const_matrix.dim[1]==1
        
        L=coef_matrix.lowtri
        U=coef_matrix.uptri
        
        #UX=Y & LY=B
        X=Matrix([L.dim[0],1],fill=0,decimal=8)
        Y=Matrix([U.dim[0],1],fill=0,decimal=8)
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
            
            Y.matrix[ys][0] = B.matrix[ys][0]
            
            for rest in range(ys):
                Y.matrix[ys][0] -= L.matrix[ys][ys-rest-1]*Y.matrix[ys-rest-1][0]
        #Solve UX=Y using backward substitution
        for xs in range(1,Y.dim[0]+1):
            
            X.matrix[-xs][0] = Y.matrix[-xs][0]/U.matrix[-xs][-xs]
            
            for rest in range(1,xs):
                X.matrix[-xs][0] -= U.matrix[-xs][-rest] * X.matrix[-rest][0]/U.matrix[-xs][-xs]
                
        for i in range(X.dim[0]):
            print("x{} =".format(i),X.matrix[i][0])
        print("")
        
        return X

def Gauss_Siebel(coef_matrix,const_matrix,init_pred=None,iterations=25):
    """
    max(coef_matrix.eigenvalues)<1 IS REQUIRED FOR THIS METHOD TO CONVERGE THEREFORE WORK AS INTENDED
    coef_matrix: Matrix that holds coefficients of the system
    const_matrix: Matrix that holds the constants/right side of the equation
    init_pred: Initial predictions for the unknowns, leave empty if not desired
    iterations: How many iterations to evaluate the roots for
    """
    def diagonalCheck(matr1,matr2):
        try:
            assert isinstance(matr1,Matrix)
        except:
            print("err")
            return False
        else:
            i=0
            while i<matr1.dim[0]:
                if matr1.matrix[i][i]==0:
                    print("0 on diagonal, swapping rows")
                    i2=i
                    while matr1.matrix[i2][i]==0:
                        i2-=1
                    temp=matr1.matrix[i][:]
                    matr1.matrix[i]=matr1.matrix[i2][:]
                    matr1.matrix[i2]=temp
                    
                    temp=matr2.matrix[i][0]
                    matr2.matrix[i][0]=matr2.matrix[i2][0]
                    matr2.matrix[i2][0]=temp
                    
                    return diagonalCheck(matr1,matr2)
                i+=1
            return True
        
    try:
        assert coef_matrix.dim[0]==coef_matrix.dim[1] and coef_matrix.dim[1]==const_matrix.dim[0] and const_matrix.dim[1]==1
        assert diagonalCheck(coef_matrix,const_matrix)==True
        assert coef_matrix.det!=0
        
        if init_pred!=None:
            assert isinstance(init_pred,list)
            assert len(init_pred)==coef_matrix.dim[1]
        else:
            init_pred=[0]*coef_matrix.dim[1]
        iters=0
        
        while iters<iterations:
            
            for i in range(coef_matrix.dim[1]):
                num=0
                
                for j in range(coef_matrix.dim[1]):
                    if j==i:
                        continue
                    else:
                        num+=coef_matrix.matrix[i][j]*init_pred[j]
                        
                init_pred[i]=(const_matrix.matrix[i][0]-num)/coef_matrix.matrix[i][i]
                
            iters+=1
            
    except Exception as err:
        print(err)
    else:
        print("")
        total=0
        for i in range(const_matrix.dim[0]):
            print("x{} =".format(i),init_pred[i])
            total += init_pred[i]*coef_matrix.matrix[0][i]
            
        print("Error rate:{}%".format(abs((total-const_matrix.matrix[0][0])/total)*100))
        print("")
        
        return Matrix([const_matrix.dim[0],1],listed=[init_pred])

#EXAMPLES
e1=linearEq(Matrix(5,ranged=[-15,15]),Matrix([5,1],ranged=[-55,50]))
e1.eq   
LU(e1.coefs,e1.consts)

print("################################################\n")
e2=linearEq(Matrix(3,listed="4 -1 -1 -2 6 1 -1 1 7"),Matrix([3,1],listed="3 9 -6"))
e2.eq

print("LU answer:")
answerLU=LU(e2.coefs,e2.consts)
print("GS estimation after 5 iterations:")
answerGS=Gauss_Siebel(e2.coefs,e2.consts,iterations=5)
