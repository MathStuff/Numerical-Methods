# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 13:09:52 2019

@author: Semih

Models in this module:
    
-Linear regression (intercept & no intercept)
-Power regression
-Exponential regression
    
"""
from matrices import Matrix,FMatrix,Identity

class dataTable:
    """
    xs:list if only 1 variable is given, matrix or dictionary if multiple variables are given
    ys:list or column matrix of real numbers
    Dict: data in a dictionary; keys are the variables except the last key being the value to predict
        Example = N columns of data to be treated as :
            {var1:[real numbers], var2:[real numbers], ... , var(N-1):[real numbers],value_to_pred:[real numbers]}
            
    """
    def __init__(self,xs=None,ys=None,Dict=None):
        try:
            self.xs=xs
            self.ys=ys
            self.Dict=Dict
            if isinstance(xs,list):
                if isinstance(ys,list):
                    self.ys=FMatrix([len(ys),1],listed=ys)    
                
                assert len(xs)
                self.xs
    @property
    def p(self):
        if self.Dict!=None:
            pass
        else:
            pass
                