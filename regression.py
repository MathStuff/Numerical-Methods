# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 13:09:52 2019

@author: Semih

Models in this module:
    
-Linear regression (intercept & no intercept)
-Power regression
-Exponential regression
    
"""
class dataTable:
    def __init__(self,xs=None,ys=None,Dict=None):
        self.xs=xs
        self.ys=ys
        self.Dict=Dict
    
    @property
    def p(self):
        if self.Dict!=None:
            pass
        else:
            pass
                