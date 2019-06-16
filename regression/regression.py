# -*- coding: utf-8 -*-
"""
Created on Sat Jan 26 13:09:52 2019

@author: Semih

Models in this module:
    
-Linear regression (intercept & no intercept)
-Power regression
-Exponential regression
    
"""
from MatricesM import *

def linear_reg(xs,ys):
    """
    xs: Matrix; Values to use for prediction
    ys: Matrix; Values to be predicted
    Returns in a tuple -> (Coefficients,Predictions,Errors)
    """
    #Coefficients
    coefs = (((xs.t@xs).inv)@xs.t)@ys
    #Re-name coefficient matrix's columns and also dtype (to display column names)
    coefs.dtype = dataframe
    coefs.features = [ys.features[i]+" coefficient" for i in range(coefs.dim[1])]
    
    #Predict values
    preds = xs@coefs
    #Re-name prediction matrix's columns
    preds.features = [ys.features[i]+" prediction" for i in range(preds.dim[1])]
                      
    #Errors
    err = ys-preds
    #Re-name error matrix's columns
    err.features = [ys.features[i]+" error" for i in range(err.dim[1])]
    
    return (coefs,preds,err)
