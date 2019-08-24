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

def linear_reg(x,y,intercept=False):
    """
    xs: Matrix; Values to use for prediction
    ys: Matrix; Values to be predicted
    intercept: bool: Wheter or not to add a constant to the formula
    Returns in a tuple -> (Coefficients,Predictions,Errors,Function)
    """
    from string import ascii_letters,digits
    
    xs = x.copy
    ys = y.copy
    
    feats = xs.features[:]
    validchars = ascii_letters+digits+"_"
    
    for i,name in enumerate(feats[:]):
        for char in name[:]:
            if not char in validchars:
                name = name.replace(char,"_")
        feats[i] = name
    #Coefficients
    if intercept:
        xs.add([1 for _ in range(xs.d0)],col=1,feature="constant")
    coefs = (((xs.t@xs).inv)@xs.t)@ys
    #Re-name coefficient matrix's columns, row labels as coefficient names
    coefs.features = [ys.features[i]+" coefficient" for i in range(coefs.dim[1])]
    coefs.index = xs.features[:]
    
    #Predict values
    preds = xs@coefs
    #Re-name prediction matrix's columns
    preds.features = [ys.features[i]+" prediction" for i in range(preds.dim[1])]
                      
    #Errors
    err = ys-preds
    #Re-name error matrix's columns
    err.features = [ys.features[i]+" error" for i in range(err.dim[1])]

    #Create the formula
    params = ",".join(feats)
    right_h_s = ""
    coef_list = coefs.col(1,0)
    
    if intercept:
        right_h_s += str(coef_list[0])
        del coef_list[0]
        
    for i,coef in enumerate(coef_list):
        right_h_s += "+(" + str(coef) + f")*{feats[i]}"
        
    formula = eval(f'lambda {params}:{right_h_s}')
    
    return (coefs,preds,err,formula,right_h_s)
