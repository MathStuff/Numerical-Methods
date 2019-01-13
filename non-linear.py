# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 19:21:54 2018

@author: Semih
"""

def derivative(func,n,error_percentage=1e-10):
    """
    Limit definition of derivative
    """
    return (func(n+error_percentage)-func(n))/error_percentage

def integral(func,a,b,seg=1e+5):
    """
    Riemann integral with an error rate of ~(1/seg)%
    func: function
    a: lower limit
    b: upper limit
    seg: amount of segments to split the interval
    """
    inc=(b-a)/seg
    i=a
    total=0
    while i<b:
        total+=abs(func(i))*inc
        i+=inc
    return round(total,4)

#print(derivative(lambda a:a*a -4,3))
#print(integral(lambda a:a*a -4, 2, 12)) 
    
""" SOLUTIONS FOR NON-LINEAR EQUATION """
    
# =============================================================================
# BISECTION METHOD 
# =============================================================================
def bisection(func,lower,upper,error_upper_bound=None,i=100):
    """
    func: function
    lower: lower bound
    upper: upper bound
    error_upper_bound: relative error's upper limit, stops the function when the relative error is lower than this parameter
    i: iteration amount, if relative error upper bound is already reached, iterations stop
    """
    #Initiate values
    interval=[min(lower,upper),max(lower,upper)]
    e=abs((interval[1]-interval[0])*100/interval[1])
    
    #Correct error bounds
    if error_upper_bound==None:
        error_upper_bound=0
        
    #Start iterating until the bound is achieved or iterations are done
    while i>0 and not e<=error_upper_bound:
        try:
            l=func(interval[0])
            u=func(interval[1])
            guess=(sum(interval))/2
            g=func(guess)
            if 0 in (l,g,u):
                ind=[l,g,u].index(0)
                return [lower,guess,upper][ind]
            elif l*g<0:
                i-=1
                interval[1]=guess
            elif g*u<0:
                i-=1
                interval[0]=guess
            else:
                #print("No root in given range!")
                return None
        except:
            print("Bad arguments")
        else:
            #Recalculate the relative error
            e=abs((interval[1]-interval[0])*100/interval[1])
            
    print("Relative error rate:{}%".format(e))
    return guess

#print(bisection(lambda a: (a*a)-5,0,4,1e-5)

# =============================================================================
# NEWTON-RAPHSON METHOD
# =============================================================================
def newton_raphson(func,init,error_upper_bound=None,derivative_increment=1e-10):
    """
    func: function
    init: initial guess for the root
    error_upper_bound: relative error's upper limit, function stops when relative error is smaller than this parameter
    derivative_increment: amount of incrementation while taking derivative by limit definition, for better results decrease the value 
    """
    try:
        #Assertion of the parameters
        assert isinstance(error_upper_bound,float) or isinstance(error_upper_bound,int)
    except AssertionError:
        print("Error upper bound must be a float or integer value")
        return None
    else:
        try:
            #Initial calculations 
            e=100
            prev=init-(func(init))/derivative(func,init,error_percentage=derivative_increment)
            g1=func(prev)
            
            #Start the loop and don't stop until root is found or the relative error is small enough
            while e>error_upper_bound and g1!=0:
                guess=prev-(func(prev)/derivative(func,prev,error_percentage=derivative_increment))
                g1=func(guess)
                e=abs((guess-prev)*100/guess)
                prev=guess
        #If the incrementation amount, while taking derivative, was too small increase the value    
        except ZeroDivisionError:
            print("Lowering derivative incrementation amount")
            return newton_raphson(func,init,error_upper_bound,derivative_increment=derivative_increment*10)
        else:
            print("Error rate:{}%".format(e))
            return prev

#print(newton_raphson(lambda a:(a**4)-4*(a**3)-2*a+6, 0, 1e-4, 1e-30))         

# =============================================================================
# SECANT METHOD           
# =============================================================================
def secant():
    pass        
            
    
    




