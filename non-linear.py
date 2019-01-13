# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 19:21:54 2018

@author: Semih
"""

def derivative(func,n,error_percentage=1e-10):
    return (func(n+error_percentage)-func(n))/error_percentage

def integral(func,a,b,segments=1e+5):
    inc=(b-a)/segments
    i=a+inc
    total=0
    while i<=b:
        total+=abs(func(i))*inc
        i+=inc
    return round(total,4)
     
# SOLUTIONS FOR NON-LINEAR EQUATION

def bisection(func,lower,upper,error_upper_bound=None,i=20):
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
            
    print("Error rate:{}%".format(e))
    return guess

#print(bisection(lambda a: (a*a)-5,0,4,0.01))

def newton_raphson():
    pass