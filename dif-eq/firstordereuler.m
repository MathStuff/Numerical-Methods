function [pred,real_sol,abs_err,rel_err] = firstordereuler(a,b,n)
    % Fill in the blanks (...) 
    % Applicable format : y' = f(t,y), y(a) = c, y = y(t), t>=a & t<=b
    % linear differential equation format :  y' + P(t)*y = Q(t)
    
    h = (b-a)/n;     % length of each sub-interval
    
    f = @(t,y) ... ;     % y' = f(t,y)
    
    pred = zeros(n+1,1);
    pred(1) = ... ;      % y(t=a)
    
    t = a:h:b;         % [a,b] divided into sub-intervals with length h
    
    for i = 2:n+1 
        prev = pred(i-1);
        pred(i) = prev + h*f(t(i-1),prev);  %y1 = y0 + h*f(t0,y0)
    end
    
    pred = pred';
    
    real_sol = ... ; % Analytic solution
    
    abs_err = abs(real_sol-pred);
    
    rel_err = abs(abs_err./real_sol)*100;
    
    end