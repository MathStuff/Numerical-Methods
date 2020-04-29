function [pred,real_sol,abs_err,rel_err] = higherordereuler(a,b,n)
    % Fill in the blanks (...) 
    % y(a) = c1, y'(a) = c2, y''(a) = c3, ..., y = y(t), t>=a & t<=b

    h = (b-a)/n;     % length of each sub-interval

    %Add more functions if needed
    f1 = @(t,x,z) ... ;     % z  = y'
    f2 = @(t,x,z) ... ;     % z' = y'' = f(t,y)
    
    %Increment MANUALLY
    func_amount = 2;  %how many functions were used
        
    pred = zeros(n+1,func_amount);
    
    %Add more initial values if needed
    pred(1,1) = ... ;      % y (t=a)
    pred(1,2) = ... ;      % z = y'(t=a)
    
    t = a:h:b;         % [a,b] divided into sub-intervals with length h
    
    %Add more calculations if necessary
    for i = 2:n+1 
        
        pred(i,2) = pred(i-1,2) + h*f2(t(i-1),pred(i-1,1),pred(i-1,2));  %z1 = z0 + h*f(t0,y0,z0)
        
        pred(i,1) = pred(i-1,1) + h*f1(t(i-1),pred(i-1,1),pred(i-1,2));  %y1 = y0 + h*f(t0,y0,z0)
        
    end
    
    pred = (pred(:,1:1))'; % get the actual prediction columns and transpose
    
    real_sol = ... ; % Analytic solution
    
    abs_err = abs(real_sol-pred);
    
    rel_err = abs(abs_err./real_sol)*100;
    
    end
