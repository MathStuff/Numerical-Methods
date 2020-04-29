function [pred,real_sol,abs_err,rel_err] = milnesimpson(a,b,n)
    % Fill in the blanks (...) 
    % Applicable format : y' = f(t,y), y(a) = c, y = y(t), t>=a & t<=b
    % (n >= 4) <=> (b-a/h >= 4)
    if (n<4)
        error("This method requires at least 4 sub-intervals");
    end
    
    h = (b-a)/n;     % length of each sub-interval
    
    f = @(t,y) ...;     % y' = f(t,y)
    
    pred = zeros(n+1,1);
    pred(1) = ...;      % y(t=a)
    
    t = a:h:b;         % [a,b] divided into sub-intervals with length h
    
    
    %RK-4 for initial values
    for i = 2:4
        prev_y = pred(i-1);
        prev_t = t(i-1);
        
        k1 = f(prev_t,prev_y);
        k2 = f(prev_t+(h/2),prev_y+(h/2)*k1);
        k3 = f(prev_t+(h/2),prev_y+(h/2)*k2);
        k4 = f(prev_t+h,prev_y+h*k3);
        
        pred(i) = prev_y + (h/6)*(k1+2*k2+2*k3+k4);
    end
    
    %Milne predictions and Simpson corrections
    for i = 5:n+1
        f1 = f(t(i-3),pred(i-3));
        f2 = f(t(i-2),pred(i-2));
        f3 = f(t(i-1),pred(i-1));
        p = pred(i-4) +  (4*h/3)*(2*f1 - f2 + 2*f3);
        
        pred(i) = pred(i-2) + (h/3)*(f(t(i-2),pred(i-2))+4*f(t(i-1),pred(i-1))+f(t(i),p));
    end
    
    pred = pred';
    
    real_sol = ...; % Analytic solution
    
    abs_err = abs(real_sol-pred);
    
    rel_err = abs(abs_err./real_sol)*100;
    
    end