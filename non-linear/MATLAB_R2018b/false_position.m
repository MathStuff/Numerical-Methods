function root = false_position(fun,low,up,iters)
%Regula-falsi(False Position) method to find an estimate root of a function
%fun:funtion
%low:lower bound
%up:upper bound
%iters:iteration amount
%EXAMPLE:
%   false_position(@(a) a*a*a+5*a*a-7*a-5,1,3,14) == 1.5572

    %%%CHECK_INTERVAL%%%
    if fun(low)*fun(up)>0
        fprintf("Not in the interval\n");
        iters=0;
    end
    
    while iters>0
        %%%ESTIMATE_ROOT%%%
        root = up - (fun(up)*(low-up)/(fun(low)-fun(up)));
        
        %%%CHANGE_INTERVAL%%%
        if fun(low)*fun(root)<0
            up = root;
        elseif fun(low)*fun(root)>0
            low = root;
        %%%FOUND_ROOT%%%
        else
            iters = 1;
        end
        
        iters = iters-1;
    end
end