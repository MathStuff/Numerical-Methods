function root = bisection(fun,low,up,iters)
%Bisection method to find the root of a function
%   fun:function
%   low:lower bound
%   upper:upper bound
%   iters:iteration amount  
%   EXAMPLE:
%       bisection_answer=bisection(@(a) a*a*a-29,0,4,16) == 3.0723
    %%%%%CHECK_INTERVAL%%%%%
    if fun(low)*fun(up)>0
        fprintf("Not in the interval\n");
        iters=0;
    end

    while iters>0
        %%%%%FOUND%%%%%%
        if fun(low)==0
            root=low;
            break
        end

        if fun(up)==0
            root=up;
            break
        end
        %%%%%%%%%%%%%%%%

        %%%%ESTIMATE%%%%
        root=(low+up)/2;

        if fun(low)*fun(root)<0
            up=root;
        end

        if fun(root)*fun(up)<0
            low=root;
        end
        iters=iters-1;
        %%%%%%%%%%%%%%%%
    end
end