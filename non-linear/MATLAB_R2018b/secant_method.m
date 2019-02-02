function root = secant_method(fun,g1,g2,iters)
%Newton-Raphson method to find an estimation of the root of a function
%fun:function
%g1:initial guess 1
%g2:initial guess 2
%iters:iteration amount
%EXAMPLE:
%   secant_method(@(a) a*a*a+5*a*a-7*a-5,3,5,7) == 1.5572
    while iters>0
        root = g2 - (fun(g2)*(g2-g1)/(fun(g2)-fun(g1)));
        g1 = g2;
        g2 = root;
        iters = iters-1;
    end
end
