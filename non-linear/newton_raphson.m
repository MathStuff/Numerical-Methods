function root = Newton_Raphson(fun,init,iters)
%Newton-Raphson method to find an estimation of the root of a function
%fun:function
%init:initial guess
%iters:iteration amount
%EXAMPLE:
%   newton_raphson(@(a) a*a*a+5*a*a-7*a-5, 4, 5) == 1.5572
%   newton_raphson(@(a) a*a*a+5*a*a-7*a-5, 0, 5) == -0.5330
%   newton_raphson(@(a) a*a*a+5*a*a-7*a-5, -4, 10) == -6.0242

    function num = derivative(fun,x)
    %Derivative of a function around some point x 
        num = (fun(x+1e-5)-fun(x))/(1e-5);
    end

    while iters>0
        root = init - (fun(init)/derivative(fun,init));
        init = root;
        iters = iters-1;
    end
end
