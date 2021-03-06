function [pred] = newtons_method(pred,hessian,grad)
    % Newton's algorithm for minimizing functions
    
    % pred is a vector of initial guesses for each variable as rows
    % hessian is the Hessian matrix
    % grad is the gradient of the function
    
    % gradient_descent(guess,H,[wrt_dx1;wrt_dx2;...])
    % Example of a function with 3 variables: "x(1)^3/3 + x(2)^2 / 2 - x(3)";
    
    %FUNCTIONS SHOULD BE PASSED AS STRINGS
    %EXAMPLE:
    %
    %   [pred] = newtons_method([1;1],[4 2;2 4],["4*x(1)+2*x(2)-4" ;"2*x(1)+4*x(2)-6"])
    %
    
    temp = size(pred);
    amount_of_vars = temp(1,1);
    x = sym('x', [amount_of_vars 1]);

    %Create the Hessian matrix with functions
    H = cell(amount_of_vars);
    
    col = 1;
    for i=1:amount_of_vars^2
        md = mod(i,amount_of_vars);
        func = eval(string(hessian(md+1,col)+"+0*x(1)"));
        H{i} = matlabFunction(func , 'vars', {x});
        if not(md)
            col = col + 1;
        end
    end
    
    %Create the gradiant vector with functions
    gradF = cell(amount_of_vars,1);
    
    for i=1:amount_of_vars
        func = eval(grad(i));
        gradF{i} = matlabFunction(func , 'vars', {x});     
    end

    g = zeros([amount_of_vars 1]);
    Q = zeros([amount_of_vars amount_of_vars]);
    
    % When to stop for gradient
    c = 1e-8;
    no_des = ones([amount_of_vars 1]).*c;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Calculate gradient with guess
    for i=1:amount_of_vars
        func = gradF{i};
        g(i) = func(pred);
    end
    
    %Calculate Hessian with guess
    col = 1;
    for i=1:amount_of_vars^2
        md = mod(i,amount_of_vars);
        func = H{i};
        Q(md+1,col) = func(pred);
        if not(md)
            col = col + 1;
        end
    end

    %Update guess
    pred = pred - inv(Q)*g;
    
    %Condition to decide gradient is 0
    cond = false;
    for i=1:amount_of_vars
        if abs(g(i))>no_des(i)
            cond = true;
            break
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Loop until 0 vector gradient
    while cond
        
        %Calculate gradient with guess
        for i=1:amount_of_vars
            func = gradF{i};
            g(i) = func(pred);
        end
        
        %Calculate hessian with guess
        col = 1;
        for i=1:amount_of_vars^2
            md = mod(i,amount_of_vars);
            func = H{i};
            Q(md+1,col) = func(pred);
            if not(md)
                col = col + 1;
            end
        end

        %Update guess
        pred = pred - inv(Q)*g;
        
        %Update for stopping
        for i=1:amount_of_vars
            if abs(g(i))>no_des(i)
                cond = true;
                break
            end
            cond = false;
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end