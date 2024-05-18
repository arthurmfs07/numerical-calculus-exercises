x0 = input('');
m = input('');
s = input('');

A = lotka_leslie(m, s);
[x, k] = potencia(A, x0, 1e-6);


disp('x = ');
disp(x);
disp(k);

function [A] = lotka_leslie(m, s)
    n = length(m);
    A = zeros(n);
    A(1,:) = m';
    for i = 2:n
        A(i,i-1) = s(i-1);
    end
end

function [x, k] = potencia(A, x0, tol)
    x = x0 / norm(x0);
    k = 0;
    erro = inf;
    
    while erro > tol
        x_1 = A * x;
        x_1 = x_1 / norm(x_1);
        
        cos = dot(x, x_1) / (norm(x) * norm(x_1));
        
        erro = 1 - cos;
        
        x = x_1;
        k = k+1;
    end
end

