A = input("");
b = input("");
x_0 = input("");
tol = input("");

[x,k] = gj(A,b,x_0,tol);

disp("k =");
disp(k);
disp("x =");
disp(x);

function [x, k] = gj(A, b, x_0, tol)
    n = size(A, 1); % number of lines in A
    D = diag(diag(A));
    C = eye(n) - D\A;
    g = D\b;
    kmax = 1000;
    k = 0;
    erro = Inf;
    
    while (erro > tol && k < kmax)
        k = k + 1;
        x_old = x_0;
        x_0 = C * x_0 + g;
        x = x_0;
        erro = norm(x - x_old) / norm(x); 
    end

    if (k == kmax)
        disp('Erro: método não converge');
        return;
    end
end