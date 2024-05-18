A = input("");
b = input("");
w = input("");
tol = input("");
max_iter = input("");

[phi, k, res] = SOR(A, b, w, tol, max_iter);

disp(phi);
disp(k);
disp(res);  % Display all error residue values until convergence

function [phi, k, res] = SOR(A, b, w, tol, max_iter)
    % Input:
    % A: coefficient matrix
    % b: right-hand side vector
    % w: relaxation parameter
    % tol: tolerance for convergence
    % max_iter: maximum number of iterations allowed

    % Initializations
    n = size(A, 1);
    phi = zeros(n, 1);  % initial guess
    res = zeros(max_iter, 1);  % error residue vector
    k = 0;  % iteration counter

    % Main loop
    while true
        k = k + 1;
        for i = 1:n
            sigma = 0;
            for j = 1:n
                if j ~= i
                    sigma = sigma + A(i, j) * phi(j);
                end
            end
            phi(i) = (1 - w) * phi(i) + (w / A(i, i)) * (b(i) - sigma);
        end

        % Calculate error residue
        res(k) = norm(A * phi - b);

        % Check convergence
        if res(k) < tol || k >= max_iter
            break;
        end
    end
    
    % Resize the res vector to contain only the non-zero error residue values
    res = res(1:k);
end
