A = input("");
tol = input("");

[V,D,k,erros] = francis(A,tol);
[Q,R] = mfgs(A);

disp('D = ');
disp(sort(D, 'descend')); % a variável D deve ser um vetor coluna que contém os autovalores em ordem decrescente
disp('k = ');
disp(k);
disp('erros = ');
disp(erros); % a variável erros deve ser um vetor coluna

%if ~isequal(A, A') || any(eig(A) <= 0)

if ~isequal(A, A') || any(eig(A) <= 0) %verificar se A SPD
    r = 'não é';
else
    r = 'é';
end

disp(r);

function [Q,R] = mfgs(A)
[m,n] = size(A);
V = A;
Q = zeros(m,n);
R = zeros(n,n);

for i=1:n
    R(i,i) = norm(V(:,i));
    Q(:,i) = V(:,i)/R(i,i);
    for j=i+1:n
        R(i,j) = Q(:,i)' * V(:, j);
        V(:, j) = V(:, j) - R(i, j) * Q(:, i);
    end
end
end

function [V, D, k, erros] = francis(A, tol)
    n = size(A, 1);
    V = eye(n);
    erro = inf;
    k = 0;
    erros = []; % Initialize as an empty column vector
    
    while erro > tol
        k = k + 1;
        [Q, R] = mfgs(A);
        A = R * Q;
        V = V * Q;
        erro = sqrt(norm(A, 'fro')^2 - sum(diag(A).^2));
        erros(end+1, 1) = erro; % Concatenate errors as a column vector
    end
    
    D = diag(A);
end

