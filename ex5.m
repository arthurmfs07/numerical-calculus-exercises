A = input("");
sass(A);

function sass(A)
n = length(A);
beta = zeros(n, 1);

for i = 1:n
    b = 0;
    for j = 1:n
        if (i ~= j && i == 0) || i < j
            b = A(i,j) + b;
        end
        if i ~= j && i ~= 0
            b = b + A(i,j) * beta(j);
        end
    end
    b = b / A(i, i);
    beta(i) = b;
end

if max(beta) < 1
    disp("SATISFAZ");
else
    disp("NÃO SATISFAZ");
end
end
