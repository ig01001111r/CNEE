%   Resolve sistema triangular superior Ux = y pelas substituições retroativas
%   U = matriz_triangular_superior
%   x = vetor_solução
%   d = vetor_independente

function [x] = subst_retroativas(U,y)
    n = length(U);
    x = zeros(1,n); % Aloca um vetor de zeros 
    x(n) = y(n)/U(n,n);
    for i = n-1:-1:1
        soma = 0;
        for j = i+1:n
            soma = soma + U(i,j)*x(j);
        end
        x(i) = (y(i) - soma)/U(i,i);
    end
end