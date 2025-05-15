% Obs: Resolver o sistema triangular inferior Ly = Pb pela substituicões
% sucessivas, com a matrix L obtida de decomposição LU  com pivotação
% parcial
% Parâmetros de entrada L,b,pivot { matriz triangular inferior unitária,vetor independente 
% e posição dos pivôs}
% Parâmetros de saida [y] { solução do sistema triangular inferior }

function [y] = Subs_Sucessivas_pivotal(L,b,Pivotal)
n = length(L);
k = Pivotal(1); y(1) = b(k);
    for i=2:n
        Soma = 0;
        for j=1:i-1
            Soma = Soma + L(i ,j)*y(j);
        end
        k = Pivotal(i); y(i)= b(k) - Soma;
    end
end