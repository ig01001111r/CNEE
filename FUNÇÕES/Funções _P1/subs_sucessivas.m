function [y] = subs_sucessivas(L,c) 

% Obj: Resolver o sistema triangular inferior Lx = c pelas subs. sucessivas
% Parâmetros de entrada L,c { matriz triangular infeiror e vetor independente}
% Parâmetros de saida [y] { solução do sistema triangular infeior }

n = length(L);
y = zeros(1,n);
y(1) = c(1)/L(1,1);
for i=2:n
    soma=0;
    for j=1:i-1
        soma = soma + L(i,j)*y(j);
    end
    y(i) = (c(i) - soma)/L(i,i);
end