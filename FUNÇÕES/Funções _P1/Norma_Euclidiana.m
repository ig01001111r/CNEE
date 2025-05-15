function Norma = Norma_Euclidiana(v)
% Obj: Calcular a norma-2 de um vetor

% Parâmetros de entrada v
v = [1,2,3,4];
n = length(v); %Encontra o número de elementos do vetor.
soma = 0;
for i = 1:n
    soma = soma+abs(v(i))^2;
end
% Parãmetro de saida norma-2
Norma = sqrt(soma);
 
 
 
