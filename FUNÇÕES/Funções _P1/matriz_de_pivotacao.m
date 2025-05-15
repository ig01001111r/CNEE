%   Calcula a matriz pivot a partir do vetor pivotação 

function [matriz_pivot] = matriz_de_pivotacao(pivot)
    n = length(pivot);
    matriz_I = eye(n);
    % Reordena as linhas da matriz identidade de acordo com o vetor pivot
    matriz_pivot = matriz_I(pivot, :);
end