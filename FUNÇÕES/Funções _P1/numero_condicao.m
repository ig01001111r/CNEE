%MAL-DIMENSINAMENTO
%Interpretação dos valores de k:
% k∼1: A matriz está bem condicionada. A solução do sistema é estável e não sensível a pequenas variações nas entradas de 
% k≫1: A matriz é mal condicionada. Pequenas alterações em 
% b podem causar grandes variações na solução, tornando a solução numérica instável.
% k muito grande (por exemplo, k∼10^(6)ou mais): O sistema é extremamente mal condicionado e a solução será muito imprecisa.
function K = numero_condicao(A)
    % Número de Condição usando norma 1
    % A - Matriz de entrada
    
    % Calcula a norma infinita de A
    norma_A =  Norma_1(A);
    
    % Calcula a norma infinita da inversa de A
    A_inv = inver(A);
    norma_A_inv =   Norma_1(A_inv);
    
    % Número de Condição K = norma(A) * norma(A^-1)
    K = norma_A * norma_A_inv;
end

