%   O algoritmo decomp_LU decompõe a matriz A e escreve L\U sobre A
%   Este algoritmo NÃO realiza a decomposição LU, ele meramente
%   separa L de U da matriz L\U já decomposta sobre A

function [L,U] = extrair_LU_da_decomposta(A_LU)
    n = length(A_LU);
    L = zeros(n,n);
    U = zeros(n,n);
    L(1,1) = 1;
    for i = 2:n
        for j = 1:i-1
            L(i,j) = A_LU(i,j);
        end
        L(i,j+1) = 1;
    end
    for i = n:-1:1
        for j = i:n
            U(i,j) = A_LU(i,j);
        end
    end
end