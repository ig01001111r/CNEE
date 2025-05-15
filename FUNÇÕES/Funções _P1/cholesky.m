%Obj: Fazer a decomposição LL^(t) de uma matriz (A)
% simétrica e definada positiva
% Parâmetros de entrada:A {matriz a ser decomposta}
% Parâmetros de saida : L,det,codErro {fator L escrito sobre A, determinante, e condição de erro}
% Ax=B => (L*L')x=B => L'x=y e L*y=B

function [L, det, condErro] = cholesky(A)

    n = length(A);
    
    % Verifica se A é simétrica
    if ~eh_simetrica(A)
        error('A matriz não é simétrica. O método de Cholesky exige simetria.');
    end
    
    condErro = 0; det = 1;
    for j = 1:n
        soma = 0;
        for k = 1:j-1
            soma = soma + A(j,k)^2;
        end
        t = A(j,j) - soma;
        if t > 0
            A(j,j) = sqrt(t); r = 1/A(j,j); det = det*t;
        else
            condErro = 1;
            disp("A matriz não é definida positiva")
            return
        end
        for i = j+1:n
            soma = 0;
            for k = 1:j-1
                soma = soma + A(i,k)*A(j,k);
            end
            A(i,j) = (A(i,j) - soma)*r;
        end
    end
    
    for i = 1:n %Zera toda a parte não inferior da matriz
        for j = i+1:n
            A(i,j) = 0;
        end
    end
    L = A; %Copia a matriz decomposta para L
end

function simetrica = eh_simetrica(M)
    [l, c] = size(M);
    simetrica = true;
    for i = 1:l
        for j = 1:c
            if M(i,j) ~= M(j,i)
                simetrica = false;
                return;
            end
        end
    end
end