function Norma = Norma_1(X)
%Essa função calcula a norma-1 de um VETOR ou MATRIZ.

    [m,n] = size(X);
    Norma = 0;
    for j=1:n
        soma = 0;
        for i=1:m
            soma = soma + abs(X(i,j));
        end
        if soma>Norma
            Norma = soma;
        end
    end
end