function Norma = norma_infinita_Matrix(A)

[m, n] = size(A); 
Norma = 0; 

    for i = 1:m
        soma = 0;
        for j = 1:n
            soma = soma + abs(A(i,j)); 
        end
        if soma > Norma
            Norma = soma;
        end
    end
end
