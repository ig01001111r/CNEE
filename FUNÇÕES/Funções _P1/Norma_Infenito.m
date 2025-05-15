function Norma = Norma_Infinito(c)
    
% Norma_Infinito - Calcula a norma infinito (máximo valor absoluto) de um vetor.
% Entrada:
%   c - vetor numérico
% Saída:
%   Norma - maior valor absoluto dos elementos do vetor

    n = length(c);
    Norma = abs(c(1));
    
    for i = 2:n
        if abs(c(i)) > Norma
            Norma = abs(c(i));
        end
    end
end
