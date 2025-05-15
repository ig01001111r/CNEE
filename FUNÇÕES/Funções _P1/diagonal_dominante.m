function [possivelSeDominante, matrizPermutacao] = diagonal_dominante(A)
    % Verifica se A é quadrada
    [n, m] = size(A);
    if n ~= m
        error('A matriz deve ser quadrada.');
    end

    % Inicializa saída
    possivelSeDominante = false;
    matrizPermutacao = eye(n);  % identidade por padrão

    % Testa se A já é diagonalmente dominante
    if matriz_diagonalmente_dominante(A)
        possivelSeDominante = true;
        return;
    end

    % Testa todas as permutações possíveis de linhas
    permsLinhas = perms(1:n);  % todas as permutações de linhas
    for k = 1:size(permsLinhas, 1)
        P = eye(n);
        P = P(permsLinhas(k, :), :);  % gera matriz de permutação
        A_permutada = P * A;

        if matriz_diagonalmente_dominante(A_permutada)
            possivelSeDominante = true;
            matrizPermutacao = P;
            return;
        end
    end
end

% Função auxiliar: verifica se matriz é diagonalmente dominante
function ehDominante = matriz_diagonalmente_dominante(M)
    n = size(M, 1);
    ehDominante = true;
    for i = 1:n
        diagonal = abs(M(i, i));
        somaLinha = sum(abs(M(i, :))) - diagonal;
        if diagonal < somaLinha
            ehDominante = false;
            break;
        end
    end
end
