function A_inversa = inver(A)
    % Verifica se a matriz eh quadrada
    [m, n] = size(A);
    if m ~= n
        error('A matriz nao eh quadrada. A inversao nao eh possivel.');
    end
    
    % Verifica se o determinante eh zero
    if det(A) == 0
        error('A matriz tem determinante igual a zero, nao eh possivel inverter.');
    end
    
    % Numero de linhas ou colunas (A eh quadrada)
    n = length(A);
    
    % Cria a matriz aumentada [A | I], onde I eh a identidade
    matriz_aumentada = [A eye(n)];
    
    % Processo de Eliminacao de Gauss-Jordan
    for col = 1:n
        % Pivo eh o elemento diagonal
        pivo = matriz_aumentada(col, col);
        
        % Normaliza a linha do pivo
        matriz_aumentada(col, :) = matriz_aumentada(col, :) / pivo;
        
        % Elimina os elementos acima e abaixo do pivo
        for linha = 1:n
            if linha ~= col
                fator = matriz_aumentada(linha, col);
                matriz_aumentada(linha, :) = matriz_aumentada(linha, :) - fator * matriz_aumentada(col, :);
            end
        end
    end
    
    % A inversa de A esta na parte direita da matriz aumentada
    A_inversa = matriz_aumentada(:, n+1:end);
end
