function f = interNewton(m, x, y, z)
    
    % newton_interp: interpola o valor f(z) usando o polinômio de Newton
    % m: número de pontos
    % x: vetor com as abscissas (tamanho m)
    % y: vetor com as ordenadas (tamanho m)
    % z: ponto no qual se deseja interpolar
    
    Dely = y;  % Inicializa com os y dados

    % Construção das diferenças divididas
    for k = 1 : m - 1
        for i = m : -1 : k + 1
            Dely(i) = (Dely(i) - Dely(i - 1)) / (x(i) - x(i - k));
        end
    end

    % Avaliação usando Horner
    r = Dely(m);
    for i = m - 1 : -1 : 1
        r = r * (z - x(i)) + Dely(i);
    end

    f = r;
end
