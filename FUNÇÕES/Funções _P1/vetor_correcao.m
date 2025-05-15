function c = vetor_correcao(L, U, r, pivot, precisao)
    % Calcula o vetor de correção no refinamento sucessivo
    % L, U   - matrizes da decomposição LU
    % r      - vetor resíduo (b - A*x)
    % pivot  - vetor de permutação (índices)
    % precisao - número de casas decimais para arredondar

    % Resolver L * y = P * r (substituição progressiva)
    y = subst_sucessivas_pivotal(L, r, pivot);
    
    % Resolver U * c = y (substituição regressiva)
    c = subst_retroativas(U, y);
    
    % Arredondar a correção com a precisão especificada
    c = round(c, precisao);
end
