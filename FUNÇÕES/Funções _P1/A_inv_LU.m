function A_inv = A_inv_LU(A)
    % Função para calcular a inversa de uma matriz A usando decomposição LU
    % A: Matriz quadrada de entrada
    % A_inv: Matriz inversa de A

    % Verificar se a matriz A é quadrada
    [m, n] = size(A);
    if m ~= n
        error('A matriz A deve ser quadrada!');
    end

    % Decomposição LU com matriz de permutação P
    [A_LU, P] = decomp_LU(A);

    % Extrair L e U da decomposição LU
    [L, U] = extrair_LU_da_decomposta(A_LU);

    % Inicializar a matriz inversa
    A_inv = zeros(n);

    % Resolver para cada coluna da matriz identidade I
    for i = 1:n
        e = zeros(n, 1);
        e(i) = 1;

        Y = L \ (P * e);
        A_inv(:, i) = U \ Y;
    end
end
