
function matrix_transposta = transposta(any)
    [m, n] = size(any);
    matrix_transposta = zeros(n, m);
    for i = 1:m
        for j = 1:n
            matrix_transposta(j, i) = any(i, j);
        end
    end
end