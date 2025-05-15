function A_inv = inv_cholesky(A)
%INVERSA_CHOLESKY Inverte uma matriz simétrica e definida positiva usando a fatoração de Cholesky
%   Entrada:
%     A - matriz simétrica e definida positiva
%   Saída:
%     A_inv - inversa da matriz A

    % Verifica se A é quadrada
    [l, c] = size(A);
    if l ~= c
        error('A matriz deve ser quadrada.');
    end

    % Verifica se A é simétrica
    if ~eh_simetrica(A)
        error('A matriz não é simétrica. O método de Cholesky exige simetria.');
    end

    n = l;  % Tamanho da matriz
    L = zeros(n, n);  % Matriz triangular inferior

    % Fatoração de Cholesky: A = L * L'
    for i = 1:n
        for j = 1:i
            soma = A(i,j);
            for k = 1:j-1
                soma = soma - L(i,k) * L(j,k);
            end
            if i == j
                if soma <= 0
                    error('A matriz não é definida positiva. Cholesky não é aplicável.');
                end
                L(i,j) = sqrt(soma);
            else
                L(i,j) = soma / L(j,j);
            end
        end
    end
 
    % Resolvendo L * Y = I (substituição sucessiva)
    Y = zeros(n);
    I = eye(n);
    for i = 1:n
        Y(:,i) = subs_sucessivas(L, I(:,i));
    end

    % Resolvendo L' * A_inv = Y (substituição retroativa)
    A_inv = zeros(n);
    for i = 1:n
        A_inv(:,i) = subs_retroativas(L', Y(:,i));
    end
end

function x = subs_sucessivas(L, b)
    n = length(b);
    x = zeros(n,1);
    for i = 1:n
        x(i) = (b(i) - L(i,1:i-1)*x(1:i-1)) / L(i,i);
    end
end

function x = subs_retroativas(U, b)
    n = length(b);
    x = zeros(n,1);
    for i = n:-1:1
        x(i) = (b(i) - U(i,i+1:end)*x(i+1:end)) / U(i,i);
    end
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
