% PROVA 1.1: CÁLCULO NUMÉRICO
% ALUNO: IGOR RODRIGUES
%FUNCÔES BASEADAS DO LIVRO FREDERICO FERREIRA 

fprintf('\n=================== QUESTÃO 1 ====================\n');

% DADOS DE ENTRADA:

A = [ 2 1 -8 1 ; 1 -1 2 -5; 10 2 1 2; -1 -9 -2 1];

b = [5,0,13,-2];

%% Letra (a)
fprintf('\n------------------ Letra (a) -----------------------\n');

% Decomposição LU com pivotação parcial e determinante
[A_LU, det, pivot] = decomp_LU(A);
[L, U] = extrair_LU_da_decomposta(A_LU);
fprintf("\n A matriz decomposta (LU):\n");disp(A_LU);
fprintf('\nMatriz triãngular inferior (L):\n'); disp(L);
fprintf('\nMatriz triãngular superior (U):\n'); disp(U)
fprintf('\nVetor de pivotação:\n'); disp(matriz_de_pivotacao(pivot));
fprintf('\nDeterminante (det): '); disp(det);

%% Letra (b)
fprintf('\n------------------ Letra (b) -----------------------\n');

% Solução dos sistemas triangulares L y = b -> U x = y
y = subst_sucessivas_pivotal(L, b, pivot);
x = subst_retroativas(U, y);

fprintf('\nVetor solução do sistema triãngular inferior (y):\n'); disp(transposta(y));
fprintf('\nVetor solução do sistema triãngular superior (x):\n'); disp(transposta(x));

%% Letra (c)
fprintf('\n------------------ Letra (c) -----------------------\n');

% Verificação se o sistema pode ser resolvido por métodos iterativos
[possivelSeDominante, matrizPermutacao] = diagonal_dominante(A);

if possivelSeDominante
    fprintf("\nO sistema pode ser resolvido por um método iterativo,pois")
    fprintf("\nexistem permutações que torna a matrix DIAGONALMENTE DOMINANTE!\n");
    %A_P =A*P
    A_permutada = matrizPermutacao * A;
    fprintf("\nMatriz A permutada:\n"); disp(A_permutada);

    % Método de Gauss-Seidel
    toler = 1e-5;
    iterMax = 30;
    [x_it, iter, ~] = gauss_seidel(A_permutada,b, toler, iterMax);
    fprintf("\nQuantidade de iterações para tolerância de %.1e: %d\n", toler, iter);
else
    fprintf("\nO sistema não pode ser resolvido por métodos iterativos diretamente,\npois a matriz não é diagonalmente dominante.\n");
end

%% Letra (d)
fprintf('\n------------------ Letra (d) -----------------------\n');

% Número de condição e inversa da matriz
k = numero_condicao(A);
fprintf("\nNúmero de condição: %.6f\n", k);

% Cálculo da inversa via decomposição LU
A_inv = A_inv_LU(A);
fprintf("\nMatriz inversa de A:\n"); disp(A_inv);

%% Letra (e)
fprintf('\n------------------ Letra (e) -----------------------\n');

% Refinamento sucessivo com arredondamento de dois dígitos

% DADOS DE ENTRADA:
criterio = 1e-4;
precisao = 2;

[x_refinado, interacoes] = refinar_LU(A, b, criterio, precisao);
fprintf("\nVetor solução após refinamento:\n"); disp(transposta(x_refinado));
fprintf("Quantidade de iterações: %d\n", interacoes);

% Vetor resíduo
r = vetor_residuo(A, b, x_refinado);
fprintf("\nVetor resíduo:\n"); disp(transposta(r));

% vetor correção
c = vetor_correcao(L, U, r, pivot, precisao);
fprintf("\nO vetor correção:\n"); disp(transposta(c));

%% Questão 2: Decomposição de Cholesky
fprintf('\n==================== QUESTÃO 2 ====================\n');

% DADOS DE ENTRADA:

B = [4 1 1; 1 5 2; 1 2 6];
c = [6 13 16];

%% Letra (a)
fprintf('\n------------------ Letra (a) -----------------------\n');

% Decomposição de Cholesky: B = L * L'
[L, det_B, erro_cond] = cholesky(B);
L_transposta = transposta(L);

fprintf("\nMatriz L:\n"); disp(L);
fprintf("\nMatriz L transposta (L'):\n"); disp(L_transposta);

% Resolvendo os sistemas triangulares
fprintf("\nResolvendo L * y = c (substituição sucessiva)\ny_chol:\n");
y_chol = subs_sucessivas(L, c); disp(transposta(y_chol));

fprintf("\nResolvendo L' * x = y (substituição retroativa):\n x_chol:\n");
x_chol = subst_retroativas(L_transposta, y_chol); disp(transposta(x_chol));

% Vetor resíduo
r_chol = vetor_residuo(B, c, x_chol);
fprintf("\nVetor resíduo (Cholesky):\n"); disp(transposta(r_chol));

%% FUNÇÕES UTILIZADAS:

%% 1
function [A, det, pivot] = decomp_LU(A)
    
    n = length(A);
    pivot = zeros(1,n);
    for i = 1:n
        pivot(i) = i; det = 1;
    end
    for j = 1:n-1
        % Escolha do elemento pivô
        p = j;
        A_max = abs(A(j,j));
        for k = j+1:n
            if abs(A(k,j)) > A_max
                A_max = abs(A(k,j)); p = k;
            end
        end
        if p ~= j
            % Troca de linhas
            for k = 1:n
                t = A(j,k); A(j,k) = A(p,k); A(p,k) = t;
            end
            m = pivot(j); pivot(j) = pivot(p); pivot(p) = m;
            det = -det;
        end
        det = det * A(j,j);
        if abs(A(j,j)) ~= 0
            % Eliminação de Gauss
            r = 1/A(j,j);
            for i = j+1:n
                mult = A(i,j) * r; A(i,j) = mult;
                for k = j+1:n
                    A(i,k) = A(i,k) - mult * A(j,k);
                end
            end
        end
    end
    det = det * A(n,n);
end

%% 2
function [L,U] = extrair_LU_da_decomposta(A_LU)
    n = length(A_LU);
    L = zeros(n,n);
    U = zeros(n,n);
    L(1,1) = 1;
    for i = 2:n
        for j = 1:i-1
            L(i,j) = A_LU(i,j);
        end
        L(i,j+1) = 1;
    end
    for i = n:-1:1
        for j = i:n
            U(i,j) = A_LU(i,j);
        end
    end
end

%% 3
function [matriz_pivot] = matriz_de_pivotacao(pivot)
    n = length(pivot);
    matriz_I = eye(n);
    % Reordena as linhas da matriz identidade de acordo com o vetor
    % pivotação
    matriz_pivot = matriz_I(pivot, :);
end


%% 4
function [y] = subst_sucessivas_pivotal(L, b, pivot)
    n = length(L);
    k = pivot(1); y(1) = b(k);
    for i = 2:n
        soma = 0;
        for j = 1:i-1
            soma = soma + L(i,j)*y(j);
        end
        k = pivot(i); y(i) = b(k) - soma;
    end
end


%% 5
function [x] = subst_retroativas(U,y)
    n = length(U);
    x = zeros(1,n); % Vetor de zeros 
    x(n) = y(n)/U(n,n);
    for i = n-1:-1:1
        soma = 0;
        for j = i+1:n
            soma = soma + U(i,j)*x(j);
        end
        x(i) = (y(i) - soma)/U(i,i);
    end
end




%% 6
function [possivelSeDominante, matrizPermutacao] = diagonal_dominante(A)
    % Verifica se A é quadrada
    [n, m] = size(A);
    if n ~= m
        error('A matriz deve ser quadrada.');
    end

    % Inicializa saída
    possivelSeDominante = false;
    matrizPermutacao = eye(n);  % identidade 

    % Testa se A já é diagonalmente dominante
    if matriz_diagonalmente_dominante(A)
        possivelSeDominante = true;
        return;
    end

    % Testa todas as permutações possíveis de linhas!

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

%% 7
function [x, iter, condErro] = gauss_seidel(A, b, toler, iterMax)
    n = length(A);
    x = zeros(1, n);
    v = zeros(1, n);
    for i = 1:n
        r = 1/A(i,i);
        for j = 1:n
            if i ~= j
                A(i,j) = A(i,j) * r;
            end
        end
        b(i) = b(i) * r; x(i) = b(i);
    end
    iter = 0;
    
    % Representação das soluções
    fprintf("\nSolução do sistema linear pelo método de Gauss-Seidel:\n" + ...
        "Iter    ")
    for i = 1:n
        if i == 2
            fprintf("x%d          ", i)
        elseif i == n
            fprintf("x%d       ", i)
        else
            fprintf("x%d         ", i)
        end
    end
    fprintf("normaRelativa\n")

    % Iterações de Gauss-Seidel
    while true
        iter = iter + 1;
        for i = 1:n
            soma = 0;
            for j = 1:n
                if i ~= j
                    soma = soma + A(i,j) * x(j);
                end
            end
            v(i) = x(i); x(i) = b(i) - soma;
        end
        normaNum = 0; normaDen = 0;
        for i = 1:n
            t = abs(x(i) - v(i));
            if t > normaNum
                normaNum = t;
            end
            if abs(x(i)) > normaDen
                normaDen = abs(x(i));
            end
        end
        normaRel = normaNum/normaDen;
        % Mostra valores
        disp([sprintf(' %d  ', iter), sprintf('  %.5f  ', x), sprintf('  %.6e', normaRel)]);

        % Teste de convergência
        if normaRel <= toler || iter >= iterMax
            break
        end
    end
    if normaRel <= toler
        condErro = 0;
    else
        condErro = 1;
    end
end
%% 8
function k = numero_condicao(A)
    % Número de Condição usando norma_1
    % A - Matriz de entrada
    
    % Calcula a norma_1 de A
    norma_A = Norma_1(A);
    
    % Calcula a norma_1 da inversa de A
    A_inv = inver(A);
    norma_A_inv = Norma_1(A_inv);
    
    % Número de Condição k = norma(A) * norma(A^-1)
    k = norma_A * norma_A_inv;
end

%% 9

function A_inv = A_inv_LU(A)
    % Função para calcular a inversa de uma matriz A usando decomposição LU
   
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

    % Resolver para cada coluna da matriz identidade (I)
    for i = 1:n
        e = zeros(n, 1);
        e(i) = 1;

        Y = L \ (P * e);
        A_inv(:, i) = U \ Y;
    end
end
%% 10
function [x, interacoes] = refinar_LU(A, b, criterio, precisao)

    % Refinamento de solução de sistema linear Ax = b usando decomposição LU
    % A - matriz do sistema
    % b - vetor do termo independente
    % criterio - critério de parada baseado na norma do vetor de correção
    % precisao - número de casas decimais usadas nas operações

    [A_LU, ~, pivot] = decomp_LU(A); 
    interacoes = 0;
    [L, U] = extrair_LU_da_decomposta(A_LU);

    % Solução inicial
    y = subst_sucessivas_pivotal(L, b, pivot);
    x = subst_retroativas(U, y);
    x = round(x, precisao);  % Aplica precisão
    fprintf("\nx0 =\n\n"); 
    disp(transposta(x))

    % Refinamento iterativo
    flag = true;
    while flag
        r = vetor_residuo(A, b, x);
        r = round(r, precisao);

        y = subst_sucessivas_pivotal(L, r, pivot);
        c = subst_retroativas(U, y);
        c = round(c, precisao);

        x = x + c;
        x = round(x, precisao);

    % Representação
        fprintf("\nr%d =\n\n", interacoes);
        disp(transposta(r))
        fprintf("\nc%d =\n\n", interacoes); 
        disp(transposta(c))
        fprintf("\nx%d =\n\n", interacoes+1); 
        disp(transposta(x))     

        interacoes = interacoes + 1;

        if Norma_Infinito(c) < criterio
            flag = false;
        end
    end
end
%% 11
function r = vetor_residuo(A, b, x)
  n = length(b);
  r = zeros(1, n);
  for i = 1:n
    soma = 0;
    for j = 1:n
      soma = soma + A(i, j) * x(j);
    end
    r(i) = b(i) - soma;
  end
end

%% 12
function matrix_transposta = transposta(any)
    [m, n] = size(any);
    matrix_transposta = zeros(n, m);
    for i = 1:m
        for j = 1:n
            matrix_transposta(j, i) = any(i, j);
        end
    end
end

%% 13
function Norma = Norma_1(X)
%Essa função calcula a norma-1 de um VETOR ou MATRIZ.
%Obs:Norma de soma máxima de coluna.

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
%% 14
function [L, det, condErro] = cholesky(A)

    n = length(A);
    
    % Verifica se A é simétrica
    if ~eh_simetrica(A)
        error('A matriz não é simétrica. O método de Cholesky exige simetria.');
    end
    
    condErro = 0; det = 1;
    for j = 1:n
        soma = 0;
        for k = 1:j-1
            soma = soma + A(j,k)^2;
        end
        t = A(j,j) - soma;
        if t > 0
            A(j,j) = sqrt(t); r = 1/A(j,j); det = det*t;
        else
            condErro = 1;
            disp("A matriz não é definida positiva")
            return
        end
        for i = j+1:n
            soma = 0;
            for k = 1:j-1
                soma = soma + A(i,k)*A(j,k);
            end
            A(i,j) = (A(i,j) - soma)*r;
        end
    end
    
    for i = 1:n %Zera toda a parte não inferior da matriz
        for j = i+1:n
            A(i,j) = 0;
        end
    end
    L = A; %Copia a matriz decomposta para L
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

%% 15
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

%% 16
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
%% 17
function c = vetor_correcao(L, U, r, pivot, precisao)

    y = subst_sucessivas_pivotal(L, r, pivot);
    c = subst_retroativas(U, y);
    c = round(c, precisao);
end

%% FIM
