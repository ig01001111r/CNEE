% PARTE 2 DA PROVA 2 DE CÁLCULO NUMÉRICO
% ALUNO: IGOR RODRIGUES 

minimo1 = [0 0 0 0]; % mínimo local inicial (4 variáveis)
tolerancia1 = 0.05;
iterMAX1 = 100;

fprintf("\nLetra (a):\n")

[aproximacao, ~, ~] = metGradienteDescendente([], minimo1, tolerancia1, iterMAX1);

fprintf("\n\nA solução inicial, ao sair do método gradiente, é:\n");
disp(aproximacao);

fprintf("\nLetra (b)\n")

tolerancia2 = 10^(-6);
iterMAX2 = 100;

%{

syms x y z w
simbolos = [x, y, z, w];

f1 = sin(x) + cos(y) + z - 2 + w;
f2 = cos(x) - sin(y) - w + z;
f3 = z^2 + w^2 - 2;
f4 = x + y - pi/2;

%}

syms x1 x2 x3 x4
simbolos = [x1, x2, x3, x4];

f1 = sin(x1) + cos(x2) + x3 - 2 + x4;
f2 = cos(x1) - sin(x2) - x4 + x3;
f3 = x3^2 + x4^2 - 2;
f4 = x1 + x2 - pi/2;

vet_func = [f1; f2; f3; f4];

[solucao, ~] = metNewtonRaphsonN(vet_func, simbolos, 4, aproximacao, tolerancia2, iterMAX2);

fprintf("Usando o método de Newton-Raphson para refinar a solução temos:\n");
disp(solucao);


% Essa função implementa o método do gradiente descendente para minimizar a função escalar
% G = f1^2 + f2^2 + f3^2 + f4^2,
% onde cada fi é uma função não linear de 4 variáveis.
% O objetivo é encontrar valores das variáveis que aproximem a solução do sistema fi = 0.

function [aproximacao,gNew,k] = metGradienteDescendente(~, minimo, tolerancia, iterMAX)
  
    syms f1(x1,x2,x3,x4) f2(x1,x2,x3,x4) f3(x1,x2,x3,x4) f4(x1,x2,x3,x4)

    f1(x1,x2,x3,x4) = sin(x1) + cos(x2) + x3 - 2 + x4;
    f2(x1,x2,x3,x4) = cos(x1) - sin(x2) - x4 + x3;
    f3(x1,x2,x3,x4) = x3^2 + x4^2 - 2;
    f4(x1,x2,x3,x4) = x1 + x2 - pi/2;
 
    % O gradiente da função G(x) = f1^2 + f2^2 + f3^2 + f4^2 é dado por:
    % ∇G(x) = 2*f1*∇f1 + 2*f2*∇f2 + 2*f3*∇f3 + 2*f4*∇f4

    syms gs(x1,x2,x3,x4)
    gs(x1,x2,x3,x4) = f1^2 + f2^2 + f3^2 + f4^2;
    %gradiente = gradient(gs, [x1 x2 x3 x4]); % vetor gradiente de G.
    
    % deriva parcialmente 
    gradiente(x1,x2,x3,x4) = [
        2*f1*diff(f1,x1) + 2*f2*diff(f2,x1) + 2*f3*diff(f3,x1) + 2*f4*diff(f4,x1),
        2*f1*diff(f1,x2) + 2*f2*diff(f2,x2) + 2*f3*diff(f3,x2) + 2*f4*diff(f4,x2),
        2*f1*diff(f1,x3) + 2*f2*diff(f2,x3) + 2*f3*diff(f3,x3) + 2*f4*diff(f4,x3),
        2*f1*diff(f1,x4) + 2*f2*diff(f2,x4) + 2*f3*diff(f3,x4) + 2*f4*diff(f4,x4)
    ];

    aproximacao = minimo(:);
    k = 1;

    while k <= iterMAX
        g1 = double(gs(aproximacao(1),aproximacao(2),aproximacao(3),aproximacao(4)));
        z = double(gradiente(aproximacao(1),aproximacao(2),aproximacao(3),aproximacao(4)));
        z0 = norm(z,2);

        if z0 == 0
            fprintf("\nGradiente zero. Procedimento concluído!\n");
            break;
        end

        z = z/z0;
        alpha3 = 1;
        aux3 = aproximacao - alpha3*z(:);
        g3 = double(gs(aux3(1),aux3(2),aux3(3),aux3(4)));

        while g3 >= g1
            alpha3 = alpha3/2;
            aux3 = aproximacao - alpha3*z(:);
            g3 = double(gs(aux3(1),aux3(2),aux3(3),aux3(4)));

            if alpha3 < tolerancia/2
                fprintf("\nNenhuma melhora significativa.\n");
                return;
            end
        end

        alpha2 = alpha3/2;
        aux2 = aproximacao - alpha2*z(:);
        g2 = double(gs(aux2(1),aux2(2),aux2(3),aux2(4)));

        h1 = (g2 - g1)/alpha2;
        h2 = (g3 - g2)/(alpha3 - alpha2);
        h3 = (h2 - h1)/alpha3;

        alpha0 = 0.5*(alpha2 - h1/h3);
        aux0 = aproximacao - alpha0*z(:);
        g0 = double(gs(aux0(1),aux0(2),aux0(3),aux0(4)));

        if g0 <= g3 && g0 <= g2
            alpha = alpha0;
        elseif g3 <= g2
            alpha = alpha3;
        else
            alpha = alpha2;
        end

        aproximacao = aproximacao - alpha*z(:);

        fprintf("Iteração %d:\n", k);
        fprintf("Aproximação: %.6f %.6f %.6f %.6f\n", aproximacao);
        fprintf("G(aprox): %.10f\n\n", double(gs(aproximacao(1),aproximacao(2),aproximacao(3),aproximacao(4))));

        gNew = double(gs(aproximacao(1),aproximacao(2),aproximacao(3),aproximacao(4)));

        if abs(gNew - g1) < tolerancia
            fprintf("Convergência alcançada.\n");
            return;
        end

        k = k + 1;
    end

    if k > iterMAX
        fprintf("Número máximo de iterações excedido!\n");
    end
end

% ---------------- FUNÇÃO: Newton-Raphson ----------------
function [x, Iter] = metNewtonRaphsonN(vet_func, simbolos, n, solucao, tolerancia, maxIter)
    Iter = 0;

    while true
        Iter = Iter + 1;
        if Iter > maxIter
            fprintf("Número máximo de iterações excedido!\n");
            break
        end

        vetN = zeros(n,1);
        for i = 1:n
            vetN(i) = double(subs(vet_func(i), simbolos, num2cell(solucao(:)')));
        end

        A = zeros(n,n);
        for i = 1:n
            for j = 1:n
                A(i,j) = double(subs(diff(vet_func(i), simbolos(j)), simbolos, num2cell(solucao(:)')));
            end
        end

        vet_b = -vetN;

        [A_LU, ~, pivot] = decomp_LU(A);
        [L, U] = extrair_LU_da_decomposta(A_LU);
        Y = subst_sucessivas_pivotal(L, vet_b, pivot);
        incremento = subst_retroativas(U, Y);

        solucao = solucao + incremento';

        if norm(incremento, inf) <= tolerancia
            break
        end
    end
    x = solucao;
end

% ---------------- AUXILIARES LU E SUBSTITUIÇÃO ----------------
function [A, detA, pivot] = decomp_LU(A)
    n = length(A);
    pivot = 1:n;
    detA = 1;
    for j = 1:n-1
        [~, p] = max(abs(A(j:n,j)));
        p = p + j - 1;
        if p ~= j
            A([j p], :) = A([p j], :);
            pivot([j p]) = pivot([p j]);
            detA = -detA;
        end
        detA = detA * A(j,j);
        if abs(A(j,j)) > 0
            for i = j+1:n
                A(i,j) = A(i,j) / A(j,j);
                A(i,j+1:n) = A(i,j+1:n) - A(i,j) * A(j,j+1:n);
            end
        end
    end
    detA = detA * A(n,n);
end

function [L,U] = extrair_LU_da_decomposta(A_LU)
    n = length(A_LU);
    L = eye(n);
    U = zeros(n);
    for i = 1:n
        for j = 1:i-1
            L(i,j) = A_LU(i,j);
        end
        for j = i:n
            U(i,j) = A_LU(i,j);
        end
    end
end

function y = subst_sucessivas_pivotal(L, b, pivot)
    n = length(L);
    y = zeros(n,1);
    y(1) = b(pivot(1));
    for i = 2:n
        soma = 0;
        for j = 1:i-1
            soma = soma + L(i,j)*y(j);
        end
        y(i) = b(pivot(i)) - soma;
    end
end

function x = subst_retroativas(U,y)
    n = length(U);
    x = zeros(1,n);
    x(n) = y(n)/U(n,n);
    for i = n-1:-1:1
        soma = U(i,i+1:n)*x(i+1:n)';
        x(i) = (y(i) - soma)/U(i,i);
    end
end
