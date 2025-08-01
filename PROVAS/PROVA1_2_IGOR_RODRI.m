% PROVA 1.1: CÁLCULO NUMÉRICO
% ALUNO: IGOR RODRIGUES
% FUNÇÕES BASEADAS NO LIVRO DE FREDERICO FERREIRA

%% QUESTÃO 1

fprintf('QUESTÃO 1:\n');

a = [0.1 0.15 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.5 2 3 4 5 10 20 50 100 200 500 1000]'; % coeficiente de penetração 
N = [1.1 0.98 0.9 0.88 0.87 0.84 0.83 0.8 0.79 0.79 0.78 0.73 0.72 0.69 0.65 0.61 0.55 0.48 0.4 0.36 0.3 0.27 0.23]'; % N para b = 0.1

n = length(a);

%% Regressão polinomial de grau 7 (8 parâmetros)
X = [a];
[b, R2, sigma2] = regressao_Multipla(n, 1, 8, X, N);

%% Regressão Exponencial: ln(N) = ln(a) - b*a
lnN = log(N);      % logaritmo natural
ln_a = log(a);
[b_exp, R2_exp, sigma2_exp] = regressao_Multipla(n, 1, 2, ln_a, lnN);

fprintf("\nParâmetros do modelo exponencial:\n");
disp(b_exp);

%% Comparação dos modelos
fprintf('+-------------------------------------------------------+\n');
fprintf('|       MODELOS         |   R²          |   σ²          |\n');
fprintf('+-------------------------------------------------------+\n');
fprintf('|  Melhor Modelo        | R² = %.6f | σ² = %.6f |\n', R2, sigma2);
fprintf('| Modelo Exponencial    | R² = %.6f | σ² = %.6f |\n', R2_exp, sigma2_exp);
fprintf('+--------------------------------------------------------+\n');

%% Estimativa para α = 0.337
a_est = 0.337;

M_b = b(1) + b(2)*a_est + b(3)*a_est^2 + b(4)*a_est^3 + b(5)*a_est^4 + b(6)*a_est^5 + b(7)*a_est^6 + b(8)*a_est^7;
Ln_N = b_exp(1) + b_exp(2)*log(a_est);
N_a = exp(Ln_N);

fprintf('\n--- ESTIMATIVAS PARA a = %.3f ---\n', a_est);
fprintf('Melhor Modelo Polinomial:      N = %.6f\n', M_b);
fprintf('Modelo Exponencial:            N = %.6f\n', N_a);

%% Justificativa da escolha do melhor modelo
fprintf("\n--- JUSTIFICATIVA DA ESCOLHA DO MODELO ---\n");
fprintf("O melhor modelo é o polinomial de grau 7 pois possui R² = %.4f mais próximo de 1 e σ² = %.6f menor.\n", R2, sigma2);

%% QUESTÃO 2

fprintf('\nQUESTÃO 2:\n');

% Interpolação de Lagrange para N4(0.713)
A_lag = [0.5 0.6 0.7 0.8 0.9];
N_lag = [0.87 0.84 0.83 0.80 0.79];
N4_0_713 = interpLagrange(5, A_lag, N_lag, 0.713);

% Interpolação de Gregory-Newton para N2(0.275)
A_gn = [0.15 0.2 0.3];
N_gn = [0.98 0.95 0.9];
N2_0_275 = interpGregoryNewton(3, A_gn, N_gn, 0.275);

fprintf('Interpolação de Lagrange N4(0.713): %.5f\n', N4_0_713);
fprintf('Interpolação de Gregory-Newton N2(0.275): %.5f\n', N2_0_275);

%% FUNÇÕES AUXILIARES

function [b, r2, sigma2] = regressao_Multipla(n, v, p, x, y)
    if (v > 1) && (v+1 ~= p)
        error("Modelo inválido");
    end
    x = [ones(n,1), x]; % adiciona termo constante b1
    if (v == 1) && (p > 2)
        for j = 2:p-1
            x(:,j+1) = x(:,2).^j;
        end
    end
    Sxx = x' * x;
    Sxy = x' * y;
    L = regresCholesky(p, Sxx);
    t = regresSubstSucessiva(p, L, Sxy);
    U = L';
    b = regresSubstRegresiva(U, t);
    u = x * b;
    d = y - u;
    D = sum(d.^2);
    Sy2 = sum(y.^2);
    r2 = 1 - D / (Sy2 - (sum(y)^2) / n);
    sigma2 = D / (n - p);
end

function L = regresCholesky(p, A)
    L = zeros(p);
    for j = 1:p
        soma = sum(L(j,1:j-1).^2);
        t = A(j,j) - soma;
        if t <= 0
            error("Matriz não definida positiva");
        end
        L(j,j) = sqrt(t);
        for i = j+1:p
            soma = sum(L(i,1:j-1).*L(j,1:j-1));
            L(i,j) = (A(i,j) - soma) / L(j,j);
        end
    end
end

function Y = regresSubstSucessiva(p, L, b)
    Y = zeros(p,1);
    for i = 1:p
        Y(i) = (b(i) - L(i,1:i-1)*Y(1:i-1)) / L(i,i);
    end
end

function X = regresSubstRegresiva(U, y)
    p = length(y);
    X = zeros(p,1);
    for i = p:-1:1
        X(i) = (y(i) - U(i,i+1:end)*X(i+1:end)) / U(i,i);
    end
end

function interp = interpLagrange(n, X, Y, x)
    interp = 0;
    for i = 1:n
        num = 1;
        den = 1;
        for j = 1:n
            if i ~= j
                num = num * (x - X(j));
                den = den * (X(i) - X(j));
            end
        end
        interp = interp + Y(i) * (num / den);
    end
end

function interp = interpGregoryNewton(n, X, Y, x)
    delY = Y;
    for k = 1:n-1
        for i = n:-1:k+1
            delY(i) = delY(i) - delY(i-1);
        end
    end
    h = X(2) - X(1);
    u = (x - X(1)) / h;
    interp = delY(n);
    for i = n-1:-1:1
        interp = interp * (u - (i - 1)) / i + delY(i);
    end
end
