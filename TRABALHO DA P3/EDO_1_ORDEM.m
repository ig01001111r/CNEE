% QUESTÃO: LEI DE RESFRIAMENTO DE NEWTON
%
% De acordo com a Lei de Resfriamento de Newton, a taxa de variação da
% temperatura T(t) de um corpo é proporcional à diferença entre a sua
% temperatura e a temperatura ambiente constante Ta.
%
% A equação diferencial que modela esse fenômeno é dada por:
% dT/dt = -k*(T - Ta)
%
% Onde:
% - T(t) é a temperatura do corpo no tempo t (em °C)
% - Ta é a temperatura do meio ambiente (em °C)
% - k é uma constante positiva de resfriamento (1/min)
% - t é o tempo em minutos

% Considere as seguintes condições:
% - A temperatura inicial do corpo é T(0) = 100 °C
% - A temperatura do ambiente é Ta = 20 °C
% - Após 10 minutos, a temperatura do corpo é T(10) = 60 °C
%
% Com base nesses dados:
% (a) Encontre a solução da equação diferencial que modela T(t)
% (b) Determine a constante de resfriamento k
% (c) Esboce o gráfico da temperatura T(t) ao longo do tempo (0 ≤ t ≤ 30)

% Lei de Resfriamento de Newton
% T(t) = 80 * exp(-0.0693 * t) + 20

%% TRABALHO DE CÁLCULO NUMÉRICO 2025.1
% ALUNO: IGOR RODRIGUES 

clear; clc;

% Dados:

Ta = 20;      % Temperatura ambiente (°C)
T0 = 100;     % Temperatura inicial do corpo (°C)
T10 = 60;     % Temperatura após 10 minutos (°C)
t10 = 10;     % Tempo (minutos)

fprintf("Considerando a Lei de Resfriamento de Newton:\n");
fprintf("Temperatura inicial T(0) = %d °C\n", T0);
fprintf("Temperatura ambiente Ta = %d °C\n", Ta);
fprintf("Temperatura T(10) = %d °C após 10 minutos\n\n", T10);

% Cálculo da constante k
k_val = - (1/t10) * log( (T10 - Ta) / (T0 - Ta) );
fprintf("Constante de resfriamento k calculada: %.4f (1/min)\n\n", k_val);

% Definir função f(x,y) = dy/dt = -k*(y - Ta)
f = @(t, T) -k_val * (T - Ta);

% Parâmetros do método numérico
a = 0;        % tempo inicial
b = 30;       % tempo final
m = 100;      % número de passos
y0 = T0;      % condição inicial

% Resolver com método de Euler usando sua função
[VetX, VetY] = Euler(a, b, m, y0, f);

% Solução analítica para comparação
T_analitica = Ta + (T0 - Ta)*exp(-k_val*VetX);

% Mostrar comparação
fprintf("i \t t \t\t Euler \t\t Analítica\n");
for i = 1:length(VetX)
    fprintf("%d \t %.4f \t %.4f \t %.4f\n", i-1, VetX(i), VetY(i), T_analitica(i));
end

% Gráfico 
figure
hold on
plot(VetX, T_analitica, 'r--', 'LineWidth', 2) % solução analítica
plot(VetX, VetY, 'bo-')                         % método de Euler
xlabel('Tempo (min)');
ylabel('Temperatura (°C)');
legend('Solução Analítica', 'Método de Euler', 'Location', 'northeast');
title('Lei de Resfriamento de Newton - Método de Euler');
grid on
hold off


function [VetX, VetY] = Euler(a, b, m, y0, f)

    %{ 
       Entradas:
       a  - limite inferior do intervalo 
       b  - limite superior do intervalo 
       m  - número de subintervalos 
       y0 - valor inicial da função y no ponto a, ou seja, y(a) = y0
       f  - função anônima que representa a EDO, f(x,y) = dy/dx
    
       Saídas:
       VetX - vetor contendo os valores de x nos pontos de discretização
       VetY - vetor contendo a solução aproximada y nos respectivos pontos VetX
    %}

    h = (b - a) / m;          % tamanho do passo
    VetX = zeros(1, m+1);     % inicializa vetor de x
    VetY = zeros(1, m+1);     % inicializa vetor de y
    VetX(1) = a;              % ponto inicial x0 = a
    VetY(1) = y0;             % valor inicial y0
    
    for i = 1:m
        x = VetX(i);
        y = VetY(i);
        Fxy = f(x, y);        % calcula f(x_i, y_i)
        
        % Método de Euler:
        y_next = y + h * Fxy;
        x_next = x + h;
        
        VetX(i+1) = x_next;
        VetY(i+1) = y_next;
        
        % Mostrar valores a cada passo 
        fprintf('i = %d, x = %.4f, y = %.4f, f(x,y) = %.4f\n', i, x_next, y_next, Fxy);
    end
end
