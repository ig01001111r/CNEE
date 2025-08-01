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