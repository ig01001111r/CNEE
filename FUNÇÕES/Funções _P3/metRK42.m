function [X, Y1, Y2] = metRK42(a, b, m, y10, y20, F1, F2, FA)

    % Entrada:
    % a, b  - intervalo da variável independente 
    % m     - número de passos
    % y10   - valor inicial de y1 (v(0))
    % y20   - valor inicial de y2 (v'(0))
    % F1    - função  y1'
    % F2    - função  y2' 
    % FA    - solução exata 
    
    % Saída:
    % X     - vetor com os valores de x calculados
    % Y1    - vetor com valores aproximados de y1
    % Y2    - vetor com valores aproximados de y2

    syms z k w;
    vet = [z, k, w]; 

    X = zeros(1, m+1);
    Y1 = zeros(1, m+1);
    Y2 = zeros(1, m+1);

    h = (b - a)/m;    
    xt = a; y1t = y10; y2t = y20;
    X(1) = xt; Y1(1) = y1t; Y2(1) = y2t;

    fprintf("i\t x \t\t v \t\t v' \t\t v exata\n");
    fprintf(' 0 \t %f \t %f \t %f\n', xt, y1t, y2t);

    for i = 1:m
       
        x = xt; y1 = y1t; y2 = y2t;
        k11 = vpa(subs(F1, vet, [x, y1, y2])); % calcula k1 para y1' no ponto (x,y1,y2)
        k12 = vpa(subs(F2, vet, [x, y1, y2])); % calcula k1 para y2' no ponto (x,y1,y2)

        x = xt + h/2; y1 = y1t + h/2*k11; y2 = y2t + h/2*k12;
        k21 = vpa(subs(F1, vet, [x, y1, y2]));
        k22 = vpa(subs(F2, vet, [x, y1, y2]));

        y1 = y1t + h/2*k21; y2 = y2t + h/2*k22;
        k31 = vpa(subs(F1, vet, [x, y1, y2]));
        k32 = vpa(subs(F2, vet, [x, y1, y2]));

        x = xt + h; y1 = y1t + h*k31; y2 = y2t + h*k32;
        k41 = vpa(subs(F1, vet, [x, y1, y2]));
        k42 = vpa(subs(F2, vet, [x, y1, y2]));

        xt = a + i*h;
        y1t = y1t + h/6*(k11 + 2*(k21 + k31) + k41);
        y2t = y2t + h/6*(k12 + 2*(k22 + k32) + k42);

        Fa = vpa(subs(FA, vet, [x, y1, y2]));  % calcula a solução exata no ponto (x,y1,y2)

        fprintf('%2d \t %f \t %f \t %f \t %f\n', i, xt, y1t, y2t, Fa);

        X(i+1) = xt; Y1(i+1) = y1t; Y2(i+1) = y2t;
    end
end
