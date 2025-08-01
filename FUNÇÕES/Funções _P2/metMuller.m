% a = limite inferior
% b = limite superior
% Toler = tolerancia
% IterMax = número máximo de iterações
% Raiz = raiz
% Iter = numero de iterações realizadas

% condErro = condição de erro,
% condErro = 0 se a raiz foi encontrada
% condErro = 1 se a raiz não foi encontrada


function [raiz, iter, erro] = metMuller(a, c, toler, maxIter, x, F1)
    Fa = vpa(subs(F1, x, a));
    Fc = vpa(subs(F1, x, c));
    b = (a + c)/2;
    Fb = vpa(subs(F1, x, b));
    u = b;
    Fx = Fb;
    deltaX = c - a;
    iter = 0;
    while(1)
        h1 = c - b;
        h2 = b - a;
        r = h1/h2;
        t = u;
        A = (Fc - (r + 1)*Fb + r*Fa)/(h1*(h1 + h2));
        B = (Fc - Fb)/h1 - A*h1;
        C = Fb;
        z = (-B + sign(B)*sqrt(B^2 - 4*A*C))/(2*A);
        u = b + z;
        deltaX = u - t;
        Fx = vpa(subs(F1, x, u));
        if (abs(deltaX) <= toler && abs(Fx) <= toler) || iter >= maxIter
            break;
        end
        if u > b
            a = b;
            Fa = Fb;
        else
            c = b;
            Fc = Fb;
        end
        b = u;
        Fb = Fx;
        iter = iter + 1;
    end
    raiz = u;
    if (abs(deltaX) <= toler) && (abs(Fx) <= toler)
        erro = 0;
    else
        erro = 1;
    end
end
