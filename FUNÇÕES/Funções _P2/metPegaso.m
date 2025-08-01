% a = limite inferior
% b = limite superior
% Toler = tolerancia
% IterMax = número máximo de iterações
% Raiz = raiz
% Iter = numero de iterações realizadas

% condErro = condição de erro,
% condErro = 0 se a raiz foi encontrada
% condErro = 1 se a raiz não foi encontrada

function [raiz, iter, erro] = metPegaso(a, b, toler, maxIter, x, F1)   
    Fa = vpa(subs(F1, x, a));
    Fb = vpa(subs(F1, x, b));   
    iter = 0;
    u = b;
    Fx = Fb;
    while(1)
        deltaX = -Fx/(Fb - Fa)*(b - a);
        u = u + deltaX;
        Fx = vpa(subs(F1, x, u));
        if ((abs(deltaX) <= toler) && (abs(Fx) <= toler)) || iter >= maxIter
            break;
        end        
        if (Fx * Fb) < 0
            a = b;
            Fa = Fb;
        else
            Fa = (Fa * Fb)/(Fb + Fx);
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
