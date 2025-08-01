
% x0 = chute inicial 
% Toler = tolerancia
% IterMax = número máximo de iterações
% Raiz = raiz
% Iter = numero de iterações realizadas

% condErro = condição de erro,
% condErro = 0 se a raiz foi encontrada
% condErro = 1 se a raiz não foi encontrada

function [raiz, iter, erro] = metNewton(x0, toler, maxIter, x, F1)
    Fx = vpa(subs(F1, x, x0));
    DFx = vpa(subs(diff(F1,x), x, x0));
    u = x0;
    iter = 0;
    while (1)
        deltaX = -Fx/DFx;
        u = u + deltaX;
        Fx = vpa(subs(F1, x, u));
        DFx = vpa(subs(diff(F1,x), x, u));
        iter = iter + 1;       
        if ((abs(deltaX) <= toler) && (abs(Fx) <= toler)) || (DFx == 0 || iter>= maxIter)
            break;
        end
    end
    raiz = u;
    if (abs(deltaX) <= toler) && (abs(Fx) <= toler)
        erro = 0;
    else
        erro = 1;
    end
end