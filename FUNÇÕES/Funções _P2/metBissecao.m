% a = limite inferior
% b = limite superior
% Toler = tolerancia
% IterMax = número máximo de iterações
% Raiz = raiz
% Iter = numero de iterações realizadas

% condErro = condição de erro,
% condErro = 0 se a raiz foi encontrada
% condErro = 1 se a raiz não foi encontrada


function [raiz, iter, erro] = metBissecao(a, b, toler, maxIter, x, F1)  
    Fa = vpa(subs(F1, x, a));
    Fb = vpa(subs(F1, x, b));    
    if (Fa * Fb) > 0
        error('A função não muda de sinal nos extremos do intervalo dado');
    end    
    deltaX = abs(b-a)/2;
    iter = 0;
    while(1)
        u = (a + b)/2;
        Fx = vpa(subs(F1, x, u));
        if ((deltaX <= toler) && (abs(Fx) <= toler)) || iter >= maxIter
            break;
        end      
        if (Fa * Fx) > 0
            a = u;
            Fa = Fx;
        else
            b = u;
        end       
        deltaX = deltaX/2;
        iter = iter + 1;
    end   
    raiz = u;
    if (deltaX <= toler) && (abs(Fx) <= toler)
        erro = 0;
    else
        erro = 1;
    end 
end