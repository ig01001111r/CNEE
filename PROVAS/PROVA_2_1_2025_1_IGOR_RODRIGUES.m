clc; clear; close all;
format shortG;

% ALUNO: IGOR RODRIGUES
% :3)

syms x;
F1 =  x^5 + 3.48*x^4 + 0.0155*x^3 - 4.4243*x^2 + 3.850655*x^1 - 0.891996;
symX = sym('x');

% PRIMEIRA QUESTAO
fprintf("Primeira Questão:\n\n");
coef1 = [1 3.48 0.0155 -4.4243 3.850655 -0.891996];
r = roots(coef1);
n1 = length(coef1)-1;

fprintf("Letra A:\n");
% LETRA A: Encontre o número de raízes positivas e negativas pela regra de sinais
fprintf("A regra dos sinais nos diz que o número de raízes reais positivas\n");
fprintf("é dado pela variação de sinais na sequência dos coeficientes\n");
fprintf("ou é menor que esse número por um inteiro par.\n");

fprintf ("\nNesse caso, ou há 1 raízes reais positivas,ou 0."); 

fprintf("\n\nSimilarmente, para determinar as raízes reais negativas, analisamos\n"); 
fprintf("a permanência de sinais na sequência dos coeficientes.\n");

fprintf ("\nNesse caso, ou há 2 raízes reais negativas, ou 0."); 


fprintf("\n\nLetra B:\n");
% LETRA B: Os limites de raízes positivas e negativas pelo método de Lagrange
fprintf("Como a equação dada é polinomial, podemos determinar um limite\n");
fprintf("de onde as raízes estarão (positivas e negativas) usando o Teorema de Lagrange."); 

fprintf ("\n\nPara esse polinômio, temos os limites das raízes positivas e\n");
fprintf("negativas, respectivamente, são:\n ");
L = limitesDeRaizes(n1, coef1);
disp(L);
% 0.18808       2.6417        -4.48     -0.41573

fprintf("\nLetra C:\n");
% LETRA C: Faça o gráfico da função ao longo das raízes reais e encontre o
% intervalo onde se encontram cada uma das raízes, que serão utilizados no item d
fprintf("As 3 primeiras figuras produzidas são das raízes encontradas na análise do gráfico.\n")

% Segundo meu mano (pablo negão) fazer o grafico pelo excel é mais facil kkkkk
plotF( 18,-20, 500, x, F1);% grafico da função 
plotF(0, 0.5, 500, x, F1);% raiz 1
plotF(-2.6, -2.1, 500, x, F1);% raiz 2
plotF(-2.2, -3, 500, x, F1);% raiz 3

fprintf("Analisando o gráfico do limite das raízes positivas, temos que ξ1 ∈ [0, 0.5]\n");
fprintf("\nAnalisando o gráfico do limite das raízes negativas, temos que ξ2 ∈ [-2.6, -2.1]\n");
fprintf("e que ξ3 ∈ [-3, -2.2]\n");

fprintf("\n\nLetra D:\n");
% LETRA D: Todas as raízes reais utilizando todos os métodos com precisão de 10e-5
toler = 10e-5;
maxIter = 50;
aproxRaiz = [-3.1, -2.3, 1.2]; % chute das raízes
w = length(aproxRaiz) - 1;
AB = zeros(w+1,2);
x0 = zeros(w);
m = zeros(w);

% Esse é um loop para calcular o intervalo onde a funcao troca de sinal
% para cada raiz aproximada (e guardar na matriz AB), além de
% calcular o x0 que será utilizado no método de Newton, e também o m da
% multiplicidade. Tudo para cada raiz aproximada com base no gráfico


for i=1:w+1
    [AB(i,1), AB(i,2), ~] = ondeTrocaSinalP(aproxRaiz(i), n1, coef1);
    x0(i) = AB(i,1)+(AB(i,2)-AB(i,1))/2;
    m(i) = multiplicidade(aproxRaiz(i), symX, F1);
end

fprintf("\nApós o uso da função ondeTrocaSinal para cada raiz aproximada\n");
fprintf("visualmente, pelo gráfico: ");
disp(aproxRaiz);
fprintf("Temos que os intervalos a serem utilizados (cada linha um intervalo), serão:\n");
disp(AB);

metodos = {@metBissecao, @metSecante, @metRegulaFalsi, @metPegaso, @metMuller, @metNewton, @metSchroder};
nomes = {'Bissecao'; 'Secante'; 'Regula Falsi'; 'Pegaso'; 'Muller'; 'Newton'; 'Schroder'}; 
labels = {'Método', 'Raiz encontrada', 'Iterações', 'Erro'};
resultados1 = cell(numel(metodos), w);
resultados2 = cell(numel(metodos), w);
resultados3 = cell(numel(metodos), w);
resultados4 = cell(numel(metodos), w);

% para a raiz -3.1
for i = 1:numel(metodos)    
    metodo = metodos{i};        
    if i == 6 % Se o método for Newton                 
        [raiz, iter, erro] = metodo(x0(1), toler, maxIter, symX, F1);      
    elseif i == 7 % Se o método for Schroder  
        [raiz, iter, erro] = metodo(m(1), x0(1), toler, maxIter, symX, F1); 
    else    
        [raiz, iter, erro] = metodo(AB(1,1), AB(1,2), toler, maxIter, symX, F1);        
    end           
    resultados1{i, 1} = raiz;       
    resultados1{i, 2} = iter;       
    resultados1{i, 3} = erro;
    C1 = [nomes, resultados1];
    C1 = [labels; C1];
end  
fprintf("\nPara a raiz chutada pelo gráfico da figura 4:\n");
disp(C1);
fprintf("Ou seja, temos que a raiz ou é muito próxima, ou é -2.45.\n"); 

% para a raiz -2.3
for i = 1:numel(metodos)    
    metodo = metodos{i};        
    if i == 6 % Se o método for Newton                 
        [raiz, iter, erro] = metodo(x0(2), toler, maxIter, symX, F1);      
    elseif i == 7 % Se o método for Schroder  
        [raiz, iter, erro] = metodo(m(2), x0(2), toler, maxIter, symX, F1); 
    else    
        [raiz, iter, erro] = metodo(AB(2,1), AB(2,2), toler, maxIter, symX, F1);        
    end           
    resultados2{i, 1} = raiz;       
    resultados2{i, 2} = iter;       
    resultados2{i, 3} = erro;
    C2 = [nomes, resultados2];
    C2 = [labels; C2];
end  
fprintf("\nPara a raiz chutada pelo gráficoda da figura 3:\n");
disp(C2);
fprintf("Ou seja, temos que a raiz ou é muito próxima, ou é -2.40.\n"); 

% para a raiz 1.2
for i = 1:numel(metodos)    
    metodo = metodos{i};        
    if i == 6 % Se o método for Newton                 
        [raiz, iter, erro] = metodo(x0(3), toler, maxIter, symX, F1);      
    elseif i == 7 % Se o método for Schroder  
        [raiz, iter, erro] = metodo(m(3), x0(3), toler, maxIter, symX, F1); 
    else    
        [raiz, iter, erro] = metodo(AB(3,1), AB(3,2), toler, maxIter, symX, F1);        
    end           
    resultados3{i, 1} = raiz;       
    resultados3{i, 2} = iter;       
    resultados3{i, 3} = erro;
    C3 = [nomes, resultados3];
    C3 = [labels; C3];
end  
fprintf("\nPara a raiz chutada pelo gráficoda da figura 2:\n");
disp(C3);
fprintf("Ou seja, temos que a raiz ou é muito próxima, ou é 0.37.\n"); 



fprintf("\n\nLetra E:\n");
% LETRA E: Em caso de raízes complexas, utilize o método de Muller, Secante ou Newton para 
% determinar pelo menos um par conjugado delas. Utilize o intervalo [1 - 1𝑖, 1 + 1𝑖]
fprintf("Como o polinômino é de quinto grau, e apenas encontramos 3 raízes reais,\n");
fprintf("fica claro que as duas raízes que faltam são um par complexo.\n");
fprintf("O intervalo utilizado para encontrá-los será [1 - 1𝑖, 1 + 1𝑖].\n");
ai = 1 - 1i;
bi = 1 + 1i;
[raizSI, iterSI, ~] = metSecante(ai, bi, toler, maxIter, symX, F1);   
fprintf("Pelo método da Secante, e usando %d iterações temos que a raiz procurada é: ", iterSI);
disp(raizSI);
fprintf("Logo, seu par conjugado é: 0.50000000179219180921380788546849 + 0.40000000031747202151071412791197i");

% SEGUNDA QUESTAO
% Determine, pelo método que desejar, a raiz da equação transcendental abaixo
fprintf("\n\nSegunda Questão:\n\n");

syms j;
F2 = sin(j) - log(j) -1;
symJ = sym('j');
plotF(0.5, 1.8, 500, symJ, F2);

fprintf("\nAnalizando o gráfico na (figura 5), temos que a função muda de sinal entre [0.6,0.8]\n");
fprintf("Logo, esse será o intervalo utilizado.");
[raizS, ~, ~] = metSecante(0.6, 0.8, toler, maxIter, symJ, F2); 
fprintf("\nUsando o método da Secante, temos que a raiz da função é: %.7f", raizS);

%% FUNÇÔES UTLIZADAS

function plotF(a, b, num, x, F)
    u = linspace(a, b, num);
    Fx = vpa(subs(F, x, u));
    figure;
    plot(u, Fx);
    grid on;
    hold on;
    plot(u, zeros(size(u)), 'r--');
end

function Pa = horner(n, coef, a)
    Pa = coef(1);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    for i=2:(n+1)
        Pa = Pa*a + coef(i); 
    end  
end

function [m, erro] = multiplicidade(raiz, x, F1)
    DFx = diff(F1, x);
    m = 1;
    max = 50;
    while (vpa(subs(DFx, x, raiz)) == 0)
        DFx = diff(DFx, x);
        m = m + 1;
        if m == max
            erro = 1;
            break;
        end
    end
    erro = 0;
end

function [a, b, erro] = ondeTrocaSinalP(z, n, coef) 
    if z == 0
        a = -0.05; 
        b = 0.05;
    else
        a = 0.95 * z;
        b = 1.05 * z;
    end
    iter = 0;
    aureo = 2/(sqrt(5) - 1);
    Fa = horner(n, coef, a);
    Fb = horner(n, coef, b);
    while(1)
        if ((Fa * Fb) <= 0) || (iter >= 20)
            break;
        end
        iter = iter+1;
        if(abs(Fa) < abs(Fb))
            a = a - aureo*(b-a);
            Fa = horner(n, coef, a);
        else
            b = b + aureo*(b-a);
            Fb = horner(n, coef, b);
        end
    end  
    if((Fa * Fb) <= 0)
        erro = 0;
    else
        erro = 1;
    end
end

function L = limitesDeRaizes(n, coef)
    L = zeros(1,4);
    if coef(1) == 0
        disp('O coef(1) é nulo');
        return;
    end
    t = n+1;
    coef(t+1) = 0;
    while(1)
        if coef(t) ~= 0
            break;    
        end
        t = t-1;
    end
    for i=1:4
        if (i == 2)||(i == 4)
            for j=1:(t/2)  
                aux = coef(j);
                coef(j) = coef(t-j+1);
                coef(t-j+1) = aux;
            end
        elseif i == 3
            for j=1:(t/2)  
                aux = coef(j);
                coef(j) = coef(t-j+1);
                coef(t-j+1) = aux;
            end  
            for j=(t-1):-2:1
                coef(j) = -coef(j);
            end    
        end
        if coef(1) < 0 
            for j=1:t
                coef(j) = -coef(j);
            end    
        end
        k = 2;
        while(1)
            if (coef(k) < 0)||(k > t)
                break;
            end
            k = k+1;
        end    
        if k <= t
            B = 0;
            for j=2:t
                if (coef(j) < 0) && (abs(coef(j)) > B)
                    B = abs(coef(j));
                end    
            end
            L(i) = 1 + (B/coef(1))^(1/(k-1));
        else
            L(i) = 10^100;
        end      
    end
    aux = L(1);
    L(1) = 1/L(2);
    L(2) = aux;
    L(3) = -L(3);
    L(4) = -1/L(4);
end

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

function [raiz, iter, erro] = metSecante(a, b, toler, maxIter, x, F1)  
    Fa = vpa(subs(F1, x, a));
    Fb = vpa(subs(F1, x, b));
    if abs(Fa) < abs(Fb)
        [a, b] = deal(b, a);
        [Fa, Fb] = deal(Fb, Fa);
    end   
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
        a = b;
        Fa = Fb;
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

function [raiz, iter, erro] = metRegulaFalsi(a, b, toler, maxIter, x, F1)
    Fa = vpa(subs(F1, x, a));
    Fb = vpa(subs(F1, x, b));   
    if (Fa * Fb) > 0
        error('A função não muda de sinal nos extremos do intervalo dado');
    end   
    if Fa > 0
        [a, b] = deal(b, a);
        [Fa, Fb] = deal(Fb, Fa);
    end  
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
        if Fx < 0
            a = u;
            Fa = Fx;
        else
            b = u;
            Fb = Fx;
        end        
        iter = iter + 1;
    end   
    raiz = u;
    if (abs(deltaX) <= toler) && (abs(Fx) <= toler)
        erro = 0;
    else
        erro = 1;
    end
end

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

function [raiz, iter, erro] = metSchroder(m, x0, toler, maxIter, x, F1)
    Fx = vpa(subs(F1, x, x0));
    DFx = vpa(subs(diff(F1,x), x, x0));
    u = x0;
    iter = 0;
    while (1)
        deltaX = -m*Fx/DFx;
        u = u + deltaX;
        Fx = vpa(subs(F1, x, u));
        DFx = vpa(subs(diff(F1,x), x, u));
        iter = iter + 1;
        if ((abs(deltaX) <= toler) && (abs(Fx) <= toler)) || DFx == 0 || iter>= maxIter
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