% --------------------------------------------------------------
% Modelo original  | Transformação aplicada
% --------------------------------------------------------------
% Exponencial:
% y = a * exp(b * t)
% ln(y) = ln(a) + b * t
% Y = ln(y), X = t

% Decaimento exponencial:
% y = a * exp(-b * t)
% ln(y) = ln(a) - b * t
% Y = ln(y), X = t

% Potência:
% y = a * t^b
% ln(y) = ln(a) + b * ln(t)
% Y = ln(y), X = ln(t)

% Raiz quadrada:
% y = a * sqrt(t) + b
% Linear em x = sqrt(t)
% Y = y, X = sqrt(t)

% Inverso:
% y = a / t + b
% Linear em x = 1/t
% Y = y, X = 1/t

% Logarítmico:
% y = a + b * ln(t)
% Linear
% Y = y, X = ln(t)

% Exponencial com potência:
% y = a * t^b * exp(c * t)
% ln(y) = ln(a) + b * ln(t) + c * t
% Y = ln(y), X = [ln(t), t]

% Modelo saturação (Michaelis-Menten):
% y = a * t / (b + t)
% Não linear, usa regressão não linear direta

% Crescimento logístico:
% y = L / (1 + exp(-k * (t - t0)))
% Não linear, usa regressão não linear direta

% -----------------------------
% Regressão Exponencial: R(t) = a·e^(−λ·t)
% Linearizamos: ln(R(t)) = ln(a) − λ·t
% -----------------------------
% -----------------------------
% Regressão Exponencial (Decaimento): R(t) = a * e^(−λ * t)
% Passo a passo da linearização:
% 1) Tomamos o log natural de ambos os lados:
%    ln(R(t)) = ln(a * e^(−λ * t))
% 2) Usamos propriedades do log:
%    ln(R(t)) = ln(a) + ln(e^(−λ * t)) = ln(a) − λ * t
% 3) Definimos variáveis lineares:
%    Y = ln(R), X = t
% 4) Forma linear: Y = ln(a) − λ * t
% -----------------------------

lnR = log(R);
X_exp = [ones(n,1), -t];
[b_exp, R2_exp, sigma2_exp] = regresPolinomial(n, 2, X_exp, lnR);
a_exp = exp(b_exp(1));   % a = e^(b0)
lambda_exp = b_exp(2);   % λ = b1

% -----------------------------
% Regressão Potência: R(t) = a * t^b
% Passo a passo da linearização:
% 1) Log natural dos dois lados:
%    ln(R(t)) = ln(a * t^b)
% 2) Propriedades do log:
%    ln(R(t)) = ln(a) + b * ln(t)
% 3) Variáveis lineares:
%    Y = ln(R), X = ln(t)
% 4) Forma linear: Y = ln(a) + b * ln(t)
% -----------------------------
lnR = log(R);
lnT = log(t);
X_pot = [ones(n,1), lnT];
[b_pot, R2_pot, sigma2_pot] = regresPolinomial(n, 2, X_pot, lnR);
a_pot = exp(b_pot(1));   % a = e^(b0)
b_pot_coef = b_pot(2);   % b = b1

% -----------------------------
% Regressão Linear simples: R(t) = a + b * t
% Já é linear, basta definir:
% Y = R, X = t
% -----------------------------
X_lin = [ones(n,1), t];
[b_lin, R2_lin, sigma2_lin] = regresPolinomial(n, 2, X_lin, R);
a_lin = b_lin(1);
b_lin_coef = b_lin(2);

% -----------------------------
% Regressão Polinomial grau 2: R(t) = a + b * t + c * t^2
% Já linear nos coeficientes, definimos:
% Y = R, X = [t, t^2]
% -----------------------------
X_poly2 = [ones(n,1), t, t.^2];
[b_poly2, R2_poly2, sigma2_poly2] = regresPolinomial(n, 3, X_poly2, R);
a_poly2 = b_poly2(1);
b_poly2_coef = b_poly2(2);
c_poly2_coef = b_poly2(3);

% -----------------------------
% Regressão Logarítmica: R(t) = a + b * ln(t)
% Já linear, definimos:
% Y = R, X = ln(t)
% -----------------------------
lnT = log(t);
X_log = [ones(n,1), lnT];
[b_log, R2_log, sigma2_log] = regresPolinomial(n, 2, X_log, R);
a_log = b_log(1);
b_log_coef = b_log(2);

% -----------------------------
% Regressão Raiz Quadrada: R(t) = a + b * sqrt(t)
% Linear em sqrt(t), definimos:
% Y = R, X = sqrt(t)
% -----------------------------
sqrtT = sqrt(t);
X_sqrt = [ones(n,1), sqrtT];
[b_sqrt, R2_sqrt, sigma2_sqrt] = regresPolinomial(n, 2, X_sqrt, R);
a_sqrt = b_sqrt(1);
b_sqrt_coef = b_sqrt(2);

% -----------------------------
% Regressão Inversa: R(t) = a + b / t
% Linear em 1/t, definimos:
% Y = R, X = 1/t
% -----------------------------
invT = 1 ./ t;
X_inv = [ones(n,1), invT];
[b_inv, R2_inv, sigma2_inv] = regresPolinomial(n, 2, X_inv, R);
a_inv = b_inv(1);
b_inv_coef = b_inv(2);

% -----------------------------
% Exponencial com potência: R(t) = a * t^b * e^(c * t)
% Linearização:
% 1) ln(R(t)) = ln(a) + b * ln(t) + c * t
% 2) Variáveis lineares:
%    Y = ln(R), X = [ln(t), t]
% -----------------------------
lnR = log(R);
lnT = log(t);
X_exp_pot = [ones(n,1), lnT, t];
[b_exp_pot, R2_exp_pot, sigma2_exp_pot] = regresPolinomial(n, 3, X_exp_pot, lnR);
a_exp_pot = exp(b_exp_pot(1));
b_exp_pot_coef = b_exp_pot(2);
c_exp_pot_coef = b_exp_pot(3);

