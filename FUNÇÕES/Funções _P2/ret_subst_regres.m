% U: matriz triangular superior
% vet_ans: vetor resultado (X'Y na regressão múltipla)

function X = ret_subst_regres(U, vet_ans)
    order = length(U);
    X = zeros(order, 1); % vetor coluna, padrão em MATLAB
    X(order) = vet_ans(order) / U(order, order);
    
    for i = order-1:-1:1
        soma = 0;
        for j = i+1:order
            soma = soma + U(i, j) * X(j);
        end
        X(i) = (vet_ans(i) - soma) / U(i, i);
    end
end
end