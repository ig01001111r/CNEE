% L: matriz triangular inferior
% vet_ans: vetor resultado (por exemplo, X'Y na regressão múltipla)

function Y = suc_subst_regres(L, vet_ans)
    % Resolve sistemas de equações com matriz triangular inferior: L*Y = vet_ans
    ordem = length(L);
    
    Y = zeros(ordem, 1); % vetor coluna
    
    % Primeira variável
    Y(1) = vet_ans(1) / L(1,1);
    
    % Substituição sucessiva
    for i = 2:ordem
        soma = 0;
        for j = 1:i-1
            soma = soma + L(i,j) * Y(j);
        end
        Y(i) = (vet_ans(i) - soma) / L(i,i);
    end
end
