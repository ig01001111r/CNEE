function [x, interacoes] = refinar_LU(A, b, criterio, precisao)

    % Refinamento de solução de sistema linear Ax = b usando decomposição LU
    % A - matriz do sistema
    % b - vetor do termo independente
    % criterio - critério de parada baseado na norma do vetor de correção
    % precisao - número de casas decimais usadas nas operações

    [A_LU, ~, pivot] = decomp_LU(A); 
    interacoes = 0;
    [L, U] = extrair_LU_da_decomposta(A_LU);

    % Solução inicial
    y = subst_sucessivas_pivotal(L, b, pivot);
    x = subst_retroativas(U, y);
    x = round(x, precisao);  % Aplica precisão
    fprintf("\nx0 =\n\n"); 
    disp(transposta(x))

    % Refinamento iterativo
    flag = true;
    while flag
        r = vetor_residuo(A, b, x);
        r = round(r, precisao);

        y = subst_sucessivas_pivotal(L, r, pivot);
        c = subst_retroativas(U, y);
        c = round(c, precisao);

        x = x + c;
        x = round(x, precisao);
        
        fprintf("\nr%d =\n\n", interacoes);
        disp(transposta(r))
        fprintf("\nc%d =\n\n", interacoes); 
        disp(transposta(c))
        fprintf("\nx%d =\n\n", interacoes+1); 
        disp(transposta(x))     

        interacoes = interacoes + 1;

        if Norma_Infenito(c) < criterio
            flag = false;
        end
    end
end
