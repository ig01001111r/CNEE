function [x, iter, condErro] = gauss_seidel(A, b, toler, iterMax)
    n = length(A);
    x = zeros(1, n);
    v = zeros(1, n);
    for i = 1:n
        r = 1/A(i,i);
        for j = 1:n
            if i ~= j
                A(i,j) = A(i,j) * r;
            end
        end
        b(i) = b(i) * r; x(i) = b(i);
    end
    iter = 0;
    
    % Representação das soluções baseados no livros
    fprintf("\nSolução do sistema linear pelo método de Gauss-Seidel\n" + ...
        "Iter    ")
    for i = 1:n
        if i == 2
            fprintf("x%d          ", i)
        elseif i == n
            fprintf("x%d       ", i)
        else
            fprintf("x%d         ", i)
        end
    end
    fprintf("normaRelativa\n")

    % Iterações de Gauss-Seidel
    while true
        iter = iter + 1;
        for i = 1:n
            soma = 0;
            for j = 1:n
                if i ~= j
                    soma = soma + A(i,j) * x(j);
                end
            end
            v(i) = x(i); x(i) = b(i) - soma;
        end
        normaNum = 0; normaDen = 0;
        for i = 1:n
            t = abs(x(i) - v(i));
            if t > normaNum
                normaNum = t;
            end
            if abs(x(i)) > normaDen
                normaDen = abs(x(i));
            end
        end
        normaRel = normaNum/normaDen;
        disp([sprintf(' %d  ', iter), sprintf('  %.5f  ', x), sprintf('  %.6e', normaRel)]);
        % Teste de convergência
        if normaRel <= toler || iter >= iterMax
            break
        end
    end
    if normaRel <= toler
        condErro = 0;
    else
        condErro = 1;
    end
end