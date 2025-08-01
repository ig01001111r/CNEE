function interpolado = interpolacaoGregoryNewton(numero_pontos, X, Y, valor_interpolar)
    % X, Y: vetores com os pontos conhecidos
    % numero_pontos: quantidade de pontos
    % valor_interpolar: valor no qual queremos interpolar

    % Inicializa as diferenças finitas com os valores de Y
    for i = 1:numero_pontos
        Dely(i) = Y(i);
    end

    % Construção das diferenças finitas progressivas
    for k = 1:numero_pontos - 1
        for i = numero_pontos:-1:k+1
            Dely(i) = Dely(i) - Dely(i - 1);
        end
    end

    % Avaliação pelo processo de Horner
    u = (valor_interpolar - X(1)) / (X(2) - X(1));
    interpolado = Dely(numero_pontos);
    
    for i = numero_pontos - 1:-1:1
        interpolado = interpolado * (u - i + 1) / i + Dely(i);
    end
end
