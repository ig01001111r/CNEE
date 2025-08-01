function interpolado = interpolacaoLagrange(numero_pontos, X, Y, valor_interpolar)
    % X, Y = vetores com os pontos conhecidos
    % numero_pontos = quantidade de pontos
    % valor_interpolar = valor no qual queremos interpolar

    interpolado = 0; % inicializa o valor interpolado

    for i = 1:numero_pontos
        numerador = 1;
        denominador = 1;
        
        for j = 1:numero_pontos
            if i ~= j
                numerador = numerador * (valor_interpolar - X(j));
                denominador = denominador * (X(i) - X(j));
            end
        end
        
        interpolado = interpolado + Y(i) * (numerador / denominador);
    end
end
