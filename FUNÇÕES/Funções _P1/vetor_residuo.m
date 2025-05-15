function r = vetor_residuo(A, b, x)

  n = length(b);
  r = zeros(1, n);
  for i = 1:n
    soma = 0;
    for j = 1:n
      soma = soma + A(i, j) * x(j);
    end
    r(i) = b(i) - soma;
  end
end
