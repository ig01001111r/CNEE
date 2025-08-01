function plotFF(a, b, num, x, F, X, Y1)
    figure;
    hold on;
    u = linspace(a, b, num);
    Fx = vpa(subs(F, x, u));
    plot(u, Fx, 'r-', 'LineWidth', 2);
    scatter(X, Y1, 'b', 'filled');
    grid on;
    xlabel('Tempo (s)');
    ylabel('Tensão v(t) (V)');
    %title('Comparação RK4 (pontos azuis) e Solução Analítica (linha vermelha)');
    legend('Solução analítica', 'RK4 numérico');
    hold off;
end