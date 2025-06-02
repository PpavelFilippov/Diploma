function [total_profit, total_ad_cost, u, u2, x_opt, t] = discreteProblemPoints(M, p, c, r, delta, alpha, U, s0, T, t_a, s_min, N)
    if N < T && N == 1
        error('Ошибка: Некорректное значение N')
    end
    
    dt = T / N;     % Размер шага по времени
    
    K = length(t_a);  % Число ограничений
    
    % Производные параметры
    a = alpha / M;
    pi_val = (p - c) * M;  % Доход на единицу доли рынка
    
    % Временная сетка
    t = linspace(0, T, N+1);
    
    % Индексы для t_a
    t_a_idx = zeros(1, K);
    for k = 1:K
        t_a_idx(k) = round(t_a(k) / dt) + 1;  % Индекс времени t_a(k)
    end
    
    % Дискретизация пространства состояний
    Ns = 2000;  % Увеличиваем для большей точности
    s_grid = linspace(0, 0.001, Ns);  % Ограничиваем сверху, так как s_min малы
    ds = s_grid(2) - s_grid(1);
    
    % Инициализация функции стоимости
    V = zeros(N+1, Ns);
    u_opt = zeros(N, Ns);  % Оптимальное управление
    
    % Граничное условие
    V(N+1, :) = 0;
    
    % Обратный проход
    for n = N:-1:1
        for j = 1:Ns
            s_n = s_grid(j);
            % Проверка достижимости будущих ограничений
            constraint_violated = false;
            for k = 1:K
                if n <= t_a_idx(k)
                    steps_to_k = t_a_idx(k) - n;
                    s_future = predict_s(s_n, steps_to_k, dt, a, delta, sqrt(U));
                    if s_future < s_min(k)
                        constraint_violated = true;
                        break;
                    end
                end
            end
            if constraint_violated
                V(n, j) = -Inf;
                continue;
            end
            
            % Оптимизация по u_n
            u_vals = linspace(0, sqrt(U), 2000);
            rewards = zeros(size(u_vals));
            for iu = 1:length(u_vals)
                u_n = u_vals(iu);
                % Следующее состояние
                s_next = s_n + dt * (a * u_n * (1 - s_n) - delta * s_n);
                s_next = max(0, min(1, s_next));
                % Интерполяция V_{n+1}(s_{n+1})
                if s_next <= s_grid(end) && s_next >= s_grid(1)
                    idx = floor((s_next - s_grid(1)) / ds) + 1;
                    idx = min(max(idx, 1), Ns-1);
                    w = (s_next - s_grid(idx)) / ds;
                    V_next = (1 - w) * V(n+1, idx) + w * V(n+1, idx+1);
                else
                    V_next = -Inf;
                end
                % Текущая награда
                reward = exp(-r * t(n)) * (pi_val * s_n - u_n^2) * dt + V_next;
                rewards(iu) = reward;
            end
            [V(n, j), idx_opt] = max(rewards);
            u_opt(n, j) = u_vals(idx_opt);
        end
    end
    
    % Прямой проход
    s = zeros(N+1, 1);
    u = zeros(N, 1);
    s(1) = s0;
    for n = 1:N
        if s(n) <= s_grid(end) && s(n) >= s_grid(1)
            idx = floor((s(n) - s_grid(1)) / ds) + 1;
            idx = min(max(idx, 1), Ns-1);
            w = (s(n) - s_grid(idx)) / ds;
            u_n = (1 - w) * u_opt(n, idx) + w * u_opt(n, idx+1);
        else
            u_n = sqrt(U);  % Максимальное усилие, если вне сетки
        end
        u(n) = u_n;
        s(n+1) = s(n) + dt * (a * u_n * (1 - s(n)) - delta * s(n));
        s(n+1) = max(0, min(1, s(n+1)));
    end
    u2 = u.^2;
    
    % Вычисление прибыли и расходов
    total_profit = sum(exp(-r * t(1:N)) .* (pi_val * s(1:N) - u2)' * dt);
    total_ad_cost = sum(u2 * dt);
    
    % Вывод результатов
    fprintf('Итоги:\n');
    fprintf('Общая дисконтированная прибыль: %.2f\n', total_profit);
    fprintf('Общие расходы на рекламу: %.2f\n', total_ad_cost);
    
    % Проверка ограничений
    fprintf('\nПроверка ограничений:\n');
    for k = 1:K
        s_tk = interp1(t, s, t_a(k), 'linear');
        fprintf('t_a%d = %.2f: s(t_a%d) = %.6f, s_min%d = %.6f, разница = %.6f\n', ...
            k, t_a(k), k, s_tk, k, s_min(k), s_tk - s_min(k));
    end
    
    % Визуализация
    figure;
    subplot(4,1,1);
    x_opt = s*M;

    plot(t, s * M, 'b-', 'LineWidth', 2);
    hold on;
    for k = 1:K
        plot(t_a(k), s_min(k) * M, 'ro', 'MarkerSize', 8);
    end
    title('Продажи x(t) = s(t) \cdot M');
    xlabel('t'); ylabel('x(t)');
    grid on;
    
    subplot(4,1,2);
    plot(t, s, 'b-', 'LineWidth', 2);
    hold on;
    for k = 1:K
        plot(t_a(k), s_min(k), 'ro', 'MarkerSize', 8);
    end
    title('Доля рынка s(t)');
    xlabel('t'); ylabel('s(t)');
    grid on;
    
    subplot(4,1,3);
    plot(t(1:N), u2, 'r-', 'LineWidth', 2);
    title('Расходы на рекламу u^2(t)');
    xlabel('t'); ylabel('u^2(t)');
    grid on;
    
    subplot(4,1,4);
    plot(t(1:N), u, 'g-', 'LineWidth', 2);
    title('Управление u(t)');
    xlabel('t'); ylabel('u(t)');
    grid on;
    

end

% Функция для оценки минимального s через n шагов при максимальном u
function s_future = predict_s(s_current, n_steps, dt, a, delta, u_max)
    s = s_current;
    for i = 1:n_steps
        s = s + dt * (a * u_max * (1 - s) - delta * s);
        s = max(0, min(1, s));
    end
    s_future = s;
end