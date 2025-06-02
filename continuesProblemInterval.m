function [total_profit, total_ad_cost, u, u2, x_opt, t] = continuesProblemInterval(M, p, c, r, delta, alpha, U, s0, T, intervals, s_min)
    % Определение интервалов и минимальных долей рынка
    K = size(intervals, 1);  % Число интервалов
    for k = 1:K
        if intervals(k, 1) >= intervals(k, 2)
            error('Ошибка: Начало интервала %d (%.2f) должно быть меньше конца (%.2f).', k, intervals(k, 1), intervals(k, 2));
        end
        if intervals(k, 1) < 0 || intervals(k, 2) > T
            error('Ошибка: Интервал %d выходит за пределы [0, %.2f].', k, T);
        end
    end
    for k = 1:K-1
        if intervals(k, 2) >= intervals(k+1, 1)
            error('Ошибка: Интервалы %d и %d пересекаются или касаются (%.2f >= %.2f).', ...
                k, k+1, intervals(k, 2), intervals(k+1, 1));
        end
    end

    % Добавляем последний интервал до T
    intervals = [intervals; [intervals(end, 2), T]];
    s_min = [s_min, 0];  % s_min = 0 для последнего интервала

    % Производные параметры
    a = alpha / M;
    pi_val = (p - c) * M;

    % Временная сетка
    N = 5000;
    dt = T / N;
    t = linspace(0, T, N+1);

    % Индексы интервалов
    interval_idx = zeros(K+1, 1);
    for k = 1:K+1
        if k == 1
            interval_idx(k) = 1;
        else
            interval_idx(k) = find(t >= intervals(k-1, 2), 1);
        end
    end

   

    % Инициализация
    s = zeros(N+1, 1);
    lambda = zeros(N+1, 1);
    u = zeros(N, 1);
    eta = zeros(N+1, 1);
    residual = zeros(N+1, 1);

    % Итерации
    max_iter = 200;
    tol = 1e-6;
    for iter = 1:max_iter
        s(1) = s0;
        for i = 1:N
            y = [s(i); lambda(i)];
            k1 = odes_int(t(i), y, u(i), eta(i), a, delta, pi_val, r);
            k2 = odes_int(t(i) + dt/2, y + dt*k1/2, u(i), eta(i), a, delta, pi_val, r);
            k3 = odes_int(t(i) + dt/2, y + dt*k2/2, u(i), eta(i), a, delta, pi_val, r);
            k4 = odes_int(t(i) + dt, y + dt*k3, u(i), eta(i), a, delta, pi_val, r);
            s(i+1) = s(i) + (dt/6) * (k1(1) + 2*k2(1) + 2*k3(1) + k4(1));
        end

        lambda(N+1) = 0;
        for i = N:-1:1
            y = [s(i); lambda(i+1)];
            k1 = odes_int(t(i), y, u(i), eta(i), a, delta, pi_val, r);
            k2 = odes_int(t(i) + dt/2, y - dt*k1/2, u(i), eta(i), a, delta, pi_val, r);
            k3 = odes_int(t(i) + dt/2, y - dt*k2/2, u(i), eta(i), a, delta, pi_val, r);
            k4 = odes_int(t(i) + dt, y - dt*k3, u(i), eta(i), a, delta, pi_val, r);
            lambda(i) = lambda(i+1) - (dt/6) * (k1(2) + 2*k2(2) + 2*k3(2) + k4(2));
        end

        for i = 1:N+1
            eta_i = 0;
            residual(i) = inf;
            active_interval = 0;
            for k = 1:K
                if t(i) >= intervals(k, 1) && t(i) <= intervals(k, 2)
                    residual(i) = s(i) - s_min(k);
                    eta_i = eta(i);
                    active_interval = k;
                    break;
                end
            end
        end

        step_size = 1e11 / (iter + 10);
        if max(abs(residual)) < tol
            fprintf('Сходимость достигнута на итерации %d\n', iter);
            break;
        end
        for i = 1:N+1
            for k = 1:K
                if t(i) >= intervals(k, 1) && t(i) <= intervals(k, 2)
                    eta(i) = max(0, eta(i) - step_size * residual(i));
                    break;
                end
            end
        end

        u = compute_u(s, lambda, t, a, U, r);
    end

    % Прибыль и расходы
    u2 = u.^2;
    total_profit = dt * sum(exp(-r * t(1:N)) .* (pi_val * s(1:N) - u2)');
    total_ad_cost = dt * sum(u2);
    fprintf('\nИтоги:\n');
    fprintf('Общая дисконтированная прибыль: %.2f\n', total_profit);
    fprintf('Общие расходы на рекламу: %.2f\n', total_ad_cost);

    fprintf('\nФинальная проверка ограничений:\n');
    for k = 1:K
        idx = find(t >= intervals(k, 1) & t <= intervals(k, 2));
        s_interval = s(idx);
        fprintf('Интервал [%.2f, %.2f]: min s(t)=%.6f, s_min=%.6f, разница=%.6f\n', ...
            intervals(k, 1), intervals(k, 2), min(s_interval), s_min(k), min(s_interval) - s_min(k));
    end

    % Визуализация
    figure;
    subplot(5,1,1);
    x_opt = s * M;
    plot(t, s * M, 'b-', 'LineWidth', 2); hold on;
    for k = 1:K
        idx = find(t >= intervals(k, 1) & t <= intervals(k, 2));
        plot(t(idx), s_min(k) * M * ones(size(t(idx))), 'r--', 'LineWidth', 1);
    end
    title('Продажи x(t) = s(t) \cdot M'); xlabel('t'); ylabel('x(t)'); grid on;

    subplot(5,1,2);
    plot(t, s, 'b-', 'LineWidth', 2); hold on;
    for k = 1:K
        idx = find(t >= intervals(k, 1) & t <= intervals(k, 2));
        plot(t(idx), s_min(k) * ones(size(t(idx))), 'r--', 'LineWidth', 1);
    end
    title('Доля рынка s(t)'); xlabel('t'); ylabel('s(t)'); grid on;

    subplot(5,1,3);
    plot(t(1:N), u2, 'r-', 'LineWidth', 2);
    title('Расходы на рекламу u^2(t)'); xlabel('t'); ylabel('u^2(t)'); grid on;

    subplot(5,1,4);
    plot(t, lambda, 'g-', 'LineWidth', 2);
    title('Сопряженная переменная \lambda(t)'); xlabel('t'); ylabel('\lambda(t)'); grid on;

    subplot(5,1,5);
    plot(t, eta, 'm-', 'LineWidth', 2);
    title('Множитель Лагранжа \eta(t)'); xlabel('t'); ylabel('\eta(t)'); grid on;
    
end

function u = compute_u(s, lambda, t, a, U, r)
    N = length(s) - 1;
    u = zeros(N, 1);
    for i = 1:N
        u_unconstrained = (lambda(i) * a * (1 - s(i))) / (2 * exp(-r * t(i)));
        u(i) = max(0, min(sqrt(U), u_unconstrained));
    end
end

function dydt = odes_int(t, y, u, eta, a, delta, pi_val, r)
    s = y(1); lambda = y(2);
    dsdt = a * u * (1 - s) - delta * s;
    dlambda_dt = -exp(-r * t) * pi_val + lambda * (a * u + delta) - eta;
    dydt = [dsdt; dlambda_dt];
end
