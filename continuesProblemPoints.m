function [total_profit, total_ad_cost, u, u22, x_opt, t] = continuesProblemPoints(M, p, c, r, delta, alpha, U, s0, T, t_a, s_min)
    N = 5000;       % Количество временных шагов
    dt = T / N;     % Размер шага по времени

    K = length(t_a);

    % Производные параметры
    a = alpha / M;
    pi_val = (p - c) * M;

    % Временная сетка
    t = linspace(0, T, N+1);

    % Индексы для t_a
    t_a_idx = zeros(1, K);
    for k = 1:K
        t_a_idx(k) = find(t >= t_a(k), 1);
    end

    % Инициализация
    s = zeros(N+1, 1);
    lambda = zeros(N+1, 1);
    u = zeros(N, 1);
    mu = zeros(1, K);

    max_iter = 200;
    tol = 1e-6;

    for iter = 1:max_iter
        % Прямое интегрирование s(t)
        s(1) = s0;
        for i = 1:N
            y = [s(i); lambda(i)];
            k1 = odes(t(i), y, u(i), a, delta, pi_val, r);
            k2 = odes(t(i) + dt/2, y + dt*k1/2, u(i), a, delta, pi_val, r);
            k3 = odes(t(i) + dt/2, y + dt*k2/2, u(i), a, delta, pi_val, r);
            k4 = odes(t(i) + dt, y + dt*k3, u(i), a, delta, pi_val, r);
            s(i+1) = s(i) + (dt/6) * (k1(1) + 2*k2(1) + 2*k3(1) + k4(1));
        end

        % Обратное интегрирование lambda(t)
        lambda(N+1) = 0;
        for k = K:-1:1
            if k == K
                end_idx = N+1;
            else
                end_idx = t_a_idx(k+1);
            end
            for i = end_idx-1:-1:t_a_idx(k)
                y = [s(i); lambda(i+1)];
                k1 = odes(t(i), y, u(i), a, delta, pi_val, r);
                k2 = odes(t(i) + dt/2, y - dt*k1/2, u(i), a, delta, pi_val, r);
                k3 = odes(t(i) + dt/2, y - dt*k2/2, u(i), a, delta, pi_val, r);
                k4 = odes(t(i) + dt, y - dt*k3, u(i), a, delta, pi_val, r);
                lambda(i) = lambda(i+1) - (dt/6) * (k1(2) + 2*k2(2) + 2*k3(2) + k4(2));
            end
            lambda(t_a_idx(k)) = lambda(t_a_idx(k)) + mu(k);
        end
        for i = t_a_idx(1)-1:-1:1
            y = [s(i); lambda(i+1)];
            k1 = odes(t(i), y, u(i), a, delta, pi_val, r);
            k2 = odes(t(i) + dt/2, y - dt*k1/2, u(i), a, delta, pi_val, r);
            k3 = odes(t(i) + dt/2, y - dt*k2/2, u(i), a, delta, pi_val, r);
            k4 = odes(t(i) + dt, y - dt*k3, u(i), a, delta, pi_val, r);
            lambda(i) = lambda(i+1) - (dt/6) * (k1(2) + 2*k2(2) + 2*k3(2) + k4(2));
        end

        % Проверка ограничений
        residual = zeros(1, K);
        for k = 1:K
            s_tk = interp1(t, s, t_a(k), 'linear');
            residual(k) = s_tk - s_min(k);
            fprintf('Итерация %d, t_a%d = %.2f: s(t_a%d) = %.6f, s_min%d = %.6f, residual = %.6f, mu%d = %.2f\n', ...
                iter, k, t_a(k), k, s_tk, k, s_min(k), residual(k), k, mu(k));
        end

        step_size = 1e11 / (iter + 10);
        if max(abs(residual)) < tol
            fprintf('Сходимость достигнута на итерации %d\n', iter);
            break;
        end
        for k = 1:K
            mu(k) = max(0, mu(k) - step_size * residual(k));
        end

        u = compute_u(s, lambda, t, a, U, r);
    end

    % Результаты
    u2 = u.^2;
    u22 = [u2;0];
    total_profit = dt * sum(exp(-r * t(1:N)) .* (pi_val * s(1:N) - u2)');
    total_ad_cost = dt * sum(u2);
    fprintf('\nИтоги:\n');
    fprintf('Общая дисконтированная прибыль: %.2f\n', total_profit);
    fprintf('Общие расходы на рекламу: %.2f\n', total_ad_cost);

    fprintf('\nФинальная проверка ограничений:\n');
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
    plot(t, lambda, 'g-', 'LineWidth', 2);
    hold on;
    for k = 1:K
        plot(t_a(k), lambda(t_a_idx(k)), 'ko', 'MarkerSize', 8);
    end
    title('Сопряженная переменная \lambda(t)');
    xlabel('t'); ylabel('\lambda(t)');
    grid on;
    
end

function dydt = odes(t, y, u, a, delta, pi_val, r)
    s = y(1);
    lambda = y(2);
    dsdt = a * u * (1 - s) - delta * s;
    dlambda_dt = -exp(-r * t) * pi_val + lambda * (a * u + delta);
    dydt = [dsdt; dlambda_dt];
end

function u = compute_u(s, lambda, t, a, U, r)
    N = length(s) - 1;
    u = zeros(N, 1);
    for i = 1:N
        u_unconstrained = (lambda(i) * a * (1 - s(i))) / (2 * exp(-r * t(i)));
        u(i) = max(0, min(sqrt(U), u_unconstrained));
    end
end
