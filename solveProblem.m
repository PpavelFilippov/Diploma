function [total_profit, total_ad_cost, u, u2, x_opt, t] = solveProblem(M, p, c, r, delta, alpha, U, s0, T, model_type, constraint_type, constraint_set, s_min, N)
    if nargin < 13
        N = T;  % значение по умолчанию
    end

    % Проверка корректности constraint_set и s_min
    if strcmp(constraint_type, 'interval')
        if size(constraint_set, 2) ~= 2 || length(s_min) ~= size(constraint_set, 1)
            error('Для ограничения по интервалам constraint_set должен быть [K x 2], а s_min длины K.');
        end
    elseif strcmp(constraint_type, 'points')
        if length(constraint_set) ~= length(s_min)
            error('Для ограничения по точкам constraint_set и s_min должны иметь одинаковую длину.');
        end
    else
        error('constraint_type должен быть ''interval'' или ''points''.');
    end

    % Вызов нужной модели
    if strcmp(model_type, 'continuous')
        if strcmp(constraint_type, 'interval')
            [total_profit, total_ad_cost, u, u2, x_opt, t] = continuesProblemInterval(M, p, c, r, delta, alpha, U, s0, T, constraint_set, s_min);
        else
            [total_profit, total_ad_cost, u, u2, x_opt, t] = continuesProblemPoints(M, p, c, r, delta, alpha, U, s0, T, constraint_set, s_min);
        end
    elseif strcmp(model_type, 'discrete')
        if strcmp(constraint_type, 'interval')
            [total_profit, total_ad_cost, u, u2, x_opt, t] = discreteProblemInterval(M, p, c, r, delta, alpha, U, s0, T, constraint_set, s_min, N);
        else
            [total_profit, total_ad_cost, u, u2, x_opt, t] = discreteProblemPoints(M, p, c, r, delta, alpha, U, s0, T, constraint_set, s_min, N);
        end
    else
        error('model_type должен быть ''continuous'' или ''discrete''.');
    end
end
