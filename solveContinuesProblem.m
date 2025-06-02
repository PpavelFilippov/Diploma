function solveContinuesProblem(M, p, c, r, delta, alpha, U, s0, T, constraint_type, constraint_set, s_min)
    if strcmp(constraint_type, 'interval')
        if size(constraint_set, 2) ~= 2 || length(s_min) ~= size(constraint_set, 1)
            error('Для ограничения по интервалам constraint_set должен быть [K x 2], а s_min длины K.');
        end
        continuesProblemInterval(M, p, c, r, delta, alpha, U, s0, T, constraint_set, s_min);
    elseif strcmp(constraint_type, 'points')
        if length(constraint_set) ~= length(s_min)
            error('Для ограничения по точкам constraint_set и s_min должны иметь одинаковую длину.');
        end
        continuesProblemPoints(M, p, c, r, delta, alpha, U, s0, T, constraint_set, s_min);
    else
        error('constraint_type должен быть ''interval'' или ''points''.');
    end
end




