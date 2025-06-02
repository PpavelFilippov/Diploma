% % Скрипт для создания структуры параметров, запуска приложения и инициализации данных
% 
% % 1. Создание структуры параметров с полями по умолчанию
% params = struct(...
%     'modelType', 'Непрерывная', ... % Тип модели: Непрерывная или Дискретная
%     'constraintType', 'Точечные', ... % Тип ограничений: Точечные или Интервальные
%     'M', 50000, ... % Максимальный размер рынка
%     'p', 25, ... % Цена за единицу
%     'c', 15, ... % Переменные затраты на единицу
%     'r', 0.06, ... % Ставка дисконтирования
%     'delta', 0.7, ... % Коэффициент затухания
%     'alpha', 1.5, ... % Эффективность рекламы
%     'U', 200000, ... % Максимальный бюджет на рекламу
%     's0', 0, ... % Начальная доля рынка
%     'T', 15, ... % Горизонт времени
%     'N', 150, ... % Количество временных шагов
%     't_a', [], ... % Положения точек или интервалов
%     'x_min', [] ... % Значения ограничений
% );
% 
% 
% % Вывод структуры по умолчанию
% disp('Структура параметров по умолчанию:');
% disp(params);
% 
% % 2. Запуск приложения
% disp('Запуск приложения для ввода параметров...');
% app = ModelParametersApp;
% uiwait(app.UIFigure); % Ожидание закрытия окна приложения
% 
% 
% % 3. Загрузка введенных пользователем данных
% if exist('model_params.mat', 'file')
%     load('model_params.mat', 'params');
%     disp('Параметры, введенные пользователем:');
%     disp(params);
% else
%     disp('Параметры не были сохранены. Используются значения по умолчанию.');
% end
% 
% % 4. Инициализация полей структуры (уже выполнена при загрузке)
% % Дополнительно можно проверить корректность данных, если требуется
% if isempty(params.t_a) || isempty(params.x_min)
%     disp('Предупреждение: массивы t_a или x_min пусты.');
% end
% 
% 
% 
% 
% disp('Инициализация параметров завершена.');
% 
% % 5. Вызов solveProblem (заглушка, предполагается, что функция существует)
% try
%      [total_profit, total_ad_cost, u, u2, x_opt, t] = solveProblem(params.M, params.p, params.c, params.r, params.delta, params.alpha, ...
%          params.U, params.s0, params.T, params.modelType, params.constraintType, ...
%         params.t_a, params.x_min / params.M, params.N);
% 
%     % t = linspace(0, params.T, params.N);
%                 % y1 = sin(t); % Данные для первого графика
%                 % y2 = cos(t); % Данные для второго графика
%                 y1 = x_opt;
%                 y2 = u2;
%                 a1 = params.T / 2; % Половина горизонта времени
%                 t2 = t(1:params.N);
%                 tt = t;
%                 if strcmp(params.modelType, 'discrete')
%                     tt = t2;
%                 end
%                 a2 = params.N / 10; % Десятая часть шагов
%                 a3 = params.M / 1000; % Тысячная часть рынка
%                 app.updatePlot(t, tt, y1, y2, a1, a2, a3);
% catch e
%     disp(['Ошибка при вызове solveProblem или построении графика: ' e.message]);
% end
% 
% % Удаление временного файла
% if exist('model_params.mat', 'file')
%     delete('model_params.mat');
% end

% solveProblem(M, p, c, r, delta, alpha, U, s0, T, 'discrete', 'points', points, s_min_points, N);

% Скрипт для создания структуры параметров, запуска приложения и инициализации данных

% 1. Создание структуры параметров с полями по умолчанию
params = struct(...
    'modelType', 'Непрерывная', ... % Тип модели: Непрерывная или Дискретная
    'constraintType', 'Точечные', ... % Тип ограничений: Точечные или Интервальные
    'M', 50000, ... % Максимальный размер рынка
    'p', 25, ... % Цена за единицу
    'c', 15, ... % Переменные затраты на единицу
    'r', 0.06, ... % Ставка дисконтирования
    'delta', 0.7, ... % Коэффициент затухания
    'alpha', 1.5, ... % Эффективность рекламы
    'U', 200000, ... % Максимальный бюджет на рекламу
    's0', 0, ... % Начальная доля рынка
    'T', 15, ... % Горизонт времени
    'N', 150, ... % Количество временных шагов
    't_a', [], ... % Положения точек или интервалов
    'x_min', [] ... % Значения ограничений
);

% Вывод структуры по умолчанию
disp('Структура параметров по умолчанию:');
disp(params);

% 2. Запуск приложения
disp('Запуск приложения для ввода параметров...');
app = ModelParametersApp;

% Цикл для обработки ввода параметров и обновления графиков
while isvalid(app.UIFigure)
    % Ожидание нажатия кнопки "Подтвердить"
    uiwait(app.UIFigure);
    
    % Проверка, было ли окно закрыто
    if ~isvalid(app.UIFigure)
        disp('Приложение закрыто пользователем.');
        break;
    end
    
    % 3. Загрузка введенных пользователем данных
    if exist('model_params.mat', 'file')
        load('model_params.mat', 'params');
        disp('Параметры, введенные пользователем:');
        disp(params);
    else
        disp('Параметры не были сохранены. Используются значения по умолчанию.');
    end
    
    % 4. Проверка корректности данных
    if isempty(params.t_a) || isempty(params.x_min)
        disp('Предупреждение: массивы t_a или x_min пусты.');
    end
    
    disp('Инициализация параметров завершена.');
    
    % 5. Вызов solveProblem и обновление графиков
    try
        [total_profit, total_ad_cost, u, u2, x_opt, t] = solveProblem(params.M, params.p, params.c, params.r, params.delta, params.alpha, ...
            params.U, params.s0, params.T, params.modelType, params.constraintType, ...
            params.t_a, params.x_min / params.M, params.N);
        
        % Подготовка данных для графиков
                y1 = x_opt;
                y2 = u2;
                a1 = total_profit; % Половина горизонта времени
                t2 = t(1:params.N);
                tt = t;
                if strcmp(params.modelType, 'discrete')
                    tt = t2;
                else
                    if strcmp(params.constraintType, 'interval')
                        NN = 5000;
                        t2 = t(1:NN);
                        tt = t2;
                    end
                end

                a2 = total_ad_cost; % Десятая часть шагов
                a3 = params.M / 1000; % Тысячная часть рынка
                app.updatePlot(t, tt, y1, y2, a1, a2, a3);
        disp('Графики успешно обновлены.');
    catch e
        disp(['Ошибка при вызове solveProblem или построении графика: ' e.message]);
    end
    
    % Удаление временного файла
    if exist('model_params.mat', 'file')
        delete('model_params.mat');
    end
end

% Завершение работы скрипта
disp('Скрипт завершил выполнение.');