classdef ModelParametersApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        TaskTypeDropDownLabel   matlab.ui.control.Label
        TaskTypeDropDown        matlab.ui.control.DropDown
        MLabel                  matlab.ui.control.Label
        MEditField              matlab.ui.control.NumericEditField
        pLabel                  matlab.ui.control.Label
        pEditField              matlab.ui.control.NumericEditField
        cLabel                  matlab.ui.control.Label
        cEditField              matlab.ui.control.NumericEditField
        rLabel                  matlab.ui.control.Label
        rEditField              matlab.ui.control.NumericEditField
        deltaLabel              matlab.ui.control.Label
        deltaEditField          matlab.ui.control.NumericEditField
        alphaLabel              matlab.ui.control.Label
        alphaEditField          matlab.ui.control.NumericEditField
        ULabel                  matlab.ui.control.Label
        UEditField              matlab.ui.control.NumericEditField
        s0Label                 matlab.ui.control.Label
        s0EditField             matlab.ui.control.NumericEditField
        TLabel                  matlab.ui.control.Label
        TEditField              matlab.ui.control.NumericEditField
        NLabel                  matlab.ui.control.Label
        NEditField              matlab.ui.control.NumericEditField
        taLabel                 matlab.ui.control.Label
        taEditField             matlab.ui.control.EditField
        xminLabel               matlab.ui.control.Label
        xminEditField           matlab.ui.control.EditField
        ConfirmButton           matlab.ui.control.Button
        SaveParamsButton        matlab.ui.control.Button % Added
        LoadParamsButton        matlab.ui.control.Button % Added
        SaveArraysButton        matlab.ui.control.Button % Added
        MessageLabel            matlab.ui.control.Label
        PlotAxes                matlab.ui.control.UIAxes % First axes for plotting
        PlotAxes2               matlab.ui.control.UIAxes % Second axes for plotting
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: ConfirmButton
        function ConfirmButtonPushed(app, event)
            % Clear previous messages
            app.MessageLabel.Text = '';

            % Collect model parameters
            try
                params.M = app.MEditField.Value;
                params.p = app.pEditField.Value;
                params.c = app.cEditField.Value;
                params.r = app.rEditField.Value;
                params.delta = app.deltaEditField.Value;
                params.alpha = app.alphaEditField.Value;
                params.U = app.UEditField.Value;
                params.s0 = app.s0EditField.Value/params.M;
                params.T = app.TEditField.Value;
                params.N = app.NEditField.Value;

                % Validate model parameters
                if any([params.M, params.p, params.c, params.r, params.delta, ...
                        params.alpha, params.U, params.T, params.N] <= 0)
                    uialert(app.UIFigure, 'Все параметры модели должны быть положительными.', 'Ошибка');
                    return;
                end
                if params.s0 < 0
                    uialert(app.UIFigure, 'Начальная доля рынка (s0) не может быть отрицательной.', 'Ошибка');
                    return;
                end

                % Collect task type (convert Russian to English)
                switch app.TaskTypeDropDown.Value
                    case 'Непрерывная модель, точечные ограничения'
                        params.modelType = 'continuous';
                        params.constraintType = 'points';
                    case 'Непрерывная модель, интервальные ограничения'
                        params.modelType = 'continuous';
                        params.constraintType = 'interval';
                    case 'Дискретная модель, точечные ограничения'
                        params.modelType = 'discrete';
                        params.constraintType = 'points';
                    case 'Дискретная модель, интервальные ограничения'
                        params.modelType = 'discrete';
                        params.constraintType = 'interval';
                    otherwise
                        uialert(app.UIFigure, 'Ошибка: Неверный тип задачи.', 'Ошибка');
                        return;
                end

                % Process constraints
                try
                    if strcmp(params.constraintType, 'points')
                        % Parse t_a and x_min
                        t_a = str2num(app.taEditField.Value); %#ok<ST2NM>
                        x_min = str2num(app.xminEditField.Value); %#ok<ST2NM>

                        % Validate lengths
                        if length(t_a) ~= length(x_min)
                            uialert(app.UIFigure, ['Массивы t_a и x_min должны быть одинаковой длины. ' ...
                                'Пожалуйста, исправьте ввод.'], 'Ошибка');
                            return;
                        end
                        if any(t_a < 0) || any(t_a > params.T)
                            uialert(app.UIFigure, ['Значения t_a должны быть в пределах [0, T]. ' ...
                                'Пожалуйста, исправьте ввод.'], 'Ошибка');
                            return;
                        end
                        if any(x_min < 0)
                            uialert(app.UIFigure, ['Значения x_min должны быть неотрицательными. ' ...
                                'Пожалуйста, исправьте ввод.'], 'Ошибка');
                            return;
                        end

                        params.t_a = t_a;
                        params.x_min = x_min;

                    else % interval
                        % Parse t_a as a 2D array
                        t_a = str2num(app.taEditField.Value); %#ok<ST2NM>
                        x_min = str2num(app.xminEditField.Value); %#ok<ST2NM>

                        % Validate t_a structure
                        if size(t_a, 2) ~= 2
                            uialert(app.UIFigure, ['Массив t_a должен быть двумерным с двумя столбцами ' ...
                                '(начало и конец интервалов). Пожалуйста, исправьте ввод.'], 'Ошибка');
                            return;
                        end
                        if length(x_min) ~= size(t_a, 1)
                            uialert(app.UIFigure, ['Количество интервалов в t_a и значений x_min ' ...
                                'должны совпадать. Пожалуйста, исправьте ввод.'], 'Ошибка');
                            return;
                        end
                        % Check that start < end for each interval
                        if any(t_a(:,1) >= t_a(:,2))
                            uialert(app.UIFigure, ['Начало интервала должно быть меньше конца. ' ...
                                'Пожалуйста, исправьте ввод.'], 'Ошибка');
                            return;
                        end
                        % Check that intervals do not overlap
                        for i = 1:size(t_a, 1)-1
                            if t_a(i,2) >= t_a(i+1,1)
                                uialert(app.UIFigure, ['Интервалы не должны пересекаться. ' ...
                                    'Пожалуйста, исправьте ввод.'], 'Ошибка');
                                return;
                            end
                        end
                        if any(t_a(:) < 0) || any(t_a(:) > params.T)
                            uialert(app.UIFigure, ['Значения t_a должны быть в пределах [0, T]. ' ...
                                'Пожалуйста, исправьте ввод.'], 'Ошибка');
                            return;
                        end
                        if any(x_min < 0)
                            uialert(app.UIFigure, ['Значения x_min должны быть неотрицательными. ' ...
                                'Пожалуйста, исправьте ввод.'], 'Ошибка');
                            return;
                        end

                        params.t_a = t_a;
                        params.x_min = x_min;
                    end
                catch
                    uialert(app.UIFigure, ['Ошибка в формате массивов t_a или x_min. ' ...
                        'Пожалуйста, введите корректные числовые массивы.'], 'Ошибка');
                    return;
                end

                % Save parameters to file
                try
                    save('model_params.mat', 'params');
                    app.MessageLabel.Text = 'Параметры успешно сохранены! Ожидается обработка...';
                    disp('Параметры сохранены в model_params.mat'); % Debug message
                    
                    % Resume execution
                    uiresume(app.UIFigure);
                catch e
                    uialert(app.UIFigure, ['Ошибка при сохранении параметров: ' e.message], 'Ошибка');
                    disp(['Ошибка сохранения: ' e.message]); % Debug message
                end
            catch e
                uialert(app.UIFigure, ['Ошибка ввода параметров: ' e.message], 'Ошибка');
                disp(['Ошибка ввода: ' e.message]); % Debug message
            end
        end
    end

    % Public methods
    methods (Access = public)

        % Update plots with new data and display additional variables
        function updatePlot(app, t1, t2, y1, y2, a1, a2, a3)
            try
                % Clear previous plots
                cla(app.PlotAxes);
                cla(app.PlotAxes2);
                
                % Plot first data (e.g., sine wave)
                plot(app.PlotAxes, t1, y1, 'b-', 'LineWidth', 2);
                app.PlotAxes.Title.String = 'Оптимальный объем продаж';
                app.PlotAxes.XLabel.String = 't';
                app.PlotAxes.YLabel.String = 'x(t)';
                grid(app.PlotAxes, 'on');
                
                % Plot second data (e.g., cosine wave)
                plot(app.PlotAxes2, t2, y2, 'r-', 'LineWidth', 2);
                app.PlotAxes2.Title.String = 'Оптимальные рекламные расходы';
                app.PlotAxes2.XLabel.String = 't';
                app.PlotAxes2.YLabel.String = 'u^2(t)';
                grid(app.PlotAxes2, 'on');
                
                % Display message with values
                app.MessageLabel.Text = sprintf('Прибыль J: %.2f, суммарные рекламные расходы U: %.2f', a1, a2);
                disp('Графики успешно обновлены в приложении'); % Debug message
            catch e
                app.MessageLabel.Text = ['Ошибка при построении графиков: ' e.message];
                disp(['Ошибка в updatePlot: ' e.message]); % Debug message
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100, 100, 1200, 800]; % Set fixed window size
            app.UIFigure.Name = 'Параметры модели';

            % Define positions based on fixed window size
            leftMargin = 50;
            topMargin = 750;
            labelWidth = 200;
            fieldWidth = 200;
            dropdownWidth = 300;
            rowHeight = 50;
            buttonWidth = 110; % Twice as narrow as original 100
            buttonHeight = 22; % Same as ConfirmButton
            buttonGap = 10;

            % Create TaskTypeDropDownLabel
            app.TaskTypeDropDownLabel = uilabel(app.UIFigure);
            app.TaskTypeDropDownLabel.HorizontalAlignment = 'right';
            app.TaskTypeDropDownLabel.Position = [leftMargin, topMargin - rowHeight*0, labelWidth, 22];
            app.TaskTypeDropDownLabel.Text = 'Тип задачи';

            % Create TaskTypeDropDown
            app.TaskTypeDropDown = uidropdown(app.UIFigure);
            app.TaskTypeDropDown.Items = {'Непрерывная модель, точечные ограничения', ...
                                         'Непрерывная модель, интервальные ограничения', ...
                                         'Дискретная модель, точечные ограничения', ...
                                         'Дискретная модель, интервальные ограничения'};
            app.TaskTypeDropDown.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*0, dropdownWidth, 22];
            app.TaskTypeDropDown.Value = 'Непрерывная модель, точечные ограничения';

            % Create MLabel
            app.MLabel = uilabel(app.UIFigure);
            app.MLabel.HorizontalAlignment = 'right';
            app.MLabel.Position = [leftMargin, topMargin - rowHeight*1, labelWidth, 22];
            app.MLabel.Text = 'M (Размер рынка)';

            % Create MEditField
            app.MEditField = uieditfield(app.UIFigure, 'numeric');
            app.MEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*1, fieldWidth, 22];
            app.MEditField.Value = 50000;

            % Create pLabel
            app.pLabel = uilabel(app.UIFigure);
            app.pLabel.HorizontalAlignment = 'right';
            app.pLabel.Position = [leftMargin, topMargin - rowHeight*2, labelWidth, 22];
            app.pLabel.Text = 'p (Цена за единицу)';

            % Create pEditField
            app.pEditField = uieditfield(app.UIFigure, 'numeric');
            app.pEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*2, fieldWidth, 22];
            app.pEditField.Value = 25;

            % Create cLabel
            app.cLabel = uilabel(app.UIFigure);
            app.cLabel.HorizontalAlignment = 'right';
            app.cLabel.Position = [leftMargin, topMargin - rowHeight*3, labelWidth, 22];
            app.cLabel.Text = 'c (Переменные затраты)';

            % Create cEditField
            app.cEditField = uieditfield(app.UIFigure, 'numeric');
            app.cEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*3, fieldWidth, 22];
            app.cEditField.Value = 15;

            % Create rLabel
            app.rLabel = uilabel(app.UIFigure);
            app.rLabel.HorizontalAlignment = 'right';
            app.rLabel.Position = [leftMargin, topMargin - rowHeight*4, labelWidth, 22];
            app.rLabel.Text = 'r (Ставка дисконтир.)';

            % Create rEditField
            app.rEditField = uieditfield(app.UIFigure, 'numeric');
            app.rEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*4, fieldWidth, 22];
            app.rEditField.Value = 0.06;

            % Create deltaLabel
            app.deltaLabel = uilabel(app.UIFigure);
            app.deltaLabel.HorizontalAlignment = 'right';
            app.deltaLabel.Position = [leftMargin, topMargin - rowHeight*5, labelWidth, 22];
            app.deltaLabel.Text = 'delta (Коэф. затухания)';

            % Create deltaEditField
            app.deltaEditField = uieditfield(app.UIFigure, 'numeric');
            app.deltaEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*5, fieldWidth, 22];
            app.deltaEditField.Value = 0.7;

            % Create alphaLabel
            app.alphaLabel = uilabel(app.UIFigure);
            app.alphaLabel.HorizontalAlignment = 'right';
            app.alphaLabel.Position = [leftMargin, topMargin - rowHeight*6, labelWidth, 22];
            app.alphaLabel.Text = 'alpha (Эффект. рекламы)';

            % Create alphaEditField
            app.alphaEditField = uieditfield(app.UIFigure, 'numeric');
            app.alphaEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*6, fieldWidth, 22];
            app.alphaEditField.Value = 1.5;

            % Create ULabel
            app.ULabel = uilabel(app.UIFigure);
            app.ULabel.HorizontalAlignment = 'right';
            app.ULabel.Position = [leftMargin, topMargin - rowHeight*7, labelWidth, 22];
            app.ULabel.Text = 'U (Бюджет на рекламу)';

            % Create UEditField
            app.UEditField = uieditfield(app.UIFigure, 'numeric');
            app.UEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*7, fieldWidth, 22];
            app.UEditField.Value = 200000;

            % Create s0Label
            app.s0Label = uilabel(app.UIFigure);
            app.s0Label.HorizontalAlignment = 'right';
            app.s0Label.Position = [leftMargin, topMargin - rowHeight*8, labelWidth, 22];
            app.s0Label.Text = 'x0 (Начальная доля)';

            % Create s0EditField
            app.s0EditField = uieditfield(app.UIFigure, 'numeric');
            app.s0EditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*8, fieldWidth, 22];
            app.s0EditField.Value = 0;

            % Create TLabel
            app.TLabel = uilabel(app.UIFigure);
            app.TLabel.HorizontalAlignment = 'right';
            app.TLabel.Position = [leftMargin, topMargin - rowHeight*9, labelWidth, 22];
            app.TLabel.Text = 'T (Горизонт планирования)';

            % Create TEditField
            app.TEditField = uieditfield(app.UIFigure, 'numeric');
            app.TEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*9, fieldWidth, 22];
            app.TEditField.Value = 15;

            % Create NLabel
            app.NLabel = uilabel(app.UIFigure);
            app.NLabel.HorizontalAlignment = 'right';
            app.NLabel.Position = [leftMargin, topMargin - rowHeight*10, labelWidth, 22];
            app.NLabel.Text = 'N (Временные шаги)';

            % Create NEditField
            app.NEditField = uieditfield(app.UIFigure, 'numeric');
            app.NEditField.Position = [leftMargin + labelWidth + 10, topMargin - rowHeight*10, fieldWidth, 22];
            app.NEditField.Value = 150;

            % Create taLabel
            app.taLabel = uilabel(app.UIFigure);
            app.taLabel.HorizontalAlignment = 'right';
            app.taLabel.Position = [leftMargin, topMargin - rowHeight*11, 100, 22];
            app.taLabel.Text = 't_a';

            % Create taEditField
            app.taEditField = uieditfield(app.UIFigure, 'text');
            app.taEditField.Position = [leftMargin + 110, topMargin - rowHeight*11, 340, 22];
            app.taEditField.Value = '[1 5 10]';

            % Create xminLabel
            app.xminLabel = uilabel(app.UIFigure);
            app.xminLabel.HorizontalAlignment = 'right';
            app.xminLabel.Position = [leftMargin, topMargin - rowHeight*12, 100, 22];
            app.xminLabel.Text = 'x_min';

            % Create xminEditField
            app.xminEditField = uieditfield(app.UIFigure, 'text');
            app.xminEditField.Position = [leftMargin + 110, topMargin - rowHeight*12, 340, 22];
            app.xminEditField.Value = '[5 15 25]';

            % Create ConfirmButton
            app.ConfirmButton = uibutton(app.UIFigure, 'push');
            app.ConfirmButton.ButtonPushedFcn = createCallbackFcn(app, @ConfirmButtonPushed, true);
            app.ConfirmButton.Position = [leftMargin, topMargin - rowHeight*13, 200, 22];
            app.ConfirmButton.Text = 'Подтвердить';

            % Create SaveParamsButton (no functionality)
            app.SaveParamsButton = uibutton(app.UIFigure, 'push');
            app.SaveParamsButton.Position = [leftMargin + 200 + buttonGap, topMargin - rowHeight*13, buttonWidth, buttonHeight];
            app.SaveParamsButton.Text = 'Сохр. параметры';

            % Create LoadParamsButton (no functionality)
            app.LoadParamsButton = uibutton(app.UIFigure, 'push');
            app.LoadParamsButton.Position = [leftMargin + 200 + buttonGap*2 + buttonWidth, topMargin - rowHeight*13, buttonWidth, buttonHeight];
            app.LoadParamsButton.Text = 'Загр. параметры';

            % Create SaveArraysButton (no functionality)
            app.SaveArraysButton = uibutton(app.UIFigure, 'push');
            app.SaveArraysButton.Position = [leftMargin + 200 + buttonGap*3 + buttonWidth*2, topMargin - rowHeight*13, buttonWidth, buttonHeight];
            app.SaveArraysButton.Text = 'Сохр. массивы';

            % Create MessageLabel with adjusted position and larger font
            app.MessageLabel = uilabel(app.UIFigure);
            app.MessageLabel.Position = [leftMargin, topMargin - rowHeight*13.5, 1100, 22];
            app.MessageLabel.Text = '';
            app.MessageLabel.FontWeight = 'bold';
            app.MessageLabel.FontSize = 16;

            % Create PlotAxes (First plot)
            app.PlotAxes = uiaxes(app.UIFigure);
            app.PlotAxes.Position = [600, 400, 550, 350];
            app.PlotAxes.Title.String = 'Оптимальный объем продаж';
            app.PlotAxes.XLabel.String = 't';
            app.PlotAxes.YLabel.String = 'x(t)';

            % Create PlotAxes2 (Second plot)
            app.PlotAxes2 = uiaxes(app.UIFigure);
            app.PlotAxes2.Position = [600, 50, 550, 350];
            app.PlotAxes2.Title.String = 'Оптимальные рекламные расходы';
            app.PlotAxes2.XLabel.String = 't';
            app.PlotAxes2.YLabel.String = 'u^2(t)';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ModelParametersApp
            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)
        end

        % Code that executes before app deletion
        function delete(app)
            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end