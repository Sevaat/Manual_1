% Метод продолжения по параметру (гомотопия)
% Функция предполагает, что power_flow_newton.m модифицирована
% для приёма начального приближения U_init (см. раздел 8.5.1)
%
% на вход передаются:
%   Ybus     — матрица узловых проводимостей, о.е.
%   buses    — структура узлов при полной нагрузке
%   tol      — допуск сходимости
%   k_max    — макс. итераций для каждого шага
%   d_lambda — шаг изменения параметра (по умолчанию 0.1)
% на выход:
%   U_final    — напряжения при максимальной достигнутой нагрузке
%   conv_all   — 1 если все шаги сошлись, 0 если расходимость
%   lambda_max — максимальное lambda, при котором расчёт сошёлся
function [U_final, conv_all, lambda_max] = power_flow_continuation( ...
    Ybus, buses, tol, k_max, d_lambda)

  if nargin < 5
    d_lambda = 0.1;
  end

  n = length(buses);
  conv_all = 1;
  lambda_max = 0;

  % сохраняем исходные заданные мощности
  P_orig = [buses.P];
  Q_orig = [buses.Q];

  % начальное приближение: плоский старт
  U_init = ones(n, 1);
  for i = 1:n
    if strcmp(buses(i).type,'SLACK') || strcmp(buses(i).type,'PV')
      U_init(i) = buses(i).U_set;
    end
  end

  % последовательность значений параметра нагрузки
  lambdas = d_lambda : d_lambda : 1.0;

  for idx = 1:length(lambdas)
    lambda = lambdas(idx);

    % масштабируем нагрузку и генерацию
    buses_lambda = buses;
    for i = 1:n
    	if strcmp(buses(i).type, 'PQ')
	        % масштабируем только нагрузку (PQ-узлы)
	        buses_lambda(i).P = lambda * P_orig(i);
	        buses_lambda(i).Q = lambda * Q_orig(i);
    	else
	        % генерация PV и SLACK не масштабируется
	        buses_lambda(i).P = P_orig(i);
	        buses_lambda(i).Q = Q_orig(i);
    	end
	end

    % расчёт с использованием решения предыдущего шага
    [U, iter, conv] = power_flow_newton( ...
        Ybus, buses_lambda, tol, k_max, U_init);

    if conv
      lambda_max = lambda;
      U_init = U;  % решение как начальное приближение
    else
      conv_all = 0;
      fprintf('Расходимость при lambda = %.2f\n', lambda);
      break;
    end
  end

  U_final = U_init;
end
