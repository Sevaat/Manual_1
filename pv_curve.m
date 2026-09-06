% Построение P-V кривой для оценки устойчивости по напряжению
% на вход:
%   Ybus       — матрица проводимостей
%   buses      — структура узлов
%   crit_bus   — номер критического узла
%   d_lambda   — шаг изменения нагрузки
%   lambda_max — максимальное значение параметра
% на выход:
%   U_curve  — массив напряжений критического узла
%   P_curve  — массив нагрузок
%   lambda_cr — критическое значение параметра
function [U_curve, P_curve, lambda_cr] = pv_curve(Ybus, buses, crit_bus, d_lambda, lambda_max)
  lambdas = 0:d_lambda:lambda_max;
  n_pts = length(lambdas);
  U_curve = zeros(1, n_pts);
  P_curve = zeros(1, n_pts);
  lambda_cr = 0;
  U_init = [];

  for s = 1:n_pts
    lambda = lambdas(s);
    buses_mod = buses;
    for i = 1:length(buses)
      if strcmp(buses(i).type, 'PQ')
        buses_mod(i).P = lambda * buses(i).P;
        buses_mod(i).Q = lambda * buses(i).Q;
      end
    end

    if isempty(U_init)
      [U, ~, conv] = power_flow_newton(Ybus, buses_mod, 1e-8, 50);
    else
      [U, ~, conv] = power_flow_newton(Ybus, buses_mod, 1e-8, 50, U_init);
    end

    if conv
      U_curve(s) = abs(U(crit_bus));
      P_curve(s) = lambda * buses(crit_bus).P;
      lambda_cr = lambda;
      U_init = U;
    else
      fprintf('Расходимость при lambda = %.2f\n', lambda);
      break;
    end
  end

  % визуализация
  figure;
  plot(P_curve(1:s), U_curve(1:s), 'b-', 'LineWidth', 2);
  xlabel('Активная мощность нагрузки, о.е.');
  ylabel('Напряжение критического узла, о.е.');
  title('P-V кривая');
  grid on;

  fprintf('Предел устойчивости: lambda = %.2f (запас %.1f%%)\n', ...
    lambda_cr, (lambda_cr - 1)*100);
end
