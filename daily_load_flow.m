% Многорежимный расчёт по суточному графику нагрузки
% на вход:
%   Ybus      — матрица проводимостей
%   buses     — структура узлов при максимальной нагрузке
%   branches  — массив ветвей
%   k_profile — вектор из 24 значений относительной нагрузки (0..1)
% на выход:
%   P_loss_hourly — вектор потерь для каждого часа
%   W_loss_daily  — суточные потери энергии, о.е.*ч
%   U_min_hourly  — минимальные напряжения по часам
function [P_loss_hourly, W_loss_daily, U_min_hourly] = daily_load_flow(Ybus, buses, branches, k_profile)
  n_hours = length(k_profile);
  P_loss_hourly = zeros(1, n_hours);
  U_min_hourly = zeros(1, n_hours);
  W_loss_daily = 0;
  U_init = [];

  for h = 1:n_hours
    k = k_profile(h);
    buses_h = buses;
    for i = 1:length(buses)
      if strcmp(buses(i).type, 'PQ')
        buses_h(i).P = k * buses(i).P;
        buses_h(i).Q = k * buses(i).Q;
      end
    end

    if isempty(U_init)
      [U, ~, conv] = power_flow_newton(Ybus, buses_h, 1e-8, 50);
    else
      [U, ~, conv] = power_flow_newton(Ybus, buses_h, 1e-8, 50, U_init);
    end

    if conv
      [~, dP, ~] = calc_branch_flows(U, branches);
      P_loss_hourly(h) = dP;
      U_min_hourly(h) = min(abs(U));
      W_loss_daily = W_loss_daily + dP * 1;  % dt = 1 ч
      U_init = U;
    else
      fprintf('Расходимость в час %d (k = %.2f)\n', h-1, k);
      P_loss_hourly(h) = NaN;
      U_min_hourly(h) = NaN;
    end
  end

  % визуализация
  figure;
  subplot(2,1,1);
  plot(0:23, P_loss_hourly, 'r-o', 'LineWidth', 1.5);
  xlabel('Час суток'); ylabel('Потери, о.е.');
  title('Суточные потери активной мощности');
  grid on;

  subplot(2,1,2);
  plot(0:23, U_min_hourly, 'b-o', 'LineWidth', 1.5);
  xlabel('Час суток'); ylabel('U_{min}, о.е.');
  title('Минимальное напряжение по часам');
  grid on;

  fprintf('Суточные потери энергии: %.4f о.е.*ч (%.2f МВт*ч при S_б=100 МВ·А)\n', ...
    W_loss_daily, W_loss_daily * 100);
end
