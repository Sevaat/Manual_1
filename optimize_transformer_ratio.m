% Оптимизация коэффициента трансформации
% для минимизации суммарных потерь
% на вход:
%   buses    — структура узлов
%   branches — массив ветвей (трансформатор имеет поле U1, U2)
%   tr_idx   — индекс ветви-трансформатора в массиве branches
%   k_range  — диапазон коэффициентов трансформации [k_min, k_max]
%   n_steps  — число шагов перебора
% на выход:
%   k_opt    — оптимальный коэффициент трансформации
%   P_min    — минимальные потери
function [k_opt, P_min] = optimize_transformer_ratio(buses, branches, tr_idx, k_range, n_steps)
  k_values = linspace(k_range(1), k_range(2), n_steps);
  P_losses = zeros(1, n_steps);

  n_bus = length(buses);
  Y_shunt = zeros(n_bus, 1);

  for s = 1:n_steps
    % модификация коэффициента трансформации
    branches_mod = branches;
    branches_mod(tr_idx).U2 = branches(tr_idx).U1 / k_values(s);

    % формирование Ybus
    Ybus = build_ybus(n_bus, branches_mod, Y_shunt);

    % расчёт режима
    [U, ~, conv] = power_flow_newton(Ybus, buses, 1e-8, 50);

    if conv
      [~, dP, ~] = calc_branch_flows(U, branches_mod);
      P_losses(s) = dP;
    else
      P_losses(s) = Inf;  % расходимость — недопустимый режим
    end
  end

  % поиск минимума
  [P_min, idx_min] = min(P_losses);
  k_opt = k_values(idx_min);

  fprintf('Оптимальный k_тр = %.4f, потери dP = %.6f о.е.\n', k_opt, P_min);
end
