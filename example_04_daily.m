% ============================================================
% ПРИМЕР 10.4: Многорежимный расчёт по суточному графику
% ============================================================
clear; clc; close all;
fprintf('=== ПРИМЕР 10.4: Суточный расчёт ===\n\n');

S_base = 100;  % МВ·А
U_base = 110;  % кВ
Z_base = U_base^2 / S_base;

% === 1. ПАРАМЕТРЫ СЕТИ (из примера 10.1) ===
R = 0.07713; X = 0.26347; B_pu = 0.02743;

n_bus = 2;
branches(1).type = 'L'; branches(1).from = 1; branches(1).to = 2;
branches(1).R = R; branches(1).X = X;
branches(1).G = 0; branches(1).B = B_pu;
Y_shunt = zeros(n_bus, 1);
Ybus = build_ybus(n_bus, branches, Y_shunt);

% === 2. СУТОЧНЫЙ ГРАФИК НАГРУЗКИ ===
% Типовой график промышленного района (24 часа)
k_profile = [0.65, 0.60, 0.58, 0.55, 0.57, 0.65, ...
             0.75, 0.85, 0.95, 1.00, 1.00, 0.98, ...
             0.95, 0.97, 1.00, 1.00, 0.95, 0.90, ...
             0.85, 0.80, 0.75, 0.72, 0.70, 0.68];

% Максимальная нагрузка
P_max = 40 / S_base;  % 40 МВт
cos_phi = 0.88;
Q_max = P_max * tan(acos(cos_phi));

fprintf('Максимальная нагрузка: P = %.0f МВт, Q = %.1f Мвар\n', ...
  P_max*S_base, Q_max*S_base);
fprintf('Число расчётных часов: %d\n\n', length(k_profile));

% === 3. МНОГОРЕЖИМНЫЙ РАСЧЁТ ===
n_hours = length(k_profile);
P_loss_hourly = zeros(1, n_hours);
U_min_hourly = zeros(1, n_hours);
W_loss_daily = 0;
U_init = [];

fprintf('%4s %6s %8s %8s %8s\n', 'Час', 'k', 'P МВт', 'dP МВт', 'U_min');
fprintf('%s\n', repmat('-', 1, 40));

for h = 1:n_hours
  k = k_profile(h);

  buses(1).type='SLACK'; buses(1).P=0; buses(1).Q=0;
  buses(1).U_set=1.05; buses(1).Q_min=-5; buses(1).Q_max=5;
  buses(2).type='PQ'; buses(2).P=k*P_max; buses(2).Q=k*Q_max;
  buses(2).U_set=1.0; buses(2).Q_min=0; buses(2).Q_max=0;

  if isempty(U_init)
    [U, ~, conv] = power_flow_newton(Ybus, buses, 1e-8, 50);
  else
    [U, ~, conv] = power_flow_newton(Ybus, buses, 1e-8, 50, U_init);
  end

  if conv
    [~, dP, ~] = calc_branch_flows(U, branches);
    P_loss_hourly(h) = dP;
    U_min_hourly(h) = abs(U(2));
    W_loss_daily = W_loss_daily + dP * 1;  % dt = 1 ч
    U_init = U;
  else
    fprintf('  Расходимость в час %d!\n', h-1);
    P_loss_hourly(h) = NaN;
    U_min_hourly(h) = NaN;
  end

  fprintf('%4d %6.2f %8.1f %8.3f %8.4f\n', ...
    h-1, k, k*P_max*S_base, P_loss_hourly(h)*S_base, U_min_hourly(h));
end

% === 4. ИТОГИ ===
fprintf('\n=== ИТОГИ СУТОЧНОГО РАСЧЁТА ===\n');
fprintf('Суточные потери энергии: %.4f о.е.*ч = %.2f МВт*ч\n', ...
  W_loss_daily, W_loss_daily * S_base);
fprintf('Максимальные потери: %.3f МВт (час %d)\n', ...
  max(P_loss_hourly)*S_base, find(P_loss_hourly==max(P_loss_hourly))-1);
fprintf('Минимальные потери: %.3f МВт (час %d)\n', ...
  min(P_loss_hourly)*S_base, find(P_loss_hourly==min(P_loss_hourly))-1);
fprintf('Минимальное напряжение: %.4f о.е. (час %d)\n', ...
  min(U_min_hourly), find(U_min_hourly==min(U_min_hourly))-1);
fprintf('Максимальное напряжение: %.4f о.е.\n', max(U_min_hourly));

% годовые потери (365 одинаковых суток)
W_annual = W_loss_daily * S_base * 365;
fprintf('\nГодовые потери энергии: %.1f МВт*ч = %.2f ГВт*ч\n', ...
  W_annual, W_annual/1000);

% === 5. ВИЗУАЛИЗАЦИЯ ===
figure('Position', [100 100 900 700]);

subplot(3,1,1);
plot(0:23, k_profile, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 5);
xlabel('Час суток'); ylabel('k(t)'); title('Суточный график нагрузки');
ylim([0 1.1]); grid on;

subplot(3,1,2);
plot(0:23, P_loss_hourly * S_base, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 5);
xlabel('Час суток'); ylabel('dP, МВт');
title('Суточные потери активной мощности');
grid on;

subplot(3,1,3);
plot(0:23, U_min_hourly, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 5);
hold on; yline(0.95, 'r--'); yline(1.05, 'r--');
xlabel('Час суток'); ylabel('U_2, о.е.');
title('Напряжение в узле нагрузки');
ylim([0.93 1.07]); grid on;

sgtitle('ПРИМЕР 10.4: Многорежимный расчёт');
fprintf('\n=== РАСЧЁТ ЗАВЕРШЁН ===\n');
