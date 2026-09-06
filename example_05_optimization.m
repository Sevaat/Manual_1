% ============================================================
% ПРИМЕР 10.5: Оптимизация коэффициента трансформации
% Минимизация суммарных потерь
% ============================================================
clear; clc; close all;
fprintf('=== ПРИМЕР 10.5: Оптимизация k_тр ===\n\n');

S_base = 100;
U_base1 = 110;
Z_base1 = U_base1^2 / S_base;

% === 1. ПАРАМЕТРЫ ЛИНИИ (фиксированные) ===
R1 = 3.5 / Z_base1;
X1 = 12.266 / Z_base1;
B1 = 3.3316e-4 * Z_base1;

% === 2. ПАРАМЕТРЫ ТРАНСФОРМАТОРА ===
R2 = 0.3546 / Z_base1;
X2 = 9.6016 / Z_base1;
G2 = 9.5622e-6 * Z_base1;
B2 = -6.5569e-5 * Z_base1;
U1_tr = 121;  % кВ (ВН, фиксировано)

% === 3. НАГРУЗКА ===
P_load = 50 / S_base;
Q_load = P_load * tan(acos(0.85));

% === 4. ПЕРЕБОР КОЭФФИЦИЕНТОВ ТРАНСФОРМАЦИИ ===
k_range = 10.5:0.25:12.5;
n_steps = length(k_range);
P_losses = zeros(1, n_steps);
U_load = zeros(1, n_steps);

fprintf('Перебор k_тр от %.2f до %.2f (%d значений)...\n', ...
  k_range(1), k_range(end), n_steps);
fprintf('%8s %10s %10s\n', 'k_тр', 'dP МВт', 'U_нагр');
fprintf('%s\n', repmat('-', 1, 30));

n_bus = 2;
for s = 1:n_steps
  k_tr = k_range(s);
  U2_tr = U1_tr / k_tr;  % напряжение НН из k

  % формирование ветвей
  branches(1).type = 'L'; branches(1).from = 1; branches(1).to = 2;
  branches(1).R = R1; branches(1).X = X1;
  branches(1).G = 0; branches(1).B = B1;

  branches(2).type = 'T2'; branches(2).from = 2; branches(2).to = 3;
  branches(2).R = R2; branches(2).X = X2;
  branches(2).G = G2; branches(2).B = B2;
  branches(2).U1 = U1_tr; branches(2).U2 = U2_tr;

  Y_shunt = zeros(3, 1);
  Ybus = build_ybus(3, branches, Y_shunt);

  % узлы
  buses(1).type='SLACK'; buses(1).P=0; buses(1).Q=0;
  buses(1).U_set=1.05; buses(1).Q_min=-5; buses(1).Q_max=5;
  buses(2).type='PQ'; buses(2).P=0; buses(2).Q=0;
  buses(2).U_set=1.0; buses(2).Q_min=0; buses(2).Q_max=0;
  buses(3).type='PQ'; buses(3).P=P_load; buses(3).Q=Q_load;
  buses(3).U_set=1.0; buses(3).Q_min=0; buses(3).Q_max=0;

  % расчёт
  [U, ~, conv] = power_flow_newton(Ybus, buses, 1e-8, 50);

  if conv
    [~, dP, ~] = calc_branch_flows(U, branches);
    P_losses(s) = dP;
    U_load(s) = abs(U(3));
  else
    P_losses(s) = Inf;
    U_load(s) = NaN;
  end

  fprintf('%8.2f %10.3f %10.4f\n', k_tr, P_losses(s)*S_base, U_load(s));
end

% === 5. ОПТИМУМ ===
[P_min, idx_opt] = min(P_losses);
k_opt = k_range(idx_opt);
fprintf('\n=== РЕЗУЛЬТАТ ОПТИМИЗАЦИИ ===\n');
fprintf('Оптимальный k_тр = %.2f (U_НН = %.2f кВ)\n', k_opt, U1_tr/k_opt);
fprintf('Минимальные потери: dP = %.3f МВт (%.2f%%)\n', ...
  P_min*S_base, P_min*S_base/(P_load*S_base)*100);
fprintf('Напряжение нагрузки: U = %.4f о.е.\n', U_load(idx_opt));

% === 6. ВИЗУАЛИЗАЦИЯ ===
figure('Position', [100 100 900 500]);

subplot(1,2,1);
plot(k_range, P_losses*S_base, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(k_opt, P_min*S_base, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('k_{тр}'); ylabel('dP, МВт');
title('Потери от коэффициента трансформации');
grid on; legend('dP(k)', sprintf('Оптимум k=%.2f', k_opt));

subplot(1,2,2);
plot(k_range, U_load, 'g-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on; yline(0.95, 'r--'); yline(1.05, 'r--');
xlabel('k_{тр}'); ylabel('U_{нагр}, о.е.');
title('Напряжение нагрузки от k_{тр}');
grid on;

sgtitle('ПРИМЕР 10.5: Оптимизация коэффициента трансформации');
fprintf('\n=== РАСЧЁТ ЗАВЕРШЁН ===\n');
