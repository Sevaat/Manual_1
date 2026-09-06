% ============================================================
% ПРИМЕР 10.2: Сеть с трансформатором, ZIP-нагрузкой и КБ
% ============================================================
clear; clc; close all;
fprintf('=== ПРИМЕР 10.2: Сеть 110/10 кВ с ZIP и КБ ===\n\n');

S_base = 100;       % МВ·А
U_base1 = 110;      % кВ (сторона ВН)
U_base2 = 10.5;     % кВ (сторона НН)
Z_base1 = U_base1^2 / S_base;  % 121 Ом
Z_base2 = U_base2^2 / S_base;  % 1.1025 Ом

% === 1. ПАРАМЕТРЫ ЛИНИИ 1-2 (таблица 4.1) ===
R1 = 3.5 / Z_base1;
X1 = 12.266 / Z_base1;
B1 = 3.3316e-4 * Z_base1;

fprintf('Линия 1-2 (о.е.): R=%.5f, X=%.5f, B=%.5f\n', R1, X1, B1);

% === 2. ПАРАМЕТРЫ ТРАНСФОРМАТОРА 2-3 (таблица 4.2) ===
R2 = 0.3546 / Z_base1;
X2 = 9.6016 / Z_base1;
G2 = 9.5622e-6 * Z_base1;
B2 = -6.5569e-5 * Z_base1;
U1_tr = 121;   % кВ (ВН)
U2_tr = 10.5;  % кВ (НН)

fprintf('Трансформатор 2-3 (о.е.): R=%.6f, X=%.6f\n', R2, X2);
fprintf('  G_хх=%.3e, B_хх=%.3e, k=%.3f\n', G2, B2, U1_tr/U2_tr);

% === 3. ФОРМИРОВАНИЕ Y_BUS ===
n_bus = 3;
branches(1).type = 'L'; branches(1).from = 1; branches(1).to = 2;
branches(1).R = R1; branches(1).X = X1;
branches(1).G = 0; branches(1).B = B1;

branches(2).type = 'T2'; branches(2).from = 2; branches(2).to = 3;
branches(2).R = R2; branches(2).X = X2;
branches(2).G = G2; branches(2).B = B2;
branches(2).U1 = U1_tr; branches(2).U2 = U2_tr;

% КБ в узле 3: Q_КБ = 8 Мвар
Q_kb = 8;  % Мвар (генерация, >0 по модулю)
B_kb = Q_kb / S_base;  % о.е. (ёмкостная, >0)

% Шунты узлов
Y_shunt = zeros(n_bus, 1);
Y_shunt(3) = 1i * B_kb;  % КБ в узле 3

Ybus = build_ybus(n_bus, branches, Y_shunt);
fprintf('\nY_bus сформирована (размер %d×%d)\n', n_bus, n_bus);

% === 4. НАГРУЗКИ ===
% Узел 2: простая PQ-нагрузка
P2 = 10 / S_base;
Q2 = P2 * tan(acos(0.9));

% Узел 3: ZIP-нагрузка (промышленная, таблица 6.1)
P3_0 = 50 / S_base;
Q3_0 = P3_0 * tan(acos(0.85));
zip3 = [0.22, 0.02, 0.76, 0.15, 0.09, 0.76];

% проверка нормировки
if ~check_zip_coefficients(zip3)
  error('ZIP-коэффициенты не нормированы!');
end
fprintf('ZIP-коэффициенты узла 3 проверены (нормировка ОК)\n');

% Для первого расчёта: модель постоянной мощности (U=1.0)
buses(1).type='SLACK'; buses(1).P=0; buses(1).Q=0;
buses(1).U_set=1.05; buses(1).Q_min=-5; buses(1).Q_max=5;
buses(2).type='PQ'; buses(2).P=P2; buses(2).Q=Q2;
buses(2).U_set=1.0; buses(2).Q_min=0; buses(2).Q_max=0;
buses(3).type='PQ'; buses(3).P=P3_0; buses(3).Q=Q3_0;
buses(3).U_set=1.0; buses(3).Q_min=0; buses(3).Q_max=0;

% === 5. РАСЧЁТ (этап 1: постоянная мощность) ===
fprintf('\n--- Этап 1: модель постоянной мощности ---\n');
[U1, iter1, conv1, hist1] = power_flow_newton(Ybus, buses, 1e-10, 50);
fprintf('Сходимость: %d, итераций: %d\n', conv1, iter1);

% === 6. УТОЧНЕНИЕ НАГРУЗКИ ПО ZIP ===
U3_mag = abs(U1(3));
[P3_zip, Q3_zip] = get_zip_power(P3_0, Q3_0, U3_mag, 1.0, zip3);

fprintf('\nУточнение по ZIP (U_3 = %.4f о.е.):\n', U3_mag);
fprintf('  P: %.4f -> %.4f о.е. (изм. %.2f%%)\n', ...
  P3_0, P3_zip, (P3_zip/P3_0 - 1)*100);
fprintf('  Q: %.4f -> %.4f о.е. (изм. %.2f%%)\n', ...
  Q3_0, Q3_zip, (Q3_zip/Q3_0 - 1)*100);

% === 7. РАСЧЁТ (этап 2: уточнённая ZIP-нагрузка) ===
buses(3).P = P3_zip;
buses(3).Q = Q3_zip;

fprintf('\n--- Этап 2: уточнённая ZIP-нагрузка ---\n');
[U2, iter2, conv2, hist2] = power_flow_newton(Ybus, buses, 1e-10, 50, U1);
fprintf('Сходимость: %d, итераций: %d (нач. приближение из этапа 1)\n', ...
  conv2, iter2);

% === 8. ИТОГОВЫЕ РЕЗУЛЬТАТЫ ===
U_nom = [110; 110; 10.5];
fprintf('\n=== ИТОГОВЫЕ РЕЗУЛЬТАТЫ ===\n');
fprintf('%5s %7s %8s %9s %8s\n', 'Узел', 'Тип', 'U о.е.', 'U кВ', 'Угол°');
for i = 1:3
  fprintf('%5d %7s %8.5f %9.2f %8.3f\n', ...
    i, buses(i).type, abs(U2(i)), abs(U2(i))*U_nom(i), angle(U2(i))*180/pi);
end

% потоки и потери
[flows, P_loss, Q_loss] = calc_branch_flows(U2, branches);

fprintf('\nПотоки и потери:\n');
fprintf('%8s %12s %12s %8s %8s\n', 'Ветвь', 'P_нач МВт', 'P_кон МВт', 'dP МВт', 'dQ Мвар');
for l = 1:2
  fprintf('%4d-%-3d %12.2f %12.2f %8.3f %8.3f\n', ...
    branches(l).from, branches(l).to, ...
    real(flows(l).S_from)*S_base, real(flows(l).S_to)*S_base, ...
    real(flows(l).dS)*S_base, imag(flows(l).dS)*S_base);
end
fprintf('Суммарные потери: dP = %.2f МВт (%.2f%%), dQ = %.2f Мвар\n', ...
  P_loss*S_base, P_loss*S_base/(10+50)*100, Q_loss*S_base);

% баланс
S_slack = U2(1) * conj(Ybus(1,:) * U2);
fprintf('\nБаланс:\n');
fprintf('  Генерация: P = %.2f МВт, Q = %.2f Мвар\n', ...
  real(S_slack)*S_base, imag(S_slack)*S_base);
fprintf('  Нагрузка:  P = %.2f МВт, Q = %.2f Мвар\n', ...
  (P2+P3_zip)*S_base, (Q2+Q3_zip)*S_base);
fprintf('  Потери:    P = %.2f МВт\n', P_loss*S_base);
fprintf('  Невязка:   dP = %.2e\n', ...
  real(S_slack) - (P2+P3_zip) - P_loss);

% === 9. ВИЗУАЛИЗАЦИЯ ===
figure('Position', [100 100 1000 500]);

subplot(1,2,1);
bar([abs(U2(1)), abs(U2(2)), abs(U2(3))], 0.5);
hold on; yline(1.0,'k--'); yline(0.95,'r--'); yline(1.05,'r--');
set(gca, 'XTickLabel', {'1 (110 кВ)', '2 (110 кВ)', '3 (10 кВ)'});
ylabel('U, о.е.'); title('Профиль напряжений');
ylim([0.9 1.1]); grid on;

subplot(1,2,2);
% сравнение нагрузки: номинал и ZIP
U_range = linspace(0.8, 1.2, 50);
P_zip_curve = zeros(size(U_range));
Q_zip_curve = zeros(size(U_range));
for k = 1:length(U_range)
  [P_zip_curve(k), Q_zip_curve(k)] = get_zip_power(P3_0, Q3_0, U_range(k), 1.0, zip3);
end
plot(U_range, P_zip_curve*S_base, 'b-', 'LineWidth', 2); hold on;
yline(P3_0*S_base, 'b--', 'P = const');
xline(U3_mag, 'r--', 'U факт.');
xlabel('U, о.е.'); ylabel('P, МВт');
title('ZIP-характеристика нагрузки узла 3');
grid on; legend('ZIP', 'P=const', 'U факт.');

sgtitle('ПРИМЕР 10.2: Сеть с трансформатором и ZIP-нагрузкой');
fprintf('\n=== РАСЧЁТ ЗАВЕРШЁН ===\n');
