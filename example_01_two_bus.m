% ============================================================
% ПРИМЕР 10.1: Полный расчёт двухузловой сети 110 кВ
% От параметров провода до итогового отчёта
% ============================================================
clear; clc; close all;
fprintf('=== ПРИМЕР 10.1: Двухузловая сеть 110 кВ ===\n\n');

% === 1. БАЗИСНЫЕ ВЕЛИЧИНЫ ===
S_base = 100;         % МВ·А
U_base = 110;         % кВ
Z_base = U_base^2 / S_base;  % Ом
I_base = S_base / (sqrt(3) * U_base);  % кА

fprintf('Базисные величины:\n');
fprintf('  S_б = %d МВ·А, U_б = %d кВ\n', S_base, U_base);
fprintf('  Z_б = %.1f Ом, I_б = %.4f кА\n', Z_base, I_base);

% === 2. ПАРАМЕТРЫ ЛИНИИ (формулы главы 4) ===
rho = 28;           % Ом*мм^2/км
F_al = 240;         % мм^2
F_st = 32;          % мм^2
L = 80;             % км
n_c = 1;            % число цепей
n = 1;              % число проводов в фазе
D = 5000;           % мм (для 110 кВ)
r_pr = 10.8;        % мм (справочный радиус АС 240/32)

% активное сопротивление (4.2)
R_ohm = (1/n) * (1/n_c) * rho * L / F_al;
% погонное индуктивное сопротивление (4.8)
x0 = 0.144 * log10(D / r_pr) + 0.0157;
% индуктивное сопротивление (4.5)
X_ohm = x0 * L / n_c;
% погонная ёмкостная проводимость (4.12)
b0 = 7.58e-6 / log10(D / r_pr);
% ёмкостная проводимость (4.7)
B = b0 * L * n_c;

fprintf('\nПараметры линии (именованные):\n');
fprintf('  R = %.3f Ом, X = %.3f Ом (x0 = %.4f Ом/км)\n', R_ohm, X_ohm, x0);
fprintf('  B = %.4e См (b0 = %.4e См/км)\n', B, b0);

% перевод в о.е.
R = R_ohm / Z_base;
X = X_ohm / Z_base;
B_pu = B * Z_base;

fprintf('Параметры линии (о.е.):\n');
fprintf('  R = %.5f, X = %.5f, B = %.5f\n', R, X, B_pu);

% === 3. НАГРУЗКА ===
P_load_MW = 40;
cos_phi = 0.88;
P_load = P_load_MW / S_base;
Q_load = P_load * tan(acos(cos_phi));

fprintf('\nНагрузка: P = %.1f МВт (%.4f о.е.), Q = %.2f Мвар (%.4f о.е.)\n', ...
  P_load_MW, P_load, Q_load*S_base, Q_load);

% === 4. ФОРМИРОВАНИЕ Y_BUS ===
n_bus = 2;
branches(1).type = 'L'; branches(1).from = 1; branches(1).to = 2;
branches(1).R = R; branches(1).X = X;
branches(1).G = 0; branches(1).B = B_pu;

Y_shunt = zeros(n_bus, 1);
Ybus = build_ybus(n_bus, branches, Y_shunt);

fprintf('\nМатрица Y_bus:\n');
disp(Ybus);

% === 5. ОПИСАНИЕ УЗЛОВ ===
buses(1).type = 'SLACK'; buses(1).P = 0; buses(1).Q = 0;
buses(1).U_set = 1.05; buses(1).Q_min = -5; buses(1).Q_max = 5;
buses(2).type = 'PQ'; buses(2).P = P_load; buses(2).Q = Q_load;
buses(2).U_set = 1.0; buses(2).Q_min = 0; buses(2).Q_max = 0;

% === 6. РАСЧЁТ РЕЖИМА ===
fprintf('\n--- Расчёт режима методом Ньютона-Рафсона ---\n');
[U, iter, conv, hist] = power_flow_newton(Ybus, buses, 1e-10, 50);
fprintf('Сходимость: %d, итераций: %d\n', conv, iter);
for k = 1:length(hist)
  fprintf('  Итерация %d: невязка = %.4e\n', k, hist(k));
end

% === 7. РЕЗУЛЬТАТЫ ===
fprintf('\nНапряжения узлов:\n');
for i = 1:n_bus
  fprintf('  Узел %d: |U| = %.5f о.е. (%.2f кВ), угол = %.4f град\n', ...
    i, abs(U(i)), abs(U(i))*U_base, angle(U(i))*180/pi);
end

% === 8. ПОТОКИ И ПОТЕРИ ===
[flows, P_loss, Q_loss] = calc_branch_flows(U, branches);

fprintf('\nПоток по линии:\n');
fprintf('  Начало: S = %.4f + j%.4f о.е. (%.2f + j%.2f МВ·А)\n', ...
  real(flows(1).S_from), imag(flows(1).S_from), ...
  real(flows(1).S_from)*S_base, imag(flows(1).S_from)*S_base);
fprintf('  Конец:  S = %.4f + j%.4f о.е. (%.2f + j%.2f МВ·А)\n', ...
  real(flows(1).S_to), imag(flows(1).S_to), ...
  real(flows(1).S_to)*S_base, imag(flows(1).S_to)*S_base);
fprintf('  Потери: dP = %.4f о.е. (%.2f МВт), dQ = %.4f о.е. (%.2f Мвар)\n', ...
  P_loss, P_loss*S_base, Q_loss, Q_loss*S_base);

% === 9. ПРОВЕРКА БАЛАНСА ===
S_source = U(1) * conj(Ybus(1,:) * U);
fprintf('\nБаланс мощностей:\n');
fprintf('  Генерация: P = %.2f МВт, Q = %.2f Мвар\n', ...
  real(S_source)*S_base, imag(S_source)*S_base);
fprintf('  Нагрузка:  P = %.2f МВт, Q = %.2f Мвар\n', ...
  P_load*S_base, Q_load*S_base);
fprintf('  Потери:    P = %.2f МВт, Q = %.2f Мвар\n', ...
  P_loss*S_base, Q_loss*S_base);
fprintf('  Невязка:   dP = %.2e, dQ = %.2e\n', ...
  real(S_source) - P_load - P_loss, imag(S_source) - Q_load - Q_loss);

% === 10. ПРОВЕРКА ДОПУСТИМОСТИ ===
I_line = abs(flows(1).I);
I_max_ohm = 605;  % А, допустимый ток АС 240/32
I_max_pu = I_max_ohm / (I_base * 1000);
K_load = I_line / I_max_pu;

fprintf('\nПроверка допустимости:\n');
fprintf('  Ток линии: %.1f А (%.4f о.е.)\n', I_line*I_base*1000, I_line);
fprintf('  I_доп = %d А (%.4f о.е.)\n', I_max_ohm, I_max_pu);
fprintf('  K_загр = %.3f (%s)\n', K_load, ...
  {'Норма', 'Перегрузка'}{(K_load > 1) + 1});
fprintf('  U_2 = %.4f о.е. ', abs(U(2)));
if abs(U(2)) >= 0.95 && abs(U(2)) <= 1.05
  fprintf('(в пределах 0.95...1.05)\n');
else
  fprintf('(НАРУШЕНИЕ пределов!)\n');
end

% === 11. ВИЗУАЛИЗАЦИЯ ===
figure('Position', [100 100 900 600]);

subplot(2,2,1);
plot([1 2], abs(U), 'b-o', 'LineWidth', 2, 'MarkerSize', 10);
hold on;
yline(1.05, 'r--'); yline(0.95, 'r--'); yline(1.0, 'k:');
xlabel('Узел'); ylabel('U, о.е.'); title('Профиль напряжений');
set(gca, 'XTick', [1 2]); ylim([0.9 1.1]); grid on;

subplot(2,2,2);
semilogy(hist, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Итерация'); ylabel('Невязка, о.е.');
title('Сходимость Ньютона-Рафсона'); grid on;

subplot(2,2,3);
P_vals = [P_load*S_base, P_loss*S_base];
pie(P_vals, {'Нагрузка 40 МВт', 'Потери'});
title('Распределение активной мощности');

subplot(2,2,4);
compass(U);
title('Векторная диаграмма напряжений');

sgtitle('ПРИМЕР 10.1: Двухузловая сеть 110 кВ');
fprintf('\n=== РАСЧЁТ ЗАВЕРШЁН ===\n');
