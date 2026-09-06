% ============================================================
% ГЛАВНЫЙ СКРИПТ: Полный цикл автоматизированного расчёта
% Чтение данных -> Расчёт -> Анализ -> Отчёт
% ============================================================
clear; clc; close all;

fprintf('============================================================\n');
fprintf('  АВТОМАТИЗИРОВАННЫЙ РАСЧЁТ УСТАНОВИВШЕГОСЯ РЕЖИМА\n');
fprintf('  Дата: %s\n', datestr(now, 'dd.mm.yyyy HH:MM:SS'));
fprintf('============================================================\n\n');

% === 1. ЧТЕНИЕ ИСХОДНЫХ ДАННЫХ ===
fprintf('[1/6] Чтение исходных данных...\n');

bus_raw = dlmread('data/buses.csv', ',', 1, 0);
n_bus = size(bus_raw, 1);
for i = 1:n_bus
  buses(i).id = bus_raw(i,1);
  switch bus_raw(i,2)
    case 1, buses(i).type = 'SLACK';
    case 2, buses(i).type = 'PV';
    case 3, buses(i).type = 'PQ';
  end
  buses(i).U_set  = bus_raw(i,3);
  buses(i).P      = bus_raw(i,4) - bus_raw(i,6);  % нагрузка - генерация
  buses(i).Q      = bus_raw(i,5);
  buses(i).Q_min  = bus_raw(i,7);
  buses(i).Q_max  = bus_raw(i,8);
end

br_raw = dlmread('data/branches.csv', ',', 1, 0);
n_br = size(br_raw, 1);
for l = 1:n_br
  if br_raw(l,3) == 1
    branches(l).type = 'L';
  else
    branches(l).type = 'T2';
  end
  branches(l).from = br_raw(l,1);
  branches(l).to   = br_raw(l,2);
  branches(l).R    = br_raw(l,4);
  branches(l).X    = br_raw(l,5);
  branches(l).G    = br_raw(l,6);
  branches(l).B    = br_raw(l,7);
  branches(l).U1   = br_raw(l,8);
  branches(l).U2   = br_raw(l,9);
end

load_raw = dlmread('data/loads.csv', ',', 1, 0);
n_loads = size(load_raw, 1);
for l = 1:n_loads
  loads(l).bus = load_raw(l,1);
  loads(l).P0  = load_raw(l,2);
  loads(l).Q0  = load_raw(l,3);
  loads(l).zip = load_raw(l,4:9);
  loads(l).U0  = load_raw(l,10);
end

fprintf('  Загружено: %d узлов, %d ветвей, %d нагрузок\n', ...
  n_bus, n_br, n_loads);

% === 2. ФОРМИРОВАНИЕ Y_BUS ===
fprintf('[2/6] Формирование матрицы Y_bus...\n');
Y_shunt = zeros(n_bus, 1);
% добавить шунты из данных узлов
for i = 1:n_bus
  Y_shunt(i) = 1i * bus_raw(i,9);
end
Ybus = build_ybus(n_bus, branches, Y_shunt);
fprintf('  Матрица %d×%d, ненулевых: %d (%.1f%%)\n', ...
  n_bus, n_bus, nnz(Ybus), nnz(Ybus)/n_bus^2*100);

% === 3. РАСЧЁТ РЕЖИМА ===
fprintf('[3/6] Расчёт режима методом Ньютона-Рафсона...\n');
tic;
[U, iter, conv, hist] = power_flow_newton(Ybus, buses, 1e-10, 100);
t_calc = toc;

if conv
  fprintf('  Сходимость за %d итераций (%.4f с)\n', iter, t_calc);
else
  fprintf('  ОШИБКА: расчёт не сошёлся!\n');
  return;
end

% === 4. ПОСТ-ОБРАБОТКА ===
fprintf('[4/6] Пост-обработка...\n');
[flows, P_loss, Q_loss] = calc_branch_flows(U, branches);

% проверки
U_min = min(abs(U));
U_max = max(abs(U));
U_violation = find(abs(U) < 0.95 | abs(U) > 1.05);

fprintf('  U_min = %.4f, U_max = %.4f о.е.\n', U_min, U_max);
fprintf('  Потери: dP = %.3f о.е., dQ = %.3f о.е.\n', P_loss, Q_loss);
if ~isempty(U_violation)
  fprintf('  ВНИМАНИЕ: нарушение в узлах: ');
  fprintf('%d ', U_violation); fprintf('\n');
else
  fprintf('  Все напряжения в допустимых пределах.\n');
end

% === 5. ФОРМИРОВАНИЕ ОТЧЁТА ===
fprintf('[5/6] Формирование отчёта...\n');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
report_file = sprintf('results/report_%s.txt', timestamp);

fid = fopen(report_file, 'w');
fprintf(fid, '============================================================\n');
fprintf(fid, '  ОТЧЁТ О РАСЧЁТЕ УСТАНОВИВШЕГОСЯ РЕЖИМА\n');
fprintf(fid, '  Дата: %s\n', datestr(now, 'dd.mm.yyyy HH:MM:SS'));
fprintf(fid, '============================================================\n\n');

fprintf(fid, '--- РЕЗУЛЬТАТЫ ПО УЗЛАМ ---\n');
fprintf(fid, '%5s %7s %8s %8s %10s %10s\n', ...
  'Узел', 'Тип', 'U о.е.', 'Угол°', 'P МВт', 'Q Мвар');
for i = 1:n_bus
  S_i = U(i) * conj(Ybus(i,:) * U);
  fprintf(fid, '%5d %7s %8.5f %8.3f %10.2f %10.2f\n', ...
    i, buses(i).type, abs(U(i)), angle(U(i))*180/pi, ...
    -real(S_i)*100, -imag(S_i)*100);
end

fprintf(fid, '\n--- РЕЗУЛЬТАТЫ ПО ВЕТВЯМ ---\n');
fprintf(fid, '%8s %10s %10s %8s %8s\n', ...
  'Ветвь', 'P_нач МВт', 'P_кон МВт', 'dP МВт', 'dQ Мвар');
for l = 1:n_br
  fprintf(fid, '%4d-%-3d %10.3f %10.3f %8.4f %8.4f\n', ...
    branches(l).from, branches(l).to, ...
    real(flows(l).S_from)*100, real(flows(l).S_to)*100, ...
    real(flows(l).dS)*100, imag(flows(l).dS)*100);
end

fprintf(fid, '\n--- БАЛАНС МОЩНОСТЕЙ ---\n');
fprintf(fid, 'Потери: dP = %.2f МВт, dQ = %.2f Мвар\n', ...
  P_loss*100, Q_loss*100);

fclose(fid);
fprintf('  Отчёт сохранён: %s\n', report_file);

% === 6. ВИЗУАЛИЗАЦИЯ ===
fprintf('[6/6] Визуализация...\n');
figure('Position', [100 100 800 600]);

subplot(2,2,1);
plot(1:n_bus, abs(U), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on; yline(1.0,'k:'); yline(0.95,'r--'); yline(1.05,'r--');
xlabel('Узел'); ylabel('U, о.е.'); title('Профиль напряжений');
grid on; ylim([0.9 1.1]);

subplot(2,2,2);
semilogy(hist, 'r-o', 'LineWidth', 1.5);
xlabel('Итерация'); ylabel('Невязка'); title('Сходимость');
grid on;

subplot(2,2,3);
bar(real(flows(1:n_br).S_from)*100);
ylabel('P, МВт'); title('Потоки активной мощности');
set(gca, 'XTickLabel', arrayfun(@(l) sprintf('%d-%d', ...
  branches(l).from, branches(l).to), 1:n_br, 'UniformOutput', false));
grid on;

subplot(2,2,4);
pie([sum([buses.P])*100, P_loss*100], {'Нагрузка', 'Потери'});
title('Распределение мощности');

sgtitle('Автоматизированный расчёт установившегося режима');
saveas(gcf, sprintf('results/plot_%s.png', timestamp));

fprintf('\n============================================================\n');
fprintf('  РАСЧЁТ УСПЕШНО ЗАВЕРШЁН\n');
fprintf('============================================================\n');
