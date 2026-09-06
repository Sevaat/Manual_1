% ============================================================
% ПРИМЕР 10.3: Сценарный анализ и критерий N-1
% Пятиузловая сеть 220/110/10 кВ
% ============================================================
clear; clc; close all;
fprintf('=== ПРИМЕР 10.3: Анализ надёжности (N-1) ===\n\n');

S_base = 100;  % МВ·А
n_bus = 5;

% === 1. ИСХОДНЫЕ ДАННЫЕ СЕТИ ===
% Ветви: [from, to, R, X, k, type]
% type: 1 = линия, 2 = трансформатор
branch_data = [
  1, 2, 0.01,  0.06, 1,  1;   % линия 1-2
  1, 3, 0.015, 0.08, 1,  1;   % линия 1-3
  2, 3, 0.012, 0.07, 1,  1;   % линия 2-3
  3, 4, 0.005, 0.04, 2,  2;   % АТ 3-4, k=2
  4, 5, 0.008, 0.06, 11, 2;   % Т 4-5, k=11
];

% Нагрузки
P_load = [0; 0; 0.60; 0.40; 0.25];
Q_load = [0; 0; 0.30; 0.20; 0.12];

% === 2. ФУНКЦИЯ ФОРМИРОВАНИЯ YBUS ===
function Ybus = make_ybus(n_bus, br_data)
  Ybus = zeros(n_bus);
  for m = 1:size(br_data, 1)
    i = br_data(m,1); j = br_data(m,2);
    R = br_data(m,3); X = br_data(m,4);
    k = br_data(m,5);
    y_s = 1/(R + 1i*X);
    if br_data(m,6) == 1  % линия
      Ybus(i,i) = Ybus(i,i) + y_s;
      Ybus(j,j) = Ybus(j,j) + y_s;
      Ybus(i,j) = Ybus(i,j) - y_s;
      Ybus(j,i) = Ybus(j,i) - y_s;
    else  % трансформатор
      Ybus(i,i) = Ybus(i,i) + y_s;
      Ybus(j,j) = Ybus(j,j) + y_s * k^2;
      Ybus(i,j) = Ybus(i,j) - y_s * k;
      Ybus(j,i) = Ybus(j,i) - y_s * k;
    end
  end
end

% === 3. ФОРМИРОВАНИЕ СЦЕНАРИЕВ ===
n_sc = 0;
scenarios = struct();

% базовый режим
n_sc = n_sc + 1;
scenarios(n_sc).name = 'Базовый режим';
scenarios(n_sc).scale = 1.0;
scenarios(n_sc).outage = 0;

% +20% нагрузки
n_sc = n_sc + 1;
scenarios(n_sc).name = 'Нагрузка +20%';
scenarios(n_sc).scale = 1.2;
scenarios(n_sc).outage = 0;

% -30% нагрузки
n_sc = n_sc + 1;
scenarios(n_sc).name = 'Нагрузка -30%';
scenarios(n_sc).scale = 0.7;
scenarios(n_sc).outage = 0;

% N-1: отключение каждой ветви
branch_names = {'Линия 1-2', 'Линия 1-3', 'Линия 2-3', 'АТ 3-4', 'Тр. 4-5'};
for l = 1:5
  n_sc = n_sc + 1;
  scenarios(n_sc).name = sprintf('N-1: %s', branch_names{l});
  scenarios(n_sc).scale = 1.0;
  scenarios(n_sc).outage = l;
end

% === 4. ВЫПОЛНЕНИЕ РАСЧЁТОВ ===
fprintf('Выполнение %d сценариев...\n', n_sc);
fprintf('%-25s %4s %4s %7s %7s %7s %12s\n', ...
  'Сценарий', 'Сх.', 'Ит.', 'U_min', 'U_max', 'dP МВт', 'Статус');
fprintf('%s\n', repmat('-', 1, 75));

results = struct();
for s = 1:n_sc
  sc = scenarios(s);

  % масштабирование нагрузки
  P_load_s = P_load * sc.scale;
  Q_load_s = Q_load * sc.scale;

  % модификация топологии
  if sc.outage > 0
    br_mod = branch_data([1:sc.outage-1, sc.outage+1:end], :);
  else
    br_mod = branch_data;
  end

  % формирование Ybus
  Ybus_s = make_ybus(n_bus, br_mod);

  % описание узлов
  buses_s(1).type='SLACK'; buses_s(1).P=0; buses_s(1).Q=0;
  buses_s(1).U_set=1.0; buses_s(1).Q_min=-5; buses_s(1).Q_max=5;
  buses_s(2).type='PV'; buses_s(2).P=-0.8; buses_s(2).Q=0;
  buses_s(2).U_set=1.03; buses_s(2).Q_min=-0.5; buses_s(2).Q_max=0.8;
  for i = 3:5
    buses_s(i).type='PQ';
    buses_s(i).P=P_load_s(i); buses_s(i).Q=Q_load_s(i);
    buses_s(i).U_set=1.0; buses_s(i).Q_min=0; buses_s(i).Q_max=0;
  end

  % расчёт
  [U_s, iter_s, conv_s, ~] = power_flow_newton(Ybus_s, buses_s, 1e-8, 100);

  % результаты
  results(s).name = sc.name;
  results(s).conv = conv_s;
  results(s).iter = iter_s;
  results(s).U = U_s;
  results(s).U_min = min(abs(U_s));
  results(s).U_max = max(abs(U_s));

  if conv_s
    % потери (упрощённо: через баланс)
    S_sl = U_s(1) * conj(Ybus_s(1,:) * U_s);
    P_gen_total = real(S_sl) + 0.8;  % Slack + PV
    P_load_total = sum(P_load_s(3:5));
    results(s).P_loss = P_gen_total - P_load_total;

    % проверка напряжений
    U_viol = any(abs(U_s) < 0.95 | abs(U_s) > 1.05);
    if ~U_viol
      results(s).status = 'ДОПУСТИМ';
    else
      results(s).status = 'НАРУШЕНИЕ';
    end
  else
    results(s).P_loss = NaN;
    results(s).status = 'РАСХОДИМОСТЬ';
  end

  fprintf('%-25s %4d %4d %7.4f %7.4f %7.2f %12s\n', ...
    sc.name, conv_s, iter_s, results(s).U_min, results(s).U_max, ...
    results(s).P_loss*S_base, results(s).status);
end

% === 5. СВОДНАЯ ТАБЛИЦА ===
fprintf('\n\n=== СВОДНАЯ ТАБЛИЦА ===\n');
fprintf('%-25s %8s %8s %8s %12s\n', ...
  'Сценарий', 'U_min', 'U_max', 'dP МВт', 'Статус');
fprintf('%s\n', repmat('=', 1, 65));
for s = 1:n_sc
  fprintf('%-25s %8.4f %8.4f %8.2f %12s\n', ...
    results(s).name, results(s).U_min, results(s).U_max, ...
    results(s).P_loss*S_base, results(s).status);
end

% === 6. ВИЗУАЛИЗАЦИЯ ===
figure('Position', [50 50 1100 500]);

subplot(1,2,1);
hold on;
colors = {'b', 'r', 'g', 'm', 'c', 'k'};
for s = 1:min(n_sc, 6)
  if results(s).conv
    plot(1:5, abs(results(s).U), [colors{mod(s-1,6)+1} '-o'], ...
      'LineWidth', 1.2, 'MarkerSize', 5, 'DisplayName', results(s).name);
  end
end
yline(0.95, 'r--'); yline(1.05, 'r--');
xlabel('Узел'); ylabel('U, о.е.');
title('Профили напряжений по сценариям');
ylim([0.85 1.15]); grid on;
legend('Location', 'southwest', 'FontSize', 7);

subplot(1,2,2);
P_losses_all = zeros(1, n_sc);
for s = 1:n_sc
  P_losses_all(s) = results(s).P_loss * S_base;
end
bar(P_losses_all);
set(gca, 'XTickLabel', {results.name}, 'XTickLabelRotation', 45, 'FontSize', 7);
ylabel('dP, МВт'); title('Потери по сценариям');
grid on;

sgtitle('ПРИМЕР 10.3: Сценарный анализ критерия N-1');
fprintf('\n=== АНАЛИЗ ЗАВЕРШЁН ===\n');
