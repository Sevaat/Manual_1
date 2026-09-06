% Пример расчёта режима методом Ньютона-Рафсона
% Сеть из 3 узлов

% === Формирование Ybus ===
Z12 = 0.02 + 1i*0.08;
Z13 = 0.03 + 1i*0.10;
Z23 = 0.025 + 1i*0.09;

y12 = 1/Z12; y13 = 1/Z13; y23 = 1/Z23;

Ybus = [ y12+y13,  -y12,     -y13;
         -y12,      y12+y23, -y23;
         -y13,     -y23,      y13+y23];

% === Описание узлов ===
buses(1).type = 'SLACK'; buses(1).P = 0;    buses(1).Q = 0;
buses(1).U_set = 1.05;   buses(1).Q_min = -5; buses(1).Q_max = 5;

buses(2).type = 'PV';    buses(2).P = -0.5; buses(2).Q = 0;
buses(2).U_set = 1.02;   buses(2).Q_min = -0.4; buses(2).Q_max = 0.6;

buses(3).type = 'PQ';    buses(3).P = 0.8;  buses(3).Q = 0.4;
buses(3).U_set = 1.0;    buses(3).Q_min = 0; buses(3).Q_max = 0;

% === Расчёт ===
tol = 1e-8;
k_max = 50;
[U, iter, conv, mismatch_hist] = power_flow_newton(Ybus, buses, tol, k_max);

% === Вывод результатов ===
fprintf('Сходимость: %d, итераций: %d\n', conv, iter);
fprintf('\nНапряжения узлов:\n');
for i = 1:3
  fprintf('  Узел %d: U = %.6f о.е., угол = %.4f град\n', ...
          i, abs(U(i)), angle(U(i))*180/pi);
end

% === История сходимости ===
fprintf('\nИстория невязок:\n');
for k = 1:length(mismatch_hist)
  fprintf('  Итерация %d: макс. невязка = %.2e\n', k, mismatch_hist(k));
end
