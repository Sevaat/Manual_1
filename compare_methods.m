% Сравнение методов расчёта установившегося режима
% Сеть из 3 узлов (см. main_newton.m)

% === Исходные данные (аналогично разделу 8.4.3.7) ===
Z12 = 0.02 + 1i*0.08;
Z13 = 0.03 + 1i*0.10;
Z23 = 0.025 + 1i*0.09;
y12 = 1/Z12; y13 = 1/Z13; y23 = 1/Z23;
Ybus = [ y12+y13,  -y12,     -y13;
         -y12,      y12+y23, -y23;
         -y13,     -y23,      y13+y23];

buses(1).type='SLACK'; buses(1).P=0; buses(1).Q=0;
buses(1).U_set=1.05; buses(1).Q_min=-5; buses(1).Q_max=5;
buses(2).type='PV'; buses(2).P=-0.5; buses(2).Q=0;
buses(2).U_set=1.02; buses(2).Q_min=-0.4; buses(2).Q_max=0.6;
buses(3).type='PQ'; buses(3).P=0.8; buses(3).Q=0.4;
buses(3).U_set=1.0; buses(3).Q_min=0; buses(3).Q_max=0;

tol = 1e-6;
k_max = 500;
U_init = [1.05; 1.02; 1.0];  % начальное приближение

% === Метод Ньютона-Рафсона ===
tic;
[U_NR, iter_NR, conv_NR, hist_NR] = power_flow_newton(Ybus, buses, tol, k_max, U_init);
t_NR = toc;

% === Метод Гаусса-Зейделя ===
tic;
[U_GS, iter_GS, conv_GS] = power_flow_seidel(Ybus, buses, U_init, tol, k_max);
t_GS = toc;

% === Метод Якоби ===
tic;
[U_J, iter_J, conv_J] = power_flow_jacobi(Ybus, buses, U_init, tol, k_max);
t_J = toc;

% === Вывод ===
fprintf('=== Результаты сравнения ===\n');
fprintf('Метод Ньютона:   итераций = %d, время = %.4f с, сходимость = %d\n', ...
        iter_NR, t_NR, conv_NR);
fprintf('Метод Зейделя:   итераций = %d, время = %.4f с, сходимость = %d\n', ...
        iter_GS, t_GS, conv_GS);
fprintf('Метод Якоби:     итераций = %d, время = %.4f с, сходимость = %d\n', ...
        iter_J, t_J, conv_J);

% === Проверка совпадения результатов ===
fprintf('\nМакс. разность напряжений:\n');
fprintf('  Ньютон vs Зейдель: %.2e\n', max(abs(U_NR - U_GS)));
fprintf('  Ньютон vs Якоби:   %.2e\n', max(abs(U_NR - U_J)));

% Построение графика сходимости
figure;
semilogy(hist_NR, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
% Для Зейделя и Якоби необходимо модифицировать функции
% для возврата истории невязок
xlabel('Номер итерации');
ylabel('Максимальная невязка, о.е.');
title('Сходимость метода Ньютона-Рафсона');
grid on;
