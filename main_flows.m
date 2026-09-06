% Расчёт потоков мощностей для сети из 3 узлов
% (параметры из раздела 8.4.3.7)
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

% расчёт режима
[U, iter, conv] = power_flow_newton(Ybus, buses, 1e-8, 50);

% формирование массива ветвей
branches(1).type='L'; branches(1).from=1; branches(1).to=2;
branches(1).R=0.02; branches(1).X=0.08; branches(1).G=0; branches(1).B=0;
branches(2).type='L'; branches(2).from=1; branches(2).to=3;
branches(2).R=0.03; branches(2).X=0.10; branches(2).G=0; branches(2).B=0;
branches(3).type='L'; branches(3).from=2; branches(3).to=3;
branches(3).R=0.025; branches(3).X=0.09; branches(3).G=0; branches(3).B=0;

% расчёт потоков
[flows, P_loss, Q_loss] = calc_branch_flows(U, branches);

% вывод результатов
fprintf('Потоки мощностей по ветвям:\n');
for l = 1:length(branches)
  fprintf('Ветвь %d-%d: S_нач = %.4f + j%.4f, S_кон = %.4f + j%.4f, ', ...
    branches(l).from, branches(l).to, ...
    real(flows(l).S_from), imag(flows(l).S_from), ...
    real(flows(l).S_to), imag(flows(l).S_to));
  fprintf('dS = %.4f + j%.4f\n', real(flows(l).dS), imag(flows(l).dS));
end
fprintf('Суммарные потери: dP = %.4f о.е., dQ = %.4f о.е.\n', P_loss, Q_loss);
