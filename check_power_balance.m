% Проверка баланса активной и реактивной мощностей
% на вход:
%   U        — вектор напряжений
%   buses    — структура узлов
%   branches — массив ветвей
% на выход:
%   balance — структура с результатами проверки
function balance = check_power_balance(U, buses, branches)
  n = length(buses);

  % генерация
  P_gen = 0; Q_gen = 0;
  for i = 1:n
    if strcmp(buses(i).type, 'SLACK') || strcmp(buses(i).type, 'PV')
      P_gen = P_gen + buses(i).P;
      Q_gen = Q_gen + buses(i).Q;
    end
  end

  % нагрузка
  P_load = 0; Q_load = 0;
  for i = 1:n
    if strcmp(buses(i).type, 'PQ')
      P_load = P_load + buses(i).P;
      Q_load = Q_load + buses(i).Q;
    end
  end

  % потери
  [~, dP, dQ] = calc_branch_flows(U, branches);

  % баланс (выражения 9.50, 9.52)
  balance.P_gen = -P_gen;
  balance.P_load = P_load;
  balance.P_loss = dP;
  balance.P_mismatch = -P_gen - P_load - dP;

  balance.Q_gen = -Q_gen;
  balance.Q_load = Q_load;
  balance.Q_loss = dQ;
  balance.Q_mismatch = -Q_gen - Q_load - dQ;

  fprintf('=== Баланс мощностей ===\n');
  fprintf('P: генерация = %.4f, нагрузка = %.4f, потери = %.4f, невязка = %.2e\n', ...
    balance.P_gen, balance.P_load, balance.P_loss, balance.P_mismatch);
  fprintf('Q: генерация = %.4f, нагрузка = %.4f, потери = %.4f, невязка = %.2e\n', ...
    balance.Q_gen, balance.Q_load, balance.Q_loss, balance.Q_mismatch);
end
