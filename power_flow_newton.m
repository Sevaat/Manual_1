% Расчёт установившегося режима методом Ньютона-Рафсона
% в полярных координатах
%
% Соглашение о знаках:
%   нагрузка (потребление)  > 0
%   генерация               < 0
%   buses(i).P, buses(i).Q  — заданные ПОТРЕБЛЯЕМЫЕ мощности
%   buses(i).Q_max          — максимальная ГЕНЕРАЦИЯ (> 0)
%   buses(i).Q_min          — минимальная ГЕНЕРАЦИЯ (< 0, т.е. потребление)
%
% на вход передаются:
%   Ybus   — матрица узловых проводимостей, о.е.
%   buses  — структура узлов с полями:
%            .type   — 'SLACK', 'PV', 'PQ'
%            .P      — заданная активная мощность, о.е. (нагрузка > 0)
%            .Q      — заданная реактивная мощность, о.е. (нагрузка > 0)
%            .U_set  — заданный модуль напряжения, о.е. (для PV, SLACK)
%            .Q_min, .Q_max — пределы по Q генерации (для PV)
%   tol    — допуск сходимости
%   k_max  — макс. число итераций
%   U_init — (необязательно) начальное приближение напряжений
%
% на выход:
%   U             — вектор комплексных напряжений, о.е.
%   iter          — число итераций
%   conv          — флаг сходимости
%   mismatch_hist — история макс. невязки для каждого шага

function [U, iter, conv, mismatch_hist] = power_flow_newton(Ybus, buses, tol, k_max, U_init)

  n = length(buses);
  conv = 0;
  mismatch_hist = [];

  % === Начальное приближение ===
  if nargin >= 5 && ~isempty(U_init)
    U     = abs(U_init(:));
    delta = angle(U_init(:));
  else
    % плоский старт
    U     = ones(n, 1);
    delta = zeros(n, 1);
  end

  % фиксация модулей для SLACK и PV
  for i = 1:n
    if strcmp(buses(i).type, 'SLACK') || strcmp(buses(i).type, 'PV')
      U(i) = buses(i).U_set;
    end
  end

  % === Индексы узлов по типам ===
  idx_PQ   = find(strcmp({buses.type}, 'PQ'));
  idx_PV   = find(strcmp({buses.type}, 'PV'));
  idx_PQPV = [idx_PQ, idx_PV];   % узлы с уравнением по P

  G = real(Ybus);
  B = imag(Ybus);

  for iter = 1:k_max

    % =========================================================
    % 1. РАСЧЁТ МОЩНОСТЕЙ при текущих напряжениях
    %    P_calc, Q_calc — расчётные ПОТРЕБЛЯЕМЫЕ мощности
    %    (положительны для нагрузки, отрицательны для генерации)
    % =========================================================
    P_calc = zeros(n, 1);
    Q_calc = zeros(n, 1);
    for i = 1:n
      for j = 1:n
        P_calc(i) = P_calc(i) - U(i)*U(j) * ...
            (G(i,j)*cos(delta(i)-delta(j)) + B(i,j)*sin(delta(i)-delta(j)));
        Q_calc(i) = Q_calc(i) - U(i)*U(j) * ...
            (G(i,j)*sin(delta(i)-delta(j)) - B(i,j)*cos(delta(i)-delta(j)));
      end
    end

    % =========================================================
    % 2. ПРОВЕРКА ОГРАНИЧЕНИЙ PV-УЗЛОВ
    %    Выполняется ДО вычисления невязок, чтобы при переключении
    %    узла из PV в PQ невязки формировались уже с актуальными
    %    индексами и новым заданным Q.
    %
    %    Логика проверки знаков:
    %      Q_calc(i)        — расчётное ПОТРЕБЛЕНИЕ (нагрузка > 0)
    %      генерация        = -Q_calc(i)
    %      верхний предел   : -Q_calc > Q_max  <=>  Q_calc < -Q_max
    %      нижний предел    : -Q_calc < Q_min  <=>  Q_calc > -Q_min
    % =========================================================
    for i = idx_PV
      if Q_calc(i) < -buses(i).Q_max
        % генерация превысила верхний предел:
        % фиксируем генерацию на Q_max => потребление = -Q_max
        buses(i).type = 'PQ';
        buses(i).Q    = -buses(i).Q_max;
        idx_PQ   = [idx_PQ, i];
        idx_PV   = idx_PV(idx_PV ~= i);
        idx_PQPV = [idx_PQ, idx_PV];

      elseif Q_calc(i) > -buses(i).Q_min
        % генерация ниже нижнего предела:
        % фиксируем генерацию на Q_min => потребление = -Q_min
        buses(i).type = 'PQ';
        buses(i).Q    = -buses(i).Q_min;
        idx_PQ   = [idx_PQ, i];
        idx_PV   = idx_PV(idx_PV ~= i);
        idx_PQPV = [idx_PQ, idx_PV];
      end
    end

    % =========================================================
    % 3. НЕВЯЗКИ
    %    Вычисляются ПОСЛЕ возможного переключения типа узла,
    %    поэтому для переключённого узла dQ формируется уже
    %    с новым заданным Q, а не остаётся нулевой.
    % =========================================================
    dP = zeros(n, 1);
    dQ = zeros(n, 1);
    for i = idx_PQPV
      dP(i) = buses(i).P - P_calc(i);   % нагрузка положительна
    end
    for i = idx_PQ
      dQ(i) = buses(i).Q - Q_calc(i);
    end

    % =========================================================
    % 4. ПРОВЕРКА СХОДИМОСТИ
    % =========================================================
    max_mismatch  = max([abs(dP(idx_PQPV)); abs(dQ(idx_PQ))]);
    mismatch_hist = [mismatch_hist; max_mismatch];

    if max_mismatch < tol
      conv = 1;
      break;
    end

    % =========================================================
    % 5. ФОРМИРОВАНИЕ МАТРИЦЫ ЯКОБИ
    %    Размерности вычисляются после возможного переключения,
    %    поэтому корректны для обновлённых наборов узлов.
    % =========================================================
    nP = length(idx_PQPV);   % число уравнений по P
    nQ = length(idx_PQ);     % число уравнений по Q
    nJ = nP + nQ;            % размер Якоби

    J = zeros(nJ, nJ);

    % --- Блок H = d(dP)/d(delta) и N = d(dP)/dU * U ---
    for ii = 1:nP
      i = idx_PQPV(ii);

      % столбцы по углам (блоки H)
      for jj = 1:nP
        j = idx_PQPV(jj);
        if i == j
          J(ii, jj) = Q_calc(i);                               % H_ii
        else
          J(ii, jj) = U(i)*U(j) * ...
              (G(i,j)*sin(delta(i)-delta(j)) - B(i,j)*cos(delta(i)-delta(j))); % H_ij
        end
      end

      % столбцы по модулям напряжений (блоки N)
      for jj = 1:nQ
        j = idx_PQ(jj);
        if i == j
          J(ii, nP+jj) = -P_calc(i) + G(i,i)*U(i)^2;          % N_ii
        else
          J(ii, nP+jj) = U(i)*U(j) * ...
              (G(i,j)*cos(delta(i)-delta(j)) + B(i,j)*sin(delta(i)-delta(j))); % N_ij
        end
      end
    end

    % --- Блок M = d(dQ)/d(delta) и L = d(dQ)/dU * U ---
    for ii = 1:nQ
      i = idx_PQ(ii);

      % столбцы по углам (блоки M)
      for jj = 1:nP
        j = idx_PQPV(jj);
        if i == j
          J(nP+ii, jj) = -P_calc(i);                           % M_ii
        else
          J(nP+ii, jj) = -U(i)*U(j) * ...
              (G(i,j)*cos(delta(i)-delta(j)) + B(i,j)*sin(delta(i)-delta(j))); % M_ij
        end
      end

      % столбцы по модулям напряжений (блоки L)
      for jj = 1:nQ
        j = idx_PQ(jj);
        if i == j
          J(nP+ii, nP+jj) = -Q_calc(i) - B(i,i)*U(i)^2;       % L_ii
        else
          J(nP+ii, nP+jj) = U(i)*U(j) * ...
              (G(i,j)*sin(delta(i)-delta(j)) - B(i,j)*cos(delta(i)-delta(j))); % L_ij
        end
      end
    end

    % =========================================================
    % 6. РЕШЕНИЕ СИСТЕМЫ ПОПРАВОК
    %    Элементы блока N определены как (d(dP)/dU)*U,
    %    поэтому поправка содержит нормированную величину dU/U.
    % =========================================================
    F  = [dP(idx_PQPV); dQ(idx_PQ)];
    dx = J \ (-F);

    % =========================================================
    % 7. ОБНОВЛЕНИЕ НЕИЗВЕСТНЫХ
    % =========================================================
    % обновление углов
    for ii = 1:nP
      i = idx_PQPV(ii);
      delta(i) = delta(i) + dx(ii);
    end

    % обновление модулей напряжений для PQ-узлов
    % dx(nP + ii) содержит нормированную поправку dU/U
    for ii = 1:nQ
      i = idx_PQ(ii);
      U(i) = U(i) + dx(nP + ii) * U(i);
    end

    % фиксация модулей для PV и SLACK
    idx_PV_fix  = find(strcmp({buses.type}, 'PV'));
    idx_SL_fix  = find(strcmp({buses.type}, 'SLACK'));
    for i = [idx_PV_fix, idx_SL_fix]
      U(i) = buses(i).U_set;
    end

  end   % конец цикла итераций

  % === Формирование комплексных напряжений ===
  U = U .* exp(1i * delta);

end
