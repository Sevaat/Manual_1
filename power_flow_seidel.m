% Расчёт установившегося режима методом Гаусса-Зейделя
% на вход передаются:
%   Ybus   — матрица узловых проводимостей, о.е.
%   buses  — структура узлов
%   U_init — начальное приближение
%   tol    — допуск сходимости
%   k_max  — макс. итераций
%   alpha  — коэффициент верхней релаксации (по умолчанию 1.4,
%            рекомендуемый диапазон 1.2–1.8, см. формулу (8.30))
% на выход:
%   U, iter, conv — аналогично методу Якоби
function [U, iter, conv] = power_flow_seidel(Ybus, buses, U_init, tol, k_max, alpha)
  if nargin < 6
    alpha = 1.4;  % без релаксации
  end
  n = length(buses);
  U = U_init(:);
  conv = 0;

  % балансирующий узел определяется по типу, а не по номеру
  idx_SL = find(strcmp({buses.type}, 'SLACK'));

  for iter = 1:k_max
  U_prev = U;   % для критерия сходимости
  % определение типов узлов (обновляется на каждой итерации
  % для учёта возможного переключения PV -> PQ)
  idx_PQ = find(strcmp({buses.type}, 'PQ'));
  idx_PV = find(strcmp({buses.type}, 'PV'));
    for i = 1:n   % обход всех узлов, балансирующий пропускается
      if any(i == idx_SL)
        continue; % Slack-узел: напряжение фиксировано, не обновляется
      end
      if any(i == idx_PQ)
        % PQ-узел
        S_spec = -(buses(i).P + 1i * buses(i).Q);
        sum_term = 0;
        for j = 1:n
          if j ~= i
            sum_term = sum_term + Ybus(i, j) * U(j); % U(j) уже обновлён для j < i
          end
        end
        U_new = (conj(S_spec) / conj(U(i)) - sum_term) / Ybus(i, i);
        % релаксация
        U(i) = (1 - alpha) * U(i) + alpha * U_new;

      elseif any(i == idx_PV)
        % PV-узел: вычисляем Q
        I_inj = Ybus(i, :) * U;
        S_calc = U(i) * conj(I_inj);
        Q_calc = -imag(S_calc);

        % проверка ограничений
        Q_min = buses(i).Q_min;
        Q_max = buses(i).Q_max;
        if Q_calc < -Q_max
          buses(i).type = 'PQ';
          buses(i).Q = -Q_max;
          % пересчитать как PQ
          S_spec = -(buses(i).P + 1i * buses(i).Q);
          sum_term = 0;
          for j = 1:n
            if j ~= i
              sum_term = sum_term + Ybus(i, j) * U(j);
            end
          end
          U_new = (conj(S_spec) / conj(U(i)) - sum_term) / Ybus(i, i);
          U(i) = (1 - alpha) * U(i) + alpha * U_new;
        elseif Q_calc > -Q_min
          buses(i).type = 'PQ';
          buses(i).Q = -Q_min;
          S_spec = -(buses(i).P + 1i * buses(i).Q);
          sum_term = 0;
          for j = 1:n
            if j ~= i
              sum_term = sum_term + Ybus(i, j) * U(j);
            end
          end
          U_new = (conj(S_spec) / conj(U(i)) - sum_term) / Ybus(i, i);
          U(i) = (1 - alpha) * U(i) + alpha * U_new;
        else
          % Q в пределах: поддерживаем модуль напряжения
          U(i) = buses(i).U_set * exp(1i * angle(U(i)));
        end
      end
    end

    % проверка сходимости
    dU = max(abs(U - U_prev));
    if dU < tol
      conv = 1;
      return;
    end
  end
end
