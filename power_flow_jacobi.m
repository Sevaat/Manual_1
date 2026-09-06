% Расчёт установившегося режима методом простой итерации (Якоби)
% на вход передаются:
%   Ybus   — матрица узловых проводимостей, о.е.
%   buses  — структура узлов (поля: type, P, Q, U_set)
%   U_init — начальное приближение напряжений (комплексный вектор)
%   tol    — допуск сходимости
%   k_max  — максимальное число итераций
% на выход:
%   U     — вектор узловых напряжений
%   iter  — число выполненных итераций
%   conv  — флаг сходимости (1 — сошлось, 0 — нет)
function [U, iter, conv] = power_flow_jacobi(Ybus, buses, U_init, tol, k_max)
  n = length(buses);
  U = U_init(:);            % текущие напряжения
  U_old = U;                % напряжения предыдущей итерации
  conv = 0;

  % определение типов узлов
  idx_PQ = find(strcmp({buses.type}, 'PQ'));
  idx_PV = find(strcmp({buses.type}, 'PV'));
  idx_SL = find(strcmp({buses.type}, 'SLACK'));

  for iter = 1:k_max
    U_old = U;   % сохраняем значения предыдущей итерации
    
    % определение типов узлов (обновляется на каждой итерации
    % для учёта возможного переключения PV -> PQ)
    idx_PQ = find(strcmp({buses.type}, 'PQ'));
    idx_PV = find(strcmp({buses.type}, 'PV'));

    % --- обновление PQ-узлов ---
    for ii = idx_PQ
      S_spec = -(buses(ii).P + 1i * buses(ii).Q); % нагрузка положительна
      sum_term = 0;
      for jj = 1:n
        if jj ~= ii
          sum_term = sum_term + Ybus(ii, jj) * U_old(jj);
        end
      end
      U(ii) = (conj(S_spec) / conj(U_old(ii)) - sum_term) / Ybus(ii, ii);
    end

    % --- обновление PV-узлов (проверка ограничений по Q) ---
    for ii = idx_PV
      % вычисляем Q для проверки ограничений
      I_inj = Ybus(ii, :) * U_old;
      S_calc = U_old(ii) * conj(I_inj);
      Q_calc = -imag(S_calc); % знак: нагрузка положительна

      % проверка ограничений по Q
      if Q_calc < -buses(ii).Q_max
        % генерация = -Q_calc > Q_max: превышение верхнего предела
        % генерация превысила верхний предел: переводим в PQ
        buses(ii).type = 'PQ';
        buses(ii).Q = -buses(ii).Q_max; % потребление = -генерация
        S_spec = -(buses(ii).P + 1i * buses(ii).Q);
        sum_term = 0;
        for jj = 1:n
          if jj ~= ii
            sum_term = sum_term + Ybus(ii, jj) * U_old(jj);
          end
        end
        U(ii) = (conj(S_spec) / conj(U_old(ii)) - sum_term) / Ybus(ii, ii);

      elseif Q_calc > -buses(ii).Q_min
        % генерация ниже нижнего предела: переводим в PQ
        buses(ii).type = 'PQ';
        buses(ii).Q = -buses(ii).Q_min;
        S_spec = -(buses(ii).P + 1i * buses(ii).Q);
        sum_term = 0;
        for jj = 1:n
          if jj ~= ii
            sum_term = sum_term + Ybus(ii, jj) * U_old(jj);
          end
        end
        U(ii) = (conj(S_spec) / conj(U_old(ii)) - sum_term) / Ybus(ii, ii);

      else
        % Q в пределах: поддерживаем модуль напряжения
        U(ii) = buses(ii).U_set * exp(1i * angle(U_old(ii)));
      end
    end

    % --- проверка сходимости ---
    dU = max(abs(U(idx_PQ) - U_old(idx_PQ)));
    if dU < tol
      conv = 1;
      return;
    end
  end
end
