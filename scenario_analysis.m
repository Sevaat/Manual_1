% Сценарный анализ: серия расчётов с различными условиями
% на вход передаются:
%   Ybus_base  — базовая матрица проводимостей
%   buses_base — базовая структура узлов
%   branches   — массив ветвей
%   scenarios  — массив структур сценариев:
%                .name — имя сценария
%                .scale_P — коэффициент масштабирования нагрузки
%                .scale_Q — коэффициент масштабирования нагрузки
%                .outage  — номер отключаемой ветви (0 = нет)
% на выход:
%   results — структура с результатами для каждого сценария
function results = scenario_analysis(Ybus_base, buses_base, branches, scenarios)
  n_sc = length(scenarios);
  results = struct();

  for s = 1:n_sc
    sc = scenarios(s);
    fprintf('Сценарий %d: %s\n', s, sc.name);

    % модификация нагрузки
    buses = buses_base;
    for i = 1:length(buses)
      if strcmp(buses(i).type, 'PQ')
        buses(i).P = buses_base(i).P * sc.scale_P;
        buses(i).Q = buses_base(i).Q * sc.scale_Q;
      end
    end

    % модификация топологии (отключение ветви)
    if sc.outage > 0 && sc.outage <= length(branches)
      branches_mod = branches([1:sc.outage-1, sc.outage+1:end]);
    else
      branches_mod = branches;
    end

    % перестроение Ybus (при отключении ветви)
    if sc.outage > 0
      n_bus = length(buses);
      Y_shunt = zeros(n_bus, 1);
      Ybus = build_ybus(n_bus, branches_mod, Y_shunt);
    else
      Ybus = Ybus_base;
    end

    % расчёт режима
    [U, iter, conv] = power_flow_newton(Ybus, buses, 1e-8, 50);

    % сохранение результатов
    results(s).name = sc.name;
    results(s).conv = conv;
    results(s).iter = iter;
    results(s).U = U;
    results(s).U_min = min(abs(U));
    results(s).U_max = max(abs(U));

    if conv
      [flows, dP, dQ] = calc_branch_flows(U, branches_mod);
      results(s).P_loss = dP;
      results(s).Q_loss = dQ;
      results(s).flows = flows;
    else
      results(s).P_loss = NaN;
      results(s).Q_loss = NaN;
    end

    fprintf('  Сходимость: %d, итераций: %d, U_min = %.4f, U_max = %.4f\n', ...
      conv, iter, results(s).U_min, results(s).U_max);
  end
end
