% Проверка загрузки ветвей по допустимому току
% на вход передаются:
%   flows    — структура потоков (из calc_branch_flows.m)
%   branches — массив структур ветвей
%   I_max    — вектор допустимых токов ветвей, о.е.
%   U        — вектор напряжений узлов, о.е.
% на выход:
%   loading  — вектор коэффициентов загрузки
%   status   — массив строк: 'Норма', 'Предупреждение', 'Перегрузка'
function [loading, status] = check_loading(flows, branches, I_max, U)
  n_br = length(branches);
  loading = zeros(n_br, 1);
  status = cell(n_br, 1);

  for l = 1:n_br
    % модуль тока через продольную ветвь (9.19)
    I_branch = abs(flows(l).I);

    % коэффициент загрузки (9.22)
    loading(l) = I_branch / I_max(l);

    % оценка состояния
    if loading(l) <= 1.0
      status{l} = 'Норма';
    elseif loading(l) <= 1.2
      status{l} = 'Предупреждение';
    else
      status{l} = 'Перегрузка';
    end
  end
end
