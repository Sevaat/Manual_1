% Построение профиля напряжений
% на вход:
%   U     — вектор комплексных напряжений, о.е.
%   U_nom — номинальные напряжения узлов, о.е. (обычно все = 1.0)
%   U_min, U_max — допустимые пределы
function plot_voltage_profile(U, U_nom, U_min, U_max)
  n = length(U);
  U_mag = abs(U);
  nodes = 1:n;

  figure;
  plot(nodes, U_mag, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
  hold on;
  plot(nodes, U_nom, 'k--', 'LineWidth', 1);
  plot(nodes, U_min * ones(n,1), 'r--', 'LineWidth', 1);
  plot(nodes, U_max * ones(n,1), 'r--', 'LineWidth', 1);

  % подсветка нарушений
  idx_viol = find(U_mag < U_min | U_mag > U_max);
  if ~isempty(idx_viol)
    plot(idx_viol, U_mag(idx_viol), 'rs', 'MarkerSize', 12, 'LineWidth', 2);
  end

  xlabel('Номер узла');
  ylabel('Напряжение, о.е.');
  title('Профиль напряжений');
  legend('U_i', 'U_{ном}', 'Допустимые пределы', 'Нарушения');
  grid on;
  ylim([0.9, 1.1]);
end
