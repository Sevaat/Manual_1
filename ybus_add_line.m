% Учёт линии электропередачи в матрице Y_bus
% Выражения (7.17)–(7.20)
% на вход передаются:
% Ybus — текущая матрица проводимостей
% br   — структура с полями: from, to, R, X, G, B
% на выход: обновлённая матрица Ybus
function Ybus = ybus_add_line(Ybus, br)
  i = br.from;
  j = br.to;

  % продольная проводимость (7.17)
  y_series = 1 / (br.R + 1i * br.X);

  % поперечная проводимость
  y_shunt = br.G + 1i * br.B;

  % внедиагональные элементы (7.18)
  Ybus(i, j) = Ybus(i, j) - y_series;
  Ybus(j, i) = Ybus(j, i) - y_series;

  % диагональные элементы (7.19), (7.20)
  Ybus(i, i) = Ybus(i, i) + y_series + 0.5 * y_shunt;
  Ybus(j, j) = Ybus(j, j) + y_series + 0.5 * y_shunt;
end
