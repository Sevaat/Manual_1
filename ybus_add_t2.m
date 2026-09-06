% Учёт двухобмоточного трансформатора в матрице Y_bus
% Выражения (7.21)–(7.24)
% на вход передаются:
% Ybus — текущая матрица проводимостей
% br   — структура с полями: from (ВН), to (НН), R, X, G, B, U1, U2
% на выход: обновлённая матрица Ybus
function Ybus = ybus_add_t2(Ybus, br)
  i = br.from;  % узел ВН
  j = br.to;    % узел НН

  % продольная проводимость (7.21)
  y_series = 1 / (br.R + 1i * br.X);

  % поперечная проводимость
  y_shunt = br.G + 1i * br.B;

  % коэффициент трансформации
  k = br.U1 / br.U2;

  % внедиагональные элементы (7.22)
  Ybus(i, j) = Ybus(i, j) - y_series * k;
  Ybus(j, i) = Ybus(j, i) - y_series * k;

  % диагональные элементы (7.23), (7.24)
  Ybus(i, i) = Ybus(i, i) + y_series + y_shunt;
  Ybus(j, j) = Ybus(j, j) + y_series * k^2;
end
