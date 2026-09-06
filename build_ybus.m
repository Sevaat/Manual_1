% Формирование матрицы узловых проводимостей Y_bus
% Выражения (7.15)–(7.40)
% на вход передаются:
% n_bus    — количество узлов сети
% branches — массив структур ветвей (линии, трансформаторы)
% Y_shunt  — вектор-столбец суммарных шунтирующих проводимостей узлов, См
% на выход: комплексная матрица Y_bus размера n_bus × n_bus
function Ybus = build_ybus(n_bus, branches, Y_shunt)
  % инициализация нулевой матрицы (7.15)
  Ybus = zeros(n_bus, n_bus);

  % учёт шунтов узлов (7.16)
  Ybus = Ybus + diag(Y_shunt);

  % обработка ветвей
  for m = 1:length(branches)
    br = branches(m);
    switch br.type
      case 'L'
        Ybus = ybus_add_line(Ybus, br);
      case 'T2'
        Ybus = ybus_add_t2(Ybus, br);
      case 'T3'
        Ybus = ybus_add_t3(Ybus, br);
      otherwise
        error('Неизвестный тип ветви: %s', br.type);
    end
  end
end
