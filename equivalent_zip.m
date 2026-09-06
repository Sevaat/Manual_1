%эквивалентирование ZIP-коэффициентов группы нагрузок
%на вход передаются:
% P0  — вектор номинальных активных мощностей нагрузок, Вт
% ZIP — матрица коэффициентов: каждая строка [zp, ip, pp, zq, iq, pq]
%возвращается вектор эквивалентных коэффициентов
function zip_eq = equivalent_zip(P0, ZIP)
  %суммарная активная мощность группы нагрузок
  P_sum = sum(P0);
  %взвешенное усреднение коэффициентов (выражение 6.7)
  zp_eq = sum(P0 .* ZIP(:,1)) / P_sum;
  ip_eq = sum(P0 .* ZIP(:,2)) / P_sum;
  pp_eq = sum(P0 .* ZIP(:,3)) / P_sum;
  zq_eq = sum(P0 .* ZIP(:,4)) / P_sum;
  iq_eq = sum(P0 .* ZIP(:,5)) / P_sum;
  pq_eq = sum(P0 .* ZIP(:,6)) / P_sum;
  zip_eq = [zp_eq, ip_eq, pp_eq, zq_eq, iq_eq, pq_eq];
end
