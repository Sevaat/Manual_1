%на вход необходимо передать следующий массив, содержащий данные трансформатора:
%t2 = [n_tr, dP_id, U_h_nom, U_l_nom, S_nom, u_sc, dP_1, I_id]
function [R, X, Y, k] = get_t2_parameters(t2)
  n_tr = t2(1);
  dP_id = t2(2);
  U_h_nom = t2(3);
  U_l_nom = t2(4);
  S_nom = t2(5);
  u_sc = t2(6);
  dP_1 = t2(7);
  I_id = t2(8);
  %расчет активного сопротивления в Ом
  R = 1 / n_tr * dP_id * (U_h_nom / S_nom)^2;
  %расчет реактивного сопротивления в Ом
  Z = 1 / n_tr * u_sc / 100 * U_h_nom^2 / S_nom;
  X = sqrt(max(0, Z^2 - R^2));
  %расчет комплекса полной мощности в В*А
  S_id = n_tr * (dP_1 + 1i * I_id / 100 * S_nom);
  %расчет комплекса проводимости в См
  Y = conj(S_id) / U_h_nom^2;
  %расчет коэффициента трансформации
  k = U_h_nom / U_l_nom;
end
