% Учёт трёхобмоточного трансформатора (автотрансформатора) в матрице Y_bus
% Выражения (7.25)–(7.40)
% на вход передаются:
% Ybus — текущая матрица проводимостей
% br   — структура с полями: from (ВН), mid (СН), to (НН),
%        R1, X1 (ВН), R2, X2 (СН), R3, X3 (НН),
%        G, B (поперечная проводимость на стороне ВН),
%        U1, U2, U3 (номинальные напряжения обмоток)
% на выход: обновлённая матрица Ybus
function Ybus = ybus_add_t3(Ybus, br)
  i = br.from;  % узел ВН
  j = br.mid;   % узел СН
  k = br.to;    % узел НН

  % сопротивления ветвей
  Z_v = br.R1 + 1i * br.X1;
  Z_s = br.R2 + 1i * br.X2;
  Z_n = br.R3 + 1i * br.X3;

  % поперечная проводимость холостого хода
  Y_sh = br.G + 1i * br.B;

  % коэффициенты трансформации
  k_vs = br.U1 / br.U2;
  k_vn = br.U1 / br.U3;

  % вспомогательная величина (7.25)
  Z_ijk = Z_v * Z_s + Z_v * Z_n + Z_s * Z_n;

  % собственные проводимости ветвей (7.26)–(7.28)
  Y_v_ii = 1 / Z_v - Z_s * Z_n / (Z_ijk * Z_v);
  Y_s_jj = k_vs^2 / Z_s - Z_v * Z_n * k_vs^2 / (Z_ijk * Z_s);
  Y_n_kk = k_vn^2 / Z_n - Z_v * Z_s * k_vn^2 / (Z_ijk * Z_n);

  % взаимные проводимости (7.29)–(7.31)
  Y_vs_ij = -Z_n * k_vs / Z_ijk;
  Y_vn_ik = -Z_s * k_vn / Z_ijk;
  Y_sn_jk = -Z_v * k_vs * k_vn / Z_ijk;

  % обновление диагональных элементов (7.32)–(7.34)
  Ybus(i, i) = Ybus(i, i) + Y_v_ii + Y_sh;
  Ybus(j, j) = Ybus(j, j) + Y_s_jj;
  Ybus(k, k) = Ybus(k, k) + Y_n_kk;

  % обновление внедиагональных элементов (7.35)–(7.40)
  Ybus(i, j) = Ybus(i, j) + Y_vs_ij;
  Ybus(j, i) = Ybus(j, i) + Y_vs_ij;
  Ybus(i, k) = Ybus(i, k) + Y_vn_ik;
  Ybus(k, i) = Ybus(k, i) + Y_vn_ik;
  Ybus(j, k) = Ybus(j, k) + Y_sn_jk;
  Ybus(k, j) = Ybus(k, j) + Y_sn_jk;
end
