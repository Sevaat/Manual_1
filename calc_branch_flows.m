% Расчёт потоков мощностей и потерь в ветвях
% на вход передаются:
%   U        — вектор комплексных напряжений узлов, о.е.
%   branches — массив структур ветвей (поля: type, from, to, R, X, G, B, U1, U2)
% на выход:
%   flows  — структура с полями для каждой ветви:
%            .S_from — мощность в начале ветви
%            .S_to   — мощность в конце ветви
%            .dS     — потери в ветви
%            .I      — ток через продольную ветвь
%   P_loss_total, Q_loss_total — суммарные потери
function [flows, P_loss_total, Q_loss_total] = calc_branch_flows(U, branches)
  n_br = length(branches);
  flows = struct();
  P_loss_total = 0;
  Q_loss_total = 0;

  for l = 1:n_br
    br = branches(l);
    i = br.from;
    j = br.to;

    if strcmp(br.type, 'L')
      % === Линия электропередачи (П-образная схема) ===
      y_series = 1 / (br.R + 1i * br.X);
      y_shunt  = br.G + 1i * br.B;

      % токи (выражения 9.2, 9.4)
      I_from = y_series * (U(i) - U(j)) + 0.5 * y_shunt * U(i);
      I_to   = y_series * (U(j) - U(i)) + 0.5 * y_shunt * U(j);

      % мощности (выражения 9.1, 9.3)
      S_from = U(i) * conj(I_from);
      S_to   = U(j) * conj(I_to);

      % ток через продольную ветвь
      I_series = y_series * (U(i) - U(j));

    elseif strcmp(br.type, 'T2')
      % === Двухобмоточный трансформатор ===
      k = br.U1 / br.U2;
      y_series = 1 / (br.R + 1i * br.X);
      y_shunt  = br.G + 1i * br.B;

      % токи (выражения 9.11, 9.12)
      I_from = y_series * (U(i) - k * U(j)) + y_shunt * U(i);
      I_to   = y_series * k * (k * U(j) - U(i));

      % мощности (выражение 9.13)
      S_from = U(i) * conj(I_from);
      S_to   = U(j) * conj(I_to);

      I_series = y_series * (U(i) - k * U(j));
    end

    % потери в ветви (выражение 9.5)
    dS = S_from + S_to;

    % сохранение результатов
    flows(l).S_from = S_from;
    flows(l).S_to   = S_to;
    flows(l).dS     = dS;
    flows(l).I      = I_series;

    % накопление суммарных потерь
    P_loss_total = P_loss_total + real(dS);
    Q_loss_total = Q_loss_total + imag(dS);
  end
end
