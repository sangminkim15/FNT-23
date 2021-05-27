function [Q] = Qopt (p, H, W1, W2, PHI)

H1b = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H2b = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;
Hb = [H1b H2b];

W = [W1 ; W2];

J = (Hb * (W * W') * Hb') + p.np * eye(p.N_t);

Q = inv(eye(p.N_r) - ((W' * Hb') / J) * Hb * W);

end