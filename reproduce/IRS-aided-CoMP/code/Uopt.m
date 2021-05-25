function [U] = Uopt (p, H, W1, W2, PHI)

H1b = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H2b = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;
Hb = [H1b H2b];

W = [W1 ; W2];

U = (Hb * (W * W') * Hb' + p.np.^2 * eye(size(Hb * (W * W') * Hb'))) \ (Hb * W);

end