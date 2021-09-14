function [U] = Uopt_MCEU (p, H, W)

H11 = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H21 = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;      H31 = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;
H12 = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H22 = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;      H32 = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;
H13 = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H23 = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;      H33 = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;

H1 = [H11 H21 H31];                                 H2 = [H12 H22 H32];                                 H3 = [H13 H23 H33];

W1 = [W.W11 ; W.W21 ; W.W31];                       W2 = [W.W12 ; W.W22 ; W.W32];                       W3 = [W.W13 ; W.W23 ; W.W33];

U.U1 = (H1 * (W1 * W1' + W2 * W2' + W3 * W3') * H1' + p.np * eye(p.N_t)) \ (H1 * W1);
U.U2 = (H2 * (W1 * W1' + W2 * W2' + W3 * W3') * H2' + p.np * eye(p.N_t)) \ (H2 * W2);
U.U3 = (H3 * (W1 * W1' + W2 * W2' + W3 * W3') * H3' + p.np * eye(p.N_t)) \ (H3 * W3);

end