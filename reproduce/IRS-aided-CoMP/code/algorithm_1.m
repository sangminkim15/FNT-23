function [W] = algorithm_1 (p, H, U, Q, PHI)

H1b = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H2b = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;

muprev = [100 100];

J1 = [(H1b * U * Q * U' * H1b' + muprev(1) * eye(size(H1b * H1b'))) (H1b * U * Q * U' * H2b') ; (H2b * U * Q * U' * H1b') (H2b * U * Q * U' * H2b' + muprev(2) * eye(size(H2b * H2b')))];
J2 = [H1b' * U * Q ; H2b' * U * Q];
W = J1 \ J2;

sz = size(W);
W1 = W(1:sz(1)/2, 1:sz(2));
W2 = W(sz(1)/2+1:sz(1), 1:sz(2));

L = Lagrangian (p, H, U, Q, W1, W2, PHI, muprev);

while true
    mu(1) = muprev(1) + muprev(1)/ 10000 * (norm(W1, 'fro').^2 - p.P_max);
    mu(2) = muprev(2) + muprev(2)/ 10000 * (norm(W2, 'fro').^2 - p.P_max);
    
    if mu(1) < 0
        mu(1) = 0;
    end
    
    if mu(2) < 0
        mu(2) = 0;
    end
    
    muprev = mu;
    
    J1 = [(H1b * U * Q * U' * H1b' + muprev(1) * eye(size(H1b * H1b'))) (H1b * U * Q * U' * H2b') ; (H2b * U * Q * U' * H1b') (H2b * U * Q * U' * H2b' + muprev(2) * eye(size(H2b * H2b')))];
    J2 = [H1b' * U * Q ; H2b' * U * Q];
    W = J1 \ J2;

    W1 = W(1:sz(1)/2, 1:sz(2));
    W2 = W(sz(1)/2+1:sz(1), 1:sz(2));
    
    Lprev = L;
    L = Lagrangian(p, H, U, Q, W1, W2, PHI, muprev);
    
    if abs(L - Lprev) < p.epsilon
        break;
    end
end

end