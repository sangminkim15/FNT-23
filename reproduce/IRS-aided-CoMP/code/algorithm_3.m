function [objective_value] = algorithm_3(p,H)


%% 1. Initialization \Wv_n and \mathbf{\Phi}
temp_W_1 = randn(p.N_t,p.d)+j*randn(p.N_t,p.d);
W_1 = temp_W_1/sqrt(trace(temp_W_1'*temp_W_1))*sqrt(p.P_max);

temp_W_2 = randn(p.N_t,p.d)+j*randn(p.N_t,p.d);
W_2 = temp_W_2/sqrt(trace(temp_W_2'*temp_W_2))*sqrt(p.P_max);

W_opt = [W_1; W_2];

%% Initialzization Phi
theta = rand(1,p.M)*2*pi;
Phi_opt = diag(exp(j*theta));

%% 2. Calculte \Uv^{opt} from (15)


error = 1;
objective_value = 1;

% construct \bar{Hv}
H_bar_1 = [H.bs1_ue1+H.IRS_ue1*Phi_opt*H.bs1_IRS];

H_bar_2 = [H.bs2_ue1+H.IRS_ue1*Phi_opt*H.bs2_IRS];

H_bar = [H_bar_1, H_bar_2];
% calculate \Uv^{opt}
U_opt = inv(H_bar*W_opt*W_opt'*H_bar'+p.np*eye(p.N_r))*H_bar*W_opt;

%% 3. Calculte \Qv^{opt} from (16)

Q_opt = inv(U_opt'*(H_bar*W_opt*W_opt'*H_bar'+p.np*eye(p.N_r))*U_opt-U_opt'*H_bar*W_opt-W_opt'*H_bar'*U_opt+eye(p.d));
test = [];
while true
    
    %% 4. Calculate \Wv^{opt}_n from Algorithm 1.
    [W_opt] = algorithm_1(p,H_bar,U_opt,Q_opt,H);
    
    %% 5. Calculate \Phi^{opt}_n from Algorithm 2.
    [Phi_opt] = algorithm_2(p,H,U_opt,Q_opt,W_opt,Phi_opt);
    
    % construct \bar{Hv}
    H_bar_1 = [H.bs1_ue1+H.IRS_ue1*Phi_opt*H.bs1_IRS];
    
    H_bar_2 = [H.bs2_ue1+H.IRS_ue1*Phi_opt*H.bs2_IRS];
    
    H_bar = [H_bar_1, H_bar_2];
    % calculate \Uv^{opt}
    U_opt = inv(H_bar*W_opt*W_opt'*H_bar'+p.np*eye(p.N_r))*H_bar*W_opt;

    Q_opt = inv(U_opt'*(H_bar*W_opt*W_opt'*H_bar'+p.np*eye(p.N_r))*U_opt-U_opt'*H_bar*W_opt-W_opt'*H_bar'*U_opt+eye(p.d));

    objective_temp = real(log(det(Q_opt)));
    test = [test, objective_temp];
    error = abs(objective_value-objective_temp);
    objective_value = objective_temp;
    if error < p.eta
        break;
    end
%     disp(['Algorithm 3 error : ',num2str(error)]);
end
figure
plot(test)
end
