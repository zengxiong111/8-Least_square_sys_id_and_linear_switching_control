function system_index = alg2(A_all,B_all,C_all,G_cl_all,tau_all,Xi_all,tau_f,U_all,V_all,K_all,L_all,A,B,C,h,sigma_w_2,sigma_v_2)

N = size(G_cl_all,3);
T_all = sum(tau_all);

%noise trajectory
W = mvnrnd(zeros(n,1),sigma_w_2*eye(n),T_all);
W = W';
Z = mvnrnd(zeros(m,1),sigma_v_2*eye(m),T_all);
Z = Z';

%x0 = randn(n,1);
x0 = zeros(n,1);
xnow = x0;
x_hat_now = x0;

time_index = 1;
x = zeros(n,T_all+1);
x(:,time_index) = xnow;

for i =1:N
    Know = K_all(:,:,i); 
    Lnow = L_all(:,:,i);
    y_energy = 0;
    for t=1:tau_all(i)

        u = Know * x_hat_now;  
        xnext = A*xnow+B*u + W(:,time_index);
        ynow = C*xnow + Z(:,time_index);
        x_hat_next = A_all(:,:,i)*x_hat_now+...
            B_all(:,:,i)*u+Lnow*(ynow - C_all(:,:,i)*x_hat_now);

        x(:,time_index+1) = xnext;
        xnow = xnext;
        x_hat_now = x_hat_next;
        time_index = time_index+1;
        y_energy = y_energy + norm(ynow); 
    end
    if y_energy < Xi_all(i)
        %noise trajectory
        W = mvnrnd(zeros(n,1),sigma_w_2*eye(n),tau_f(i));
        W = W';
        Z = mvnrnd(zeros(m,1),sigma_z_2*eye(m),tau_f(i));
        Z = Z';
        Ue = mvnrnd(zeros(p,1),sigma_u_2*eye(p),tau_f(i));
        Ue = Ue';
        Y_id_all = ynow;
        U_id_all = u;
        for t = 1:tau_f(i)
            u = Know * x_hat_now + Ue(:,t);  
            xnext = A*xnow+B*u + W(:,t);
            ynow = C*xnow + Z(:,t);
            x_hat_next = A_all(:,:,i)*x_hat_now+...
                B_all(:,:,i)*u+Lnow*(ynow - C_all(:,:,i)*x_hat_now);
            x(:,time_index+1) = xnext;
            xnow = xnext;
            x_hat_now = x_hat_next;
            Y_id_all = [Y_id_all ynow];
            U_id_all = [U_id_all Ue(:,t)];
        end
        G_i_all = G_cl_all(:,:,:,i);
        U_i_all = U_all(:,:,:,i);
        V_i_all = V_all(:,:,:,i);
        system_index = alg1(G_i_all,U_i_all,V_i_all,Y_id_all,U_id_all,h);
        return
    end
end

% x_l2_norm = zeros(T_all,1);
% for i = 1:T_all
%     x_l2_norm(i) = norm(x(:,i));
% end
% 
% figure;
% hold on;
% plot(0:20,x_l2_norm(1:20+1),'-x','LineWidth',3); % plot the tight bound eq (4) in the paper
% 
%  %legend('robustness margin from Hinf','LMI sufficient [0,+m]','LMI sufficient [-m,+m]','controllability margin')
% legend('Adaptively Switching CLF')
% grid on;
% ax = gca;
% ax.LineWidth = 2;
% ax.GridLineStyle = '--';
% ax.GridAlpha = 0.8;
% lgd.FontSize = 18;
% xlabel('time t','FontSize',18) ;
% ylabel('||x_t||_2','FontSize',18) ;
% set(gca,'FontSize',20)
