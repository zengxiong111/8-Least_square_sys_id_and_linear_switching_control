function [system_index, K_o] = alg2(A_all,B_all,C_all,Q,R,sigma_u_2,sigma_w_2,sigma_v_2,h,delta_p)

N = size(A_all,3);
n = size(A_all,1);
m = size(C_all,1);
p = size(B_all,2);

%all possible closed-loop systems
[K_all,G_cl_all,A_t_all,B_t_all,C_t_all] = K_L_G_computation(A_all,B_all,... 
    C_all,Q,R,sigma_w_2,sigma_v_2,h); %delta_p is the probability

%compute least distance and critical direction
[gamma,U_all,V_all] = compute_gamma_critical_direction(G_cl_all);
tau_f = zeros(N,1);
for i =1:N
    %tau_f(i) = int16(4*(sigma_w_2+sigma_v_2)/(sigma_u_2*gamma(i))*log(N^2/delta_p));
    %tau_f(i) = int16(4*(sigma_w_2+sigma_v_2)/(sigma_u_2*gamma(i)^2)*log(N^2/delta_p));
    tau_f(i) = int16(max([800/9 * log(3*N^2/delta_p), 4*800*(sigma_w_2+sigma_v_2)/(9*gamma(i)^2*sigma_u_2),...
        4*3200*(sigma_w_2+sigma_v_2)/(9 * sigma_u_2 * gamma(i)^2)* log(16000*N^4/(delta_p^2)) ] ));
end

%What are the effects of gamma(i) and N here?
 
M_all = zeros(N,1);
tau_all = zeros(N,1);
Xi_all = zeros(N,1);

[epsilon_a,epsilon_c] = compute_epsilon(A_t_all,C_t_all);

[M_all,tau_all,Xi_all] = inputs_alg2(A_t_all,B_t_all,C_t_all,epsilon_a,epsilon_c,...
    sigma_w_2,sigma_u_2,sigma_v_2,delta_p);

T_all = sum(tau_all);

A = A_all(:,:,N);
B = B_all(:,:,N);
C = C_all(:,:,N);

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

ynow = C*xnow + Z(:,time_index);

for i =1:N
    Know = K_all(:,:,i); 
    y_energy = 0;
    for t=1:tau_all(i)

        u = Know * ynow;  
        xnext = A*xnow+B*u + W(:,time_index);
        ynow = C*xnow + Z(:,time_index);

        x(:,time_index+1) = xnext;
        xnow = xnext;
        
        time_index = time_index+1;
        y_energy = y_energy + norm(ynow); 
    end
    if (y_energy < Xi_all(i) | i==N)
        %noise trajectory
        W = mvnrnd(zeros(n,1),sigma_w_2*eye(n),tau_f(i));
        W = W';
        Z = mvnrnd(zeros(m,1),sigma_v_2*eye(m),tau_f(i));
        Z = Z';
        Ue = mvnrnd(zeros(p,1),sigma_u_2*eye(p),tau_f(i));
        Ue = Ue';
        Y_id_all = ynow;
        U_id_all = u;
        for t = 1:tau_f(i)
            u = Know * ynow + Ue(:,t);  
            xnext = A*xnow + B*u + W(:,t);
            ynow = C*xnow + Z(:,t);
             
            x(:,time_index+1) = xnext;
            xnow = xnext;
            Y_id_all = [Y_id_all ynow];
            U_id_all = [U_id_all Ue(:,t)];
        end
        G_i_all = G_cl_all(:,:,:,i);
        U_i_all = U_all(:,:,:,i);
        V_i_all = V_all(:,:,:,i);
        system_index = alg1(G_i_all,U_i_all,V_i_all,Y_id_all,U_id_all,h);
        K_o = K_all(:,:,i);
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
