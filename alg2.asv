function system_index = alg2(A,B,C,G_cl_all,A_t_all,B_t_all,C_t_all,M_all,tau_all,Xi_all,K_all,L_all)

n = size(A,1);
m = size(C,1);
N = size(G_cl_all,1);
x0 = randn(n,1);
xnow = x0;
x_hat_now = x0;

T_all = sum(tau_all);

%valid_clf = 1:14;
indnow = 1;

%noise trajectory
W = mvnrnd(zeros(n,1),sigma_w*eye(n),T_all);
W = W';
Z = mvnrnd(zeros(m,1),sigma_z*eye(m),T_all);
Z = Z';



time_index = 1;
x = zeros(n,T_all+1);
x(:,time_index) = xnow;

for i =1:N
    Know = K_all(i,:,:); 
    Lnow = L_all(i,:,:);
    for t=1:tau_all(i)
        u = Know * x_hat_now;  

        xnext = A*xnow+B*u + W(:,time_index);
        ynow = C*xnow + Z(:,time_index);

        x(:,time_index+1) = xnext;
        xnow = xnext;
        time_index = time_index+1;
    end
end

x_l2_norm = zeros(T+1,1);
for i = 1:T+1
    x_l2_norm(i) = norm(x(1:2,i));
end

 

figure;
hold on;
plot(0:20,x_l2_norm(1:20+1),'-x','LineWidth',3); % plot the tight bound eq (4) in the paper

 %legend('robustness margin from Hinf','LMI sufficient [0,+m]','LMI sufficient [-m,+m]','controllability margin')
legend('Adaptively Switching CLF')
grid on;
ax = gca;
ax.LineWidth = 2;
ax.GridLineStyle = '--';
ax.GridAlpha = 0.8;
lgd.FontSize = 18;
xlabel('time t','FontSize',18) ;
ylabel('||x_t||_2','FontSize',18) ;
set(gca,'FontSize',20)
