function system_index = alg2(A,B,C,G_cl_all,A_t_all,B_t_all,C_t_all,M_all,tau_all,Xi_all,K_all,L_all)

n = size(A,1);
N = size(G_cl_all,1);
x0 = randn(n,1);
xnow = x0;

T_all = sum(tau_all);

x = zeros(n,T_all+1);
x(:,1) = xnow;

%valid_clf = 1:14;
indnow = 1;
Know = K_all(indnow,:); 
Lnow = L_all(indnow,:);

for i =1:N
    for t=1:T
        u = Know * xnow; %solve for u
        xnext = A*xnow+B1*u;
        x(:,t+1) = xnext;
        % check CLF validity
        if xnext'*Pnow*xnext >= xnow'*Pnow*xnow
            t
            indnow = indnow+1;
            Know = K_all(indnow,:); 
            Pnow = P_all( ((indnow-1)*n+1):(n*indnow),:);
            %norm(xnow)
        end
        xnow = xnext;
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
