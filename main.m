r = 1.1; %r is the spetral radius of A
m = 1; %m is the output dimension
n = 3; %n is the system state dimension
p = 1; %p is the control input dimension
N =5; % the number of system candidates
h=5; % The length of the time horizon of Markov parameter matrix 

A_all = zeros(N,n,n);
B_all = zeros(N,n,p );
C_all = zeros(N,m,n );
G_closed_loop_all = zeros(N,N,m,T*p);
K_LQG_all = zoers(N,p,n);
L_LQG_all = zoers(N,n,m);

Q = eye(n);
R = eye(p);

sigma_u = 1;
sigma_w = 0.02;
sigma_z = 0.02;

[A_all,B_all,C_all] = similar_system_generation(r,m,n,p,N);
[K_LQG_all,L_LQG_all] = K_L_computation(A_all,B_all,C_all,Q,R,sigma_w,sigma_z);
G_closed_loop_all = similar_G_generation(A_all,B_all,C_all,K_LQG_all,L_LQG_all);

%A bug in the paper is that r is not determined by your assumption and it
%is determined by your controllers and your system properties.

[K,S1,e1] = dlqr(A,B,Q,R) ;
[K_temp,S2,e2] = dlqr(A',C',sigma_w * eye(n),sigma_z * eye(p)) ;
L = K_temp';

eig(A-B*K);
eig(A-L*C);

%A_h = [A-B*K B*K;zeros(n,n) A-L*C];
A_h = [A -B*K; L*C A-L*C-B*K];
%B_h = [B;0];
B_h = [B;B];
C_h = [C zeros(m,n)];
Ob_h = obsv(A_h,C_h);
Co_h = ctrb(A_h,B_h);
rank(Ob_h)
rank(Co_h)

[T,Abar,Bbar,Cbar]=Kalman_Decomposition(A_h,B_h,C_h);

% A_hat = A + 0.01*ones(size(A));
% B_hat = B + 0.01*ones(size(B));
% C_hat = C + 0.01*ones(size(C));
% 
% [K_hat,S1_hat,e1_hat] = dlqr(A_hat,B_hat,Q,R) ;
% [K_temp_hat,S2_hat,e2_hat] = dlqr(A_hat',C_hat',sigma_w * eye(n),sigma_z * eye(p)) ;
% L_hat = K_temp_hat';
% 
% eig(A-B*K_hat)
% eig(A-L_hat*C)
% 
% %A_h = [A-B*K B*K;zeros(n,n) A-L*C];
% A_h_hat = [A -B*K_hat; L_hat*C A_hat-L_hat*C_hat-B_hat*K_hat];
% %B_h = [B;0];
% B_h_hat = [B;B_hat];
% C_h_hat = [C zeros(m,n)];
% Ob_h_hat = obsv(A_h_hat,C_h_hat);
% Co_h_hat = ctrb(A_h_hat,B_h_hat);
% rank(Ob_h_hat)
% rank(Co_h_hat)