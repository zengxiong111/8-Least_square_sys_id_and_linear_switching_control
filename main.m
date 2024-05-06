r = 1.1; %r is the spetral radius of A
m = 1; %m is the output dimension
n = 3; %n is the system state dimension
p = 1; %p is the control input dimension
N = 5; % the number of system candidates
h = 5; % The length of the time horizon of Markov parameter matrix 

A_all = zeros(N,n,n);
B_all = zeros(N,n,p );
C_all = zeros(N,m,n );
K_all = zoers(N,p,n);
L_all = zoers(N,n,m);
G_cl_all = zeros(N,N,m,h*p);

Q = eye(n);
R = eye(p);

sigma_u = 1;
sigma_w = 0.02;
sigma_z = 0.02;

[A_all,B_all,C_all] = similar_system_generation(r,m,n,p,N);
[K_all,L_all,G_cl_all] = K_L_G_computation(A_all,B_all,C_all,Q,R,sigma_w,sigma_z,h);
 
%A bug in the paper is that r is not determined by your assumption and it
%is determined by your controllers and your system properties.


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