clear
r = 1.01; %r is the spetral radius of A
n = 2; %n is the system state dimension
m = n; %m is the output dimension
p = 1; %p is the control input dimension
N = 3; % the number of system candidates
h = 4; % The length of the time horizon of Markov parameter matrix 
delta_p=0.1; %The probability of failure

% A_all = zeros(n,n,N);
% B_all = zeros(n,p,N );
% C_all = zeros(m,n,N );
% K_all = zeros(p,n,N);
% L_all = zeros(n,m,N);
% G_cl_all = zeros(m,h*p,N,N);
% A_t_all = zeros(2*n,2*n,N,N);
% B_t_all = zeros(2*n,p,N,N );
% C_t_all = zeros(m,2*n,N,N );

Q = eye(n);
R = eye(p);

sigma_u_2 = 1;
sigma_w_2 = 0.01;
sigma_v_2 = 0.01;

%similar systems generation
[A_all,B_all,C_all] = similar_system_generation(r,m,n,p,N);

% %all possible closed-loop systems
% [K_all,L_all,G_cl_all,A_t_all,B_t_all,C_t_all] = K_L_G_computation(A_all,B_all,... 
%     C_all,Q,R,sigma_w_2,sigma_v_2,h); %delta_p is the probability
% 
% %compute least distance and critical direction
% [gamma,U_all,V_all] = compute_gamma_critical_direction(G_cl_all);
% tau_f = zeros(N,1);
% for i =1:N
%     tau_f(i) = 4*(sigma_w_2+sigma_v_2)/(sigma_u_2*gamma(i))*log(N^2/delta_p);
% end
%  
% M_all = zeros(N,1);
% tau_all = zeros(N,1);
% Xi_all = zeros(N,1);
% 
% [epsilon_a,epsilon_c] = compute_epsilon(A_t_all,C_t_all);
% 
% [M_all,tau_all,Xi_all] = inputs_alg2(A_t_all,B_t_all,C_t_all,epsilon_a,epsilon_c,...
%     sigma_w_2,sigma_u_2,sigma_v_2,delta_p);

A = A_all(:,:,N);
B = B_all(:,:,N);
C = C_all(:,:,N);

[system_index, K_o] = alg2(A_all,B_all,C_all,Q,R,sigma_u_2,sigma_w_2,sigma_v_2,h,delta_p);
vrho(A-B*K_o*C)


%This is a random variable, we need to compute the estimation of its
%expectation!

%find the critical direction
%Or just compare the spectral norm distance