r = 1.1; %r is the spetral radius of A
m = 1; %m is the output dimension
n = 3; %n is the system state dimension
p = 1; %p is the control input dimension
N = 5; % the number of system candidates
h = 5; % The length of the time horizon of Markov parameter matrix 
delta=0.1; %The probability of failure

A_all = zeros(N,n,n);
B_all = zeros(N,n,p );
C_all = zeros(N,m,n );
K_all = zoers(N,p,n);
L_all = zoers(N,n,m);
G_cl_all = zeros(N,N,m,h*p);
A_t_all = zeros(N,N,2*n,2*n);
B_t_all = zeros(N,N,2*n,p );
C_t_all = zeros(N,N,m,2*n );

Q = eye(n);
R = eye(p);

sigma_u_2 = 1;
sigma_w_2 = 0.02;
sigma_z_2 = 0.02;

%similar systems generation
[A_all,B_all,C_all] = similar_system_generation(r,m,n,p,N);

%all possible closed-loop systems
[K_all,L_all,G_cl_all,A_t_all,B_t_all,C_t_all] = K_L_G_computation(A_all,B_all,... %delta is the probability
    C_all,Q,R,sigma_w_2,sigma_z_2,h);
 
M_all = zeros(N,1);
tau_all = zeros(N,1);
Xi_all = zeros(N,1);

epsilon_a = ?;
epsilon_c = ?;

[M_all,tau_all,Xi_all] = inputs_alg2(A_t_all,B_t_all,C_t_all,epsilon_a,epsilon_c);



%find the critical direction
%Or just compare the spectral norm distance