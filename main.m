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
[K_all,L_all,G_cl_all,A_t_all,B_t_all,C_t_all] = K_L_G_computation(A_all,B_all,... %delta is the probability
    C_all,Q,R,sigma_w,sigma_z,h);
 

P_all = zeros(N,N,2*n,2*n);

%find all Xi in Proposition 5
Xi_max = find_Xi(A_t_all,B_t_all,C_t_all,delta,sigma_w,sigma_z,sigma_u);
delta=0.1;
P_all(i,j,:,:) = dlyap(A_t_all(i,j,:,:),C_t_all(i,j,:,:)'*C_t_all(i,j,:,:));

P_op_all = zeros(N,N);

%x_1 is sampled from N(0,I_{d_x*t \times d_x})

m_a
m_s
m_p = max(P_op_all);
m_t
sigma_m
c_e
c_r
c_p = 2*max(1,m_p / epsilon_c / epsilon_c);
c_s

M_all = zeros(N,1);
M_all(1) = 5 * m_p * log(1/delta);

tau_all = zeros(N,1);
%tau_all(1) = n + log(1/delta);
tau_all(1) = max(1600/9*log(1/delta),...
    6400*epsilon_c/(9*sigma_w^2*(M_1 + c_r *log(2/delta)*n^2 * m_a^(4*n) + c_s * log(1/delta) )));
for j= 2:N
    tau_all(j) = (j-1)*(2/log(m_a)*n + log(c_p)) + tau_all(1);
end

for j=2:N
    M_all(j) = m_p/(epsilon_c^2)*m_a^(2*n)*(M_all(j-1)+tau_all(j-1)*c_s)...
        + c_r * log(1/delta)* n^2 * m_a^(4*n);
end


Xi_all(i,j) ;

%find the critical direction
%Or just compare the spectral norm distance