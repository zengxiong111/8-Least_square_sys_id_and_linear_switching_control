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

m_a=0;
m_s=0;
m_p=0;
m_t=0;
c_s=0;


%For m_p,m_t,c_s, we need to check the stability of the closed-loop.
for i=1:N
    for j=1:N
        m_a_temp = max([1,norm(A_t_all(i,j,:,:))]);
        m_s_temp = max([1,norm(B_t_all(i,j,:,:)), norm(C_t_all(i,j,:,:))]);
        if(m_a_temp > m_a)
            m_a = m_a_temp;
        end
        if(m_s_temp > m_s)
            m_s = m_s_temp;
        end
    end
end

for i=1:N
    for j=1:N
        if(vrho(A_t_all(i,j,:,:))<1)
            P = dlyap(A_t_all(i,j,:,:),C_t_all(i,j,:,:)'*C_t_all(i,j,:,:));
            m_p_temp = norm(P)
            m_t_temp = trace(P)
            c_s_temp = 5*(sigma_w_2*trace(P)+sigma_u_2*...
                trace(B_t_all(i,j,:,:)'*P*B_t_all(i,j,:,:))+...
                sigma_z_2*n)*log(1/delta);
        end
    end
end


sigma_m = max([sigma_w,sigma_u,sigma_z]);
c_e = max([1,1/log(1+epsilon_a)]);
c_r = m_p * (22*epsilon_c^(-2)+1)*sigma_m^2*c_e;
c_p = 2*max([1,m_p/(epsilon_c^2)]);


M_all(1) = 5 * m_p * log(1/delta);


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