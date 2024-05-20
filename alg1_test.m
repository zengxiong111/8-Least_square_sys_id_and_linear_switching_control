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
% G_cl_all = zeros(m,h*p,N,N);

Q = eye(n);
R = eye(p);

sigma_u_2 = 1;
sigma_w_2 = 0.01;
sigma_v_2 = 0.01;

%similar systems generation
[A_all,B_all,C_all] = similar_system_generation(r,m,n,p,N);



A = A_all(:,:,N);
B = B_all(:,:,N);
C = C_all(:,:,N);

G_cl_all = zeros(m,h*p,N);
 
for i = 1:N
    A_h =  A_all(:,:,i) ;
    B_h =  B_all(:,:,i) ;
    C_h =  C_all(:,:,i) ;

    G_cl_all(:,1:p,i,j) = zeros(m,p);
    for k = 1:h-1
        G_cl_all(:,k*p + 1: (k+1)*p,i,j) = C_t_all(:,:,i,j) * A_t_all(:,:,i,j)^(k-1)* B_t_all(:,:,i,j);
    end
     
end

[U_single,Y_single] = single_trajectory_generation(N,T,A,B,C,D,sigma_u_2,sigma_w_2,sigma_v_2)
G_ls = G_least_square(U_single,Y_single,h)

N=size(G_i_all,3);

G_hat = G_least_square(U_id_all,Y_id_all,h);

i=1;
for j=2:N
    if(abs(U_i_all(:,i,j)'*(G_i_all(:,:,i) - G_hat)*V_i_all(:,i,j)) <=...
            abs(U_i_all(:,i,j)'*(G_i_all(:,:,j) - G_hat)*V_i_all(:,i,j)))
        i=j;
    end
end

system_index = i;


 