clear
r = 0.2; %r is the spetral radius of A
n = 2; %n is the system state dimension
m = 1; %m is the output dimension
p = 1; %p is the control input dimension
N = 3; % the number of system candidates
h = 4; % The length of the time horizon of Markov parameter matrix 
delta_p=0.1; %The probability of failure

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
    G_cl_all(:,1:p,i) = zeros(m,p);
    for k = 1:h-1
        G_cl_all(:,k*p + 1: (k+1)*p,i) = C_all(:,:,i) * A_all(:,:,i)^(k-1)* B_all(:,:,i);
    end
     
end

[U_single,Y_single] = single_trajectory_generation(N,T,A,B,C,D,sigma_u_2,sigma_w_2,sigma_v_2);
G_ls = G_least_square(U_single,Y_single,h);

N=size(G_i_all,3);

G_hat = G_least_square(U_id_all,Y_id_all,h);

U_all = zeros(m,N,N,N);
V_all = zeros(h*p,N,N,N);

%critical direction is actually the first lest-singular vector and
%right-singular vector
for i = 1 : N
    for j = 1 : N 
        for k = 1 : N
            [U,S,V] = svd(G_cl_all(:,:,j,i)-G_cl_all(:,:,k,i));
            U_all(:,j,k,i) = U(:,1);
            V_all(:,j,k,i) = V(:,1);
        end
    end
end
 
i=1;
for j=2:N
    if(abs(U_i_all(:,i,j)'*(G_i_all(:,:,i) - G_hat)*V_i_all(:,i,j)) <=...
            abs(U_i_all(:,i,j)'*(G_i_all(:,:,j) - G_hat)*V_i_all(:,i,j)))
        i=j;
    end
end

Ali_system_index = i;

%Find the minimum operator norm!
%And change the proof!


% Estimation error
fprintf('    The index of the ture system:  %6.3E \n', N);
fprintf('    The index from the  method of Ali:  %6.3E \n', Ali_system_index);
fprintf('    The index from the minimum operator norm:  %6.3E \n', Operator_system_index);




 