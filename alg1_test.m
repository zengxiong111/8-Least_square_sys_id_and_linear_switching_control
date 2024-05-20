clear
r = 0.2; %r is the spetral radius of A
n = 2; %n is the system state dimension
m = 1; %m is the output dimension
p = 1; %p is the control input dimension
N = 3; % the number of system candidates
h = 4; % The length of the time horizon of Markov parameter matrix 
delta_p=0.1; %The probability of failure

Length_traj = 100;

sigma_u_2 = 1;
sigma_w_2 = 0.01;
sigma_v_2 = 0.01;

%similar systems generation
[A_all,B_all,C_all] = similar_system_generation(r,m,n,p,N);

A = A_all(:,:,N);
B = B_all(:,:,N);
C = C_all(:,:,N);
D = zeros(m,p);

G_all = zeros(m,h*p,N);
 
for i = 1:N
    G_all(:,1:p,i) = zeros(m,p);
    for k = 1:h-1
        G_all(:,k*p + 1: (k+1)*p,i) = C_all(:,:,i) * A_all(:,:,i)^(k-1)* B_all(:,:,i);
    end
     
end

[U_single,Y_single] = single_trajectory_generation(Length_traj,h,A,B,C,D,sigma_u_2,sigma_w_2,sigma_v_2);

G_hat = G_least_square(U_single,Y_single,h);

U_all = zeros(m,N,N);
V_all = zeros(h*p,N,N);

%critical direction is actually the first lest-singular vector and
%right-singular vector
for i = 1 : N
    for j = 1 : N 
        [U,S,V] = svd(G_all(:,:,i) - G_all(:,:,j));
        U_all(:,i,j) = U(:,1);
        V_all(:,i,j) = V(:,1);
    end
end
 
i=1;
for j=2:N
    if(abs(U_all(:,i,j)'*(G_all(:,:,i) - G_hat)*V_all(:,i,j)) <=...
            abs(U_all(:,i,j)'*(G_all(:,:,j) - G_hat)*V_all(:,i,j)))
        i=j;
    end
end

Ali_system_index = i;

%Find the minimum operator norm!
%And change the proof!

G_G_hat_norm = zeros(N,1);
Operator_system_index = 1;
min_temp = norm(G_all(:,:,1) - G_hat);

for i =2:N
    if(norm(G_all(:,:,i) - G_hat) < min_temp)
        min_temp = norm(G_all(:,:,i) - G_hat);
        Operator_system_index = i;
    end
end


% Estimation error
fprintf('    The index of the ture system:  %d\n ', N);
fprintf('    The index from the  method of Ali:  %d\n', Ali_system_index);
fprintf('    The index from the minimum operator norm:  %d\n', Operator_system_index);




 