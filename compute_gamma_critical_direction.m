function [gamma,U_all,V_all] = compute_gamma_critical_direction(G_cl_all)

N = size(G_cl_all,3);
n = size(G_cl_all,1);
m = size(G_cl_all,2);
gamma = 1e16*ones(N,1);
U_all = zeros(n,N,N,N);
V_all = zeros(m,N,N,N);

%critical direction is actually the first lest-singular vector and
%right-singular vector
 
 
for i = 1 : N
    num=1;
    for j = 1 : N-1
        for k = j+1 : N
            gamma_temp = norm(G_cl_all(:,:,j,i)-G_cl_all(:,:,k,i));
            if(gamma_temp<gamma(i))
                gamma(i) = gamma_temp;
            end
        end
    end
end

for i = 1 : N
    for j = 1 : N 
        for k = 1 : N
            [U,S,V] = svd(G_cl_all(:,:,j,i)-G_cl_all(:,:,k,i));
            U_all(:,j,k,i) = U(:,1);
            V_all(:,j,k,i) = V(:,1);
        end
    end
end
 