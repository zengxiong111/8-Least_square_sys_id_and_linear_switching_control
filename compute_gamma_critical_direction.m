function [gamma,U_all,V_all] = compute_gamma(G_cl_all)

N = size(G_cl_all,1);
n = size(G_cl_all,3);
m = size(G_cl_all,4);
gamma = 1e16*ones(N,1);
U_all = zeros(n,N,N*(N-1)/2);
V_all = zeros(m,N,N*(N-1)/2);

%critical direction is actually the first lest-singular vector and
%right-singular vector
 
 
for i = 1 : N
    num=1;
    for j = 1 : N-1
        for k = j+1 : N
            gamma_temp = norm(G_cl_all(i,j,:,:)-G_cl_all(i,k,:,:));
            if(gamma_temp<gamma(i))
                gamma(i) = gamma_temp;
            end
            [U,S,V] = svd(G_cl_all(i,j,:,:)-G_cl_all(i,k,:,:));
            U_all(:,i,num) = U(:,1);
            V_all(:,i,num) = V(:,1);
            num = num+1;
        end
    end
end
 