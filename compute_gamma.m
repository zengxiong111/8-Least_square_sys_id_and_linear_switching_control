function gamma = compute_gamma(G_cl_all)

N = size(G_cl_all);
gamma = 1e16*ones(N,1);

 
num=0;
for i = 1 : N
    for j = 1 : N-1
        for k = j+1 : N
            gamma_temp = norm(G_cl_all(i,j,:,:)-G_cl_all(i,k,:,:));
            if(gamma_temp<gamma(i))
                gamma(i) = gamma_temp;
        end
    end
end
 