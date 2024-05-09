function system_index = alg1(G_i_all,U_i_all,V_i_all,Y_id_all,U_id_all)

N=size(G_i_all,1);

G_hat = G_least_square(Y_id_all,U_id_all);

i=1;
for j=2:N
    if(absolute(U_i_all(:,i,j)'*(G_i_all(i,:,:) - G_hat)*V_i_all(:,i,j)) <=...
            absolute(U_i_all(:,i,j)'*(G_i_all(j,:,:) - G_hat)*V_i_all(:,i,j)))
        i=j;
    end
end

system_index = i;

%How about the operator norm computation?

end

