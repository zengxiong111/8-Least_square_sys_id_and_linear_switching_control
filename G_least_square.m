% function G_ls = G_least_square(U_single,Y_single,N,m,T,p)
function G_ls = G_least_square(U_single,Y_single,h)

N_hat = size(U_single,1);

N = N_hat - h;


Y = zeros(N,m);
U = zeros(N,T*p);

for i = 1:N
    Y(i,:) = Y_single(:,i+T-1);
    u_bar = U_single(:,i+T-1);
    for j = 1:T-1
        u_bar = [u_bar;U_single(:,i+T-1 - j)];
    end
    U(i,:) = u_bar';
end

    G_ls = (pinv(U)*Y)';

end
