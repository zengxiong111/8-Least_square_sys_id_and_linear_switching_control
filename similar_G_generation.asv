function G_closed_loop_all = similar_G_generation(A_all,B_all,C_all,K_LQG_all,L_LQG_all)


G = zeros(m,T*p);
G(:,1:p) = D;
for i = 1:T-1
    G(:,i*p + 1: (i+1)*p) = C*A^(i-1)*B;
end

 

end