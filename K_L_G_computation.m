function [K_all,L_all,G_cl_all] = K_L_G_computation(A_all,B_all,C_all,Q,R,sigma_w,sigma_z,h)

N = size(A_all,1);
n = size(A_all,2);
p = size(B_all,3);
m = size(C_all,2);

K_all = zoers(N,p,n);
L_all = zoers(N,n,m);

for i = 1:N
    [K_all(i,:,:),S1,e1] = dlqr(A_all(i,:,:),B_all(i,:,:),Q,R);
    [K_temp,S2,e2] = dlqr(A_all(i,:,:)',C_all(i,:,:)',sigma_w * eye(n),sigma_z * eye(p)) ;
    L_all(i,:,:) = K_temp';
end

 
for i = 1:N
    for j=1:N
        G = zeros(m,h*p);
        %A_h = [A-B*K B*K;zeros(n,n) A-L*C];
        A_h = [A -B*K; L*C A-L*C-B*K];
        %B_h = [B;0];
        B_h = [B;B];
        C_h = [C zeros(m,n)];

        G(:,1:p) = D;
        for k = 1:T-1
            G(:,k*p + 1: (k+1)*p) = C*A^(k-1)*B;
        end
    end
end




 

end