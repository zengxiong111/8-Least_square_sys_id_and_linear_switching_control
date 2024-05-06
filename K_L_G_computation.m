function [K_all,L_all,G_cl_all] = K_L_G_computation(A_all,B_all,C_all,Q,R,sigma_w,sigma_z,h)

N = size(A_all,1);
n = size(A_all,2);
p = size(B_all,3);
m = size(C_all,2);

K_all = zoers(N,p,n);
L_all = zoers(N,n,m);

G_cl_all = zeros(N,N,m,h*p);

for i = 1:N
    A = squeeze(A_all(i,:,:));
    B = squeeze(B_all(i,:,:));
    C = squeeze(C_all(i,:,:));
    [K_all(i,:,:),S1,e1] = dlqr(A, B, Q,R);
    [K_temp,S2,e2] = dlqr(A',C',sigma_w * eye(n),sigma_z * eye(p)) ;
    L_all(i,:,:) = K_temp';
end

 
for i = 1:N
    A_h = squeeze(A_all(i,:,:));
    B_h = squeeze(B_all(i,:,:));
    C_h = squeeze(C_all(i,:,:));
    K_h = squeeze(K_all(i,:,:));
    L_h = squeeze(L_all(i,:,:));
    for j=1:N
        A = squeeze(A_all(j,:,:));
        B = squeeze(B_all(j,:,:));
        C = squeeze(C_all(j,:,:));
        A_t = [A-B*K_h   B*K_h; ...
            (A_h-A)+(B_h-B)*K_h+L_h*(C-C_h)   A_h-L_h*C_h+(B_h-B)*K_h];
        B_t = [B;zeros(size(B))];
        C_t = [C zeros(size(C))];

        G_cl_all(i,j,:,1:p) = D;
        for k = 1:T-1
            G_cl_all(i,j,:,k*p + 1: (k+1)*p) = C_t*A_t^(k-1)*B_t;
        end

        %find the critical direction
        
    end
end




 

end