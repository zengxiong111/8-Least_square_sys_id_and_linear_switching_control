function [K_all,G_cl_all,A_t_all,B_t_all,C_t_all] = K_L_G_computation(A_all,B_all,C_all,Q,R,sigma_w_2,sigma_v_2,h)

%delta is the probability.

N = size(A_all,3);
n = size(A_all,1);
p = size(B_all,2);
m = size(C_all,1);

K_all = zeros(p,n,N);

G_cl_all = zeros(m,h*p,N,N);

for i = 1:N
    A =  A_all(:,:,i) ;
    B =  B_all(:,:,i) ;
    C =  C_all(:,:,i) ;
    [K_temp,S1,e1] = dlqr(A, B, Q,R);
    K_all(:,:,i) = K_temp*inv(C);
end

A_t_all = zeros( n, n,N,N);
B_t_all = zeros( n,p,N,N );
C_t_all = zeros(m, n,N,N );

 
for i = 1:N
    A_h =  A_all(:,:,i) ;
    B_h =  B_all(:,:,i) ;
    C_h =  C_all(:,:,i) ;
    K_h =  K_all(:,:,i) ;
    for j=1:N
        A =  A_all(:,:,j) ;
        B =  B_all(:,:,j) ;
        C =  C_all(:,:,j) ;
        A_t_all(:,:,i,j) = A-B*K_h*C;
        B_t_all(:,:,i,j) = B ;
        C_t_all(:,:,i,j) =  C ;

        G_cl_all(:,1:p,i,j) = zeros(m,p);
        for k = 1:h-1
            G_cl_all(:,k*p + 1: (k+1)*p,i,j) = C_t_all(:,:,i,j) * A_t_all(:,:,i,j)^(k-1)* B_t_all(:,:,i,j);
        end
    end
end

end