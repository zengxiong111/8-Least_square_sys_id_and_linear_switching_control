function [K_LQG_all,L_LQG_all] = K_L_computation(A_all,B_all,C_all,Q,R,sigma_w,sigma_z)

 
[K,S1,e1] = dlqr(A,B,Q,R);
[K_temp,S2,e2] = dlqr(A',C',sigma_w * eye(n),sigma_z * eye(p)) ;
L = K_temp';

eig(A-B*K);
eig(A-L*C);

%A_h = [A-B*K B*K;zeros(n,n) A-L*C];
A_h = [A -B*K; L*C A-L*C-B*K];
%B_h = [B;0];
B_h = [B;B];
C_h = [C zeros(m,n)];
Ob_h = obsv(A_h,C_h);
Co_h = ctrb(A_h,B_h);
rank(Ob_h)
rank(Co_h)


 

end