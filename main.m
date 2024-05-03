%r is the spetral radius of A
%m is the output dimension
%n is the system state dimension
%p is the control input dimension
r = 1.1;
m = 1;
n = 3;
p = 1;

for i = 1:5
    [A,B,C,D] = system_generation(r,m,n,p);
    Ob = obsv(A,C);
    Co = ctrb(A,B);
    % check if the genrated system is controllable and observable
    if(rank(Ob) == n & rank(Co) == n)
        break;
    end
    if(rank(Ob) < n | rank(Co) < n)
        continue;
    end

end

Q = eye(n);
R = eye(p);

sigma_u = 1;
sigma_w = 0.02;
sigma_z = 0.02;

[K,S1,e1] = dlqr(A,B,Q,R) ;
[K_temp,S2,e2] = dlqr(A',C',sigma_w * eye(n),sigma_z * eye(p)) ;
L = K_temp';

eig(A-B*K)
eig(A-L*C)

%A_h = [A-B*K B*K;zeros(n,n) A-L*C];
A_h = [A -B*K; L*C A-L*C-B*K];
%B_h = [B;0];
B_h = [B;B];
C_h = [C zeros(m,n)];
Ob_h = obsv(A_h,C_h);
    Co_h = ctrb(A_h,B_h);
    rank(Ob_h)
    rank(Co_h)