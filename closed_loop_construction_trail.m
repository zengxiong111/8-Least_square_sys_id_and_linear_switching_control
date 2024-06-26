clear
%r is the spetral radius of A
%m is the output dimension
%n is the system state dimension
%p is the control input dimension
r = 1.1;
n = 2;
m = n;
p = 1;

for i = 1:5
    [A,B,C,D] = system_generation(r,m,n,p);
    Ob = obsv(A,C);
    Co = ctrb(A,B);
    % check if the genrated system is controllable and observable
%     if(rank(Ob) == n & rank(Co) == n & rank(C) == n)
if(rank(Ob) == n & rank(Co) == n )
        break;
    end
%     if(rank(Ob) < n | rank(Co) < n | rank(C) < n)
  if(rank(Ob) < n | rank(Co) < n )
        continue;
    end
end

Q = eye(n);
R = eye(p);

sigma_u = 1;
sigma_w = 0.02;
sigma_z = 0.02;

[K,S1,e1] = dlqr(A,B,Q,R);
Ko = K*inv(C);
vrho(A-B*Ko*C)
Ob_h = obsv(A-B*Ko*C,C);
rank(Ob_h)
svd_Ob_h = svd(Ob_h)
svd(obsv(A,C))

% [K_temp,S2,e2] = dlqr(A',C',sigma_w * eye(n),sigma_z * eye(p));
% L = K_temp';
% 
% eig(A-B*K);
% eig(A-L*C);
% 
% A_h = [A-B*K B*K;zeros(n,n) A-L*C];
% B_h = [B;zeros(size(B))];
% % A_h = [A -B*K; L*C A-L*C-B*K];
% % B_h = [B;B];
% C_h = [C zeros(m,n)];
% Ob_h = obsv(A_h,C_h);
% Co_h = ctrb(A_h,B_h);
% rank(Ob_h)
% svd_Ob_h = svd(Ob_h)
% svd(obsv(A,C))
% rank(Co_h)



function [A,B,C,D] = system_generation(r,m,n,p)

A = zeros(n,n);
B = zeros(n,p);
C = zeros(m,n);
D = zeros(m,p);

B = normrnd(0,1/n,[n,p]);
%C = normrnd(0,1/m,[m,n]);
C = normrnd(0,10,[m,n]);
A = rand(n);
A = r * A/vrho(A);

% A = [r 0.5;0 r];
% B = [0;1];
% C = [1 0];
% D = [0 0];
% 
% A = [1 0.01 0;0.01 1 0.01;0 0.01 1];
% A = r*A/vrho(A);
% B=eye(3);
% C = [1 0 0];
% D = [0 0 0];



end