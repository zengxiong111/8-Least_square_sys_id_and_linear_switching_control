function [U_single,Y_single] = single_trajectory_generation(N,T,A,B,C,D,sigma_u,sigma_w,sigma_z)

n = size(A,1);
p = size(B,2);
m = size(C,1);

U_single = mvnrnd(zeros(p,1),sigma_u*eye(p),N+T-1);
U_single = U_single';
W = mvnrnd(zeros(n,1),sigma_w*eye(n),N+T-1);
W = W';
Z = mvnrnd(zeros(m,1),sigma_z*eye(m),N+T-1);
Z = Z';

Y_single = zeros(m,N+T-1);
X = zeros(n,N+T-1);

%X(1,1:n) = normrnd(zeros(n,1),sigma_x0*eye(n),n,1)';
X(1,1:n) = zeros(n,1);

for i=2:N+T-1
    X(:,i) = A * X(:,i-1) + B * U_single(:,i-1) + W(:,i-1);
    Y_single(:,i-1) = C * X(:,i-1) + D * U_single(:,i-1) + Z(:,i-1);
end


 

end