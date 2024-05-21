function [U_single,Y_single,X] = single_trajectory_generation(N,A,B,C,D,sigma_u,sigma_w,sigma_z,X_1)

n = size(A,1);
p = size(B,2);
m = size(C,1);

U_single = mvnrnd(zeros(p,1),sigma_u*eye(p),N);
U_single = U_single';
W = mvnrnd(zeros(n,1),sigma_w*eye(n),N);
W = W';
Z = mvnrnd(zeros(m,1),sigma_z*eye(m),N);
Z = Z';

Y_single = zeros(m,N);
X = zeros(n,N);

%X(1,1:n) = normrnd(zeros(n,1),sigma_x0*eye(n),n,1)';
%X(1,1:n) = zeros(n,1);
X(1,1:n) = X_1;

for i=2:N
    X(:,i) = A * X(:,i-1) + B * U_single(:,i-1) + W(:,i-1);
    Y_single(:,i-1) = C * X(:,i-1) + D * U_single(:,i-1) + Z(:,i-1);
end


 

end