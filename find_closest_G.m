function [A,B,C,D] = system_generation(r,m,n,p)

A = zeros(n,n);

B = zeros(n,p);
C = zeros(m,n);

D = zeros(m,p);

B = normrnd(0,1/n,[n,p]);
C = normrnd(0,1/m,[m,n]);

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