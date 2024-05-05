function [A,B,C,D] = system_generation(r,m,n,p)

for i = 1:5
    A = zeros(n,n);
    B = zeros(n,p);
    C = zeros(m,n);
    D = zeros(m,p);
    
    B = normrnd(0,1/n,[n,p]);
    C = normrnd(0,1/m,[m,n]);
    A = rand(n);
    A = r * A/vrho(A);

    Ob = obsv(A,C);
    Co = ctrb(A,B);
    % check if the genrated system is controllable and observable
    if(rank(Ob) == n & rank(Co) == n)
        break;
    end
    if(rank(Ob) < n | rank(Co) < n)
        if(i==5)
            disp('system generation fails!');
        end
        continue;
    end

end

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

% A = zeros(n,n);
%     B = zeros(n,1);
%     C = zeros(1,n);
%     A(2:n-1,3:n) = v * eye(n-2);
%     A(1,1) = rho;
%     A(1,2) = v;
%     B(n,1) = v;
%     B(1,1) = b1;
%     C(1,1) = 1;


 

end