function [A,B,C,D] = system_generation(r,m,n,p)

G = zeros(m,T*p);
G(:,1:p) = D;
for i = 1:T-1
    G(:,i*p + 1: (i+1)*p) = C*A^(i-1)*B;
end

 

end