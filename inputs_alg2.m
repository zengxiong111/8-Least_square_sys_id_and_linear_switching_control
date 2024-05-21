function [M_all,tau_all,Xi_all] = inputs_alg2(A_t_all,B_t_all,C_t_all,epsilon_a,epsilon_c,sigma_w_2,sigma_u_2,sigma_v_2,delta_p)

m_a=0;
m_s=0;
m_p=0;
m_t=0;
c_s=0;

n = size(A_t_all,1);

N = size(A_t_all,3);

%For m_p,m_t,c_s, we need to check the stability of the closed-loop.
for i=1:N
    for j=1:N
        m_a_temp = max([1,norm(A_t_all(:,:,i,j))]);
        m_s_temp = max([1,norm(B_t_all(:,:,i,j)), norm(C_t_all(:,:,i,j))]);
        if(m_a_temp > m_a)
            m_a = m_a_temp;
        end
        if(m_s_temp > m_s)
            m_s = m_s_temp;
        end
    end
end

for i=1:N
    for j=1:N
        if(vrho(A_t_all(:,:,i,j))<1)
            P = dlyap(A_t_all(:,:,i,j),C_t_all(:,:,i,j)'*C_t_all(:,:,i,j));
            m_p_temp = norm(P);
            m_t_temp = trace(P);
            c_s_temp = 5*(sigma_w_2*trace(P)+sigma_u_2*...
                trace(B_t_all(:,:,i,j)'*P*B_t_all(:,:,i,j))+...
                sigma_v_2*n)*log(1/delta_p);
            if(m_p_temp > m_p)
                m_p = m_p_temp;
            end
            if(m_t_temp > m_t)
                m_t = m_t_temp;
            end
            if(c_s_temp > c_s)
                c_s = c_s_temp;
            end
        end
    end
end

sigma_m = max([sigma_w_2,sigma_u_2,sigma_v_2]);
c_e = max([1,1/log(1+epsilon_a)]);
c_r = m_p * (22*epsilon_c^(-2)+1)*sigma_m^2*c_e;
c_p = 2*max([1,m_p/(epsilon_c^2)]);

M_all = zeros(N,1);

M_all(1) = 5 * m_p * log(1/delta_p);

%tau_all(1) = n + log(1/delta_p);
tau_all(1) = int16(max([1600/9*log(1/delta_p),...
   1/log(1+epsilon_a)*log(6400*epsilon_c/(9*sigma_w_2*delta_p)*(M_all(1) + ...
   c_r *log(2/delta_p)*n^2 * m_a^(4*n) + c_s * log(1/delta_p) ))]));
for j= 2:N
    tau_all(j) = int16((j-1)*(2/log(m_a)*n + log(c_p))) + tau_all(1);
end

for j=2:N
    M_all(j) = m_p/(epsilon_c^2)*m_a^(2*n)*(M_all(j-1)+double(tau_all(j-1))*c_s) + c_r * log(2/delta_p)* n^2 * m_a^(4*n);
end

Xi_all = zeros(N,1);

for k=1:N
    Xi_all(k) = 0;
    for i=1:N
        for j=1:N
            if(vrho(A_t_all(:,:,i,j))<1)
                P = dlyap(A_t_all(:,:,1,1),C_t_all(:,:,1,1)'*C_t_all(:,:,1,1));
                Xi_all_temp = 2*M_all(k)+10*tau_all(k)*log(1/delta_p)*...
                    (sigma_w_2*trace(P)+sigma_u_2*trace(B_t_all( :,:,1,1)'*P*B_t_all( :,:,1,1))+sigma_v_2*n);
                if(Xi_all_temp > Xi_all(k))
                    Xi_all(k) = Xi_all_temp;
                end
            end
        end
    end
end