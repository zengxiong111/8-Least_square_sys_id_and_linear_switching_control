function [epsilon_a,epsilon_c] = compute_epsilon(A_t_all,B_t_all,C_t_all)

 N = size(A_t_all,1);
 n = size(A_t_all,3);

 epsilon_a = 1e16;

 for i = 1:N
     for j =1:N
         if(vrho(A_t_all(i,j,:,:))>1)
             if vrho(A_t_all(i,j,:,:))-1 < epsilon_a
                 epsilon_a = vrho(A_t_all(i,j,:,:)) - 1;
             end
         end
     end
 end

 for i = 1:N
     for j =1:N 
         [V,D] = eig(A_t_all(i,j,:,:));
         for k=1:n

         end
     end
 end
