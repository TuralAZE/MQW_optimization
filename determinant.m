
%%
determin = zeros(length(first_well),length(first_barrier),length(second_well),length(second_barrier),length(third_well));
n=4;
u=1;

Isolation_barrier = 50;
first_well = [35:1:50];
first_barrier = [10:1:30];
second_well = [20:1:40];
second_barrier = [15:1:40];
third_well = [10:1:20];

num_u = 688;

first_w = zeros(num_u,1);
second_w = zeros(num_u,1);
third_w = zeros(num_u,1);
first_b =  zeros(num_u,1);
second_b = zeros(num_u,1);

zzz_neg_4V = zeros(num_u,1);
zzz_0V = zeros(num_u,1);
zzz_4V = zeros(num_u,1);
dE_neg_4V = zeros(num_u,1);
dE_4V = zeros(num_u,1);
E_0V = zeros(num_u,1);


for f_w = 1:length(first_well)
    for f_b = 1:length(first_barrier)
        for s_w = 1:length(second_well)
            for s_b = 1:length(second_barrier)
                for t_w = 1:length(third_well)
                    
                    if  ( E(f_w,f_b,s_w,s_b,t_w,2,1)<0.145 && E(f_w,f_b,s_w,s_b,t_w,2,1)>0.125 )...
                            && ( E(f_w,f_b,s_w,s_b,t_w,3,2)<0.145 && E(f_w,f_b,s_w,s_b,t_w,3,2)>0.125 )...
                            && ( abs(E(f_w,f_b,s_w,s_b,t_w,3,2) - E(f_w,f_b,s_w,s_b,t_w,2,1))<0.002 )...
                            && ( (Z(f_w,f_b,s_w,s_b,t_w,2,2) - Z(f_w,f_b,s_w,s_b,t_w,1,1)) > 0 )...
                            && ( (Z(f_w,f_b,s_w,s_b,t_w,3,3) - Z(f_w,f_b,s_w,s_b,t_w,2,2)) > 0 )...
                            && ( (E(f_w,f_b,s_w,s_b,t_w,4,3) - 1*E(f_w,f_b,s_w,s_b,t_w,2,1)) > 0 )...
                            && ( abs(Z(f_w,f_b,s_w,s_b,t_w,3,3) + Z(f_w,f_b,s_w,s_b,t_w,1,1) - 2*Z(f_w,f_b,s_w,s_b,t_w,2,2)) < abs(Z(f_w,f_b,s_w,s_b,t_w,2,2) - Z(f_w,f_b,s_w,s_b,t_w,1,1))/10 )
                        determin(f_w,f_b,s_w,s_b,t_w)=1;
                        

                  
                        
                        dipole = zeros(4,4);
                        energy = zeros(4,4);
                        dipole(:,:) = Z(f_w,f_b,s_w,s_b,t_w,:,:);
                        energy(:,:) = E(f_w,f_b,s_w,s_b,t_w,:,:);
                        [data1,data2,data3,data4,data5,data6]=perturbation(dipole,energy,n);
                        
                        if ((abs(data2) > 1E-27)...
                            && ( E(f_w,f_b,s_w,s_b,t_w,3,2)<0.145 && E(f_w,f_b,s_w,s_b,t_w,3,2)>0.125 )...
                            && ( abs(E(f_w,f_b,s_w,s_b,t_w,3,2) - E(f_w,f_b,s_w,s_b,t_w,2,1))<0.002 ))
                        first_w(u,1) = f_w;
                        second_w(u,1) = s_w;
                        third_w(u,1) = t_w;
                        first_b(u,1) = f_b;
                        second_b(u,1) = s_b;
                        
                        zzz_neg_4V(u,1) = data1;
                        zzz_0V(u,1) = data2;
                        zzz_4V(u,1) = data3;
                        dE_neg_4V(u,1) = data4;
                        dE_4V(u,1) = data5;
                        E_0V(u,1) = data6;
                        u=u+1;
                        end
                    end
                end
            end
        end
    end
end

u

