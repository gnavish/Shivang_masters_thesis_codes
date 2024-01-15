a=1;
t=0.1;
mu=-0.1;
N_x=400;
N_y=400;
%g_real_space=0.08;
g_real_space_max=0.12;
g_real_space_min=0.06;
n_g=20;
step_g=(g_real_space_max-g_real_space_min)/n_g;
BCS_sol_diff_g=zeros(n_g,1);
% g=g_real_space/((N_x*N_y)^0.5);     %interaction strength
debye_energy=0.1;
tol=10^(-3);
delta_x_max=10;
delta_x_min=10^(-5);
LHS_min=0;
LHS_max=0;
LHS_mid=0;


%pi=3.14;
k_0=1/pi;
%delta_x_max=0.002;
delta_y_max=0.002;
n_delta=1000;
step_delta_x=(delta_x_max-(0))/n_delta;
step_delta_y=(delta_y_max-(0))/n_delta;
matrix_E=zeros(n_delta,n_delta);
LHS_diff_deltas=zeros(n_delta,1);
delta_BCS_self_consis=0;
LHS_for_one_delta_at_diff_k=zeros(N_x,N_y);
e_k_at_debye_freq=zeros(N_x,N_y);

for qq=1:n_g
g_real_space=g_real_space_min+(qq-1)*step_g;
delta_x_max=10;
delta_x_min=10^(-5);
LHS_min=0;
LHS_max=0;
LHS_mid=0;




for pp=1:10000

%calc for delta_x_min
delta_x=delta_x_min;

%for mm=1:n_delta
%    delta_x=0+(mm-1)*step_delta_x;
 %   for nn=1:1%n_delta
    
%    delta_y=0+(nn-1)*step_delta_y;
    delta_y=0;
    E=0;
    LHS=0;
for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_k=calc_e(k_x,k_y,a,t,mu);                                         %epsilon_k
        
        delta_k=calc_delta(delta_x,delta_y,k_x,k_y);                        %delta_
        mm=1;
        if(mm==2 & nn==1)
%              k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
%         k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1)
%         e_k=calc_e(k_x,k_y,a,t,mu)                                         %epsilon_k
%        
%        delta_k=calc_delta(delta_x,delta_y,k_x,k_y)                        %delta_k
            LHS_for_one_delta_at_diff_k(i,j)=(1/((delta_k^2)+(e_k^2))^0.5);%*(g_real_space/(2*N_x*N_y));

        end



%          if(abs(e_k)<debye_energy)

             e_k_at_debye_freq(i,j)=e_k;
            LHS_min=LHS_min+(1/((delta_k^2)+(e_k^2))^0.5);
%          end
    end
end


%----------------------------------------------------------
%calc for delta_x_max
delta_x=delta_x_max;

%for mm=1:n_delta
%    delta_x=0+(mm-1)*step_delta_x;
 %   for nn=1:1%n_delta
%    delta_y=0+(nn-1)*step_delta_y;
    E=0;
    LHS=0;
for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_k=calc_e(k_x,k_y,a,t,mu);                                         %epsilon_k
        
        delta_k=calc_delta(delta_x,delta_y,k_x,k_y);                        %delta_
        if(mm==2 & nn==1)
%              k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
%         k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1)
%         e_k=calc_e(k_x,k_y,a,t,mu)                                         %epsilon_k
%        
%        delta_k=calc_delta(delta_x,delta_y,k_x,k_y)                        %delta_k
            LHS_for_one_delta_at_diff_k(i,j)=(1/((delta_k^2)+(e_k^2))^0.5);%*(g_real_space/(2*N_x*N_y));

        end



%          if(abs(e_k)<debye_energy)

             e_k_at_debye_freq(i,j)=e_k;
            LHS_max=LHS_max+(1/((delta_k^2)+(e_k^2))^0.5);
%          end
    end
end


%calc for mid pt

delta_x=(delta_x_min+delta_x_max)/2;
%for mm=1:n_delta
%    delta_x=0+(mm-1)*step_delta_x;
 %   for nn=1:1%n_delta
%    delta_y=0+(nn-1)*step_delta_y;
    E=0;
    LHS=0;
for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_k=calc_e(k_x,k_y,a,t,mu);                                         %epsilon_k
        
        delta_k=calc_delta(delta_x,delta_y,k_x,k_y);                        %delta_
        if(mm==2 & nn==1)
%              k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
%         k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1)
%         e_k=calc_e(k_x,k_y,a,t,mu)                                         %epsilon_k
%        
%        delta_k=calc_delta(delta_x,delta_y,k_x,k_y)                        %delta_k
            LHS_for_one_delta_at_diff_k(i,j)=(1/((delta_k^2)+(e_k^2))^0.5);%*(g_real_space/(2*N_x*N_y));

        end



%          if(abs(e_k)<debye_energy)

             e_k_at_debye_freq(i,j)=e_k;
            LHS_mid=LHS_mid+(1/((delta_k^2)+(e_k^2))^0.5);
%          end
    end
end


LHS_max=LHS_max*(g_real_space/(N_x*N_y));
LHS_min=LHS_min*(g_real_space/(N_x*N_y));
LHS_mid=LHS_mid*(g_real_space/(N_x*N_y));



if(abs(LHS_mid-1)<tol)

BCS_sol_diff_g(qq,1)=(delta_x_min+delta_x_max)/2;


break;

end

if((LHS_mid-1)>tol)

    delta_x_min=(delta_x_min+delta_x_max)/2;

end

if((1-LHS_mid)>tol)

    delta_x_max=(delta_x_min+delta_x_max)/2;

end

end


end


g_real_space_x_axis_dimensionless=g_real_space_min:step_g:g_real_space_max-step_g;
g_real_space_x_axis_dimensionless=g_real_space_x_axis_dimensionless./(8*t);
BCS_sol_diff_g_dimensionless=BCS_sol_diff_g./(8*t);
figure
scatter(g_real_space_x_axis_dimensionless,BCS_sol_diff_g_dimensionless);
hold on



% 
%   %  end
%     LHS_diff_deltas(mm,1)=LHS*(g_real_space/(N_x*N_y));
% %end
% delta_x_axis=(0+step_delta_x):step_delta_x:delta_x_max;
% 
% for i=1:(n_delta-1)
% 
%     if(LHS_diff_deltas(i)<1 & LHS_diff_deltas(i+1)>1)
% 
%         delta_BCS_self_consis=delta_x_axis(i);
% 
%     end
% 
%     if(LHS_diff_deltas(i)>1 & LHS_diff_deltas(i+1)<1)
% 
%         delta_BCS_self_consis=delta_x_axis(i);
% 
%     end
% 
% end
% delta_BCS_self_consis
% % delta_BCS_self_consis_diff_g(u)=delta_BCS_self_consis;
% matrix_ones=ones(n_delta,1);
% % g_x_axis=(min_g_real_space+step_g_real_space):step_g_real_space:(min_g_real_space+step_g_real_space*num_g);
% % plot(g_x_axis,delta_BCS_self_consis_diff_g);
% 
% figure
% plot(delta_x_axis,LHS_diff_deltas(:,1));
% hold on
% plot(delta_x_axis,matrix_ones)
% title(['BCS self consis. eqn sol for bands between -4 to +4  and chemical potential at ',num2str(mu)]);
% xlabel('delta');
% figure 
% surf(LHS_for_one_delta_at_diff_k);
% figure
% surf(e_k_at_debye_freq);




%-------------------------------------------------------------------------------------------------
a=1;
t=0.1;
mu=-0.1;
N_x=400;
N_y=400;
%g_real_space=0.05;
g_real_space_min=0.06;
g_real_space_max=0.12;
num_g=n_g;
step_g=(g_real_space_max-g_real_space_min)/num_g;




g=g_real_space/((N_x*N_y)^0.5);                                                                        %interaction strength
%pi=3.14;
k_0=1/pi;

delta_x_max=0.003;
delta_x_min=0.000001;
tol=10^(-5);
opt_delta_prev=10000;



%delta_x_max=0.0001;
delta_y_max=0.0001;
n_delta=1000;
step_delta_x=(delta_x_max-(0))/n_delta;
step_delta_y=(delta_y_max-(0))/n_delta;
matrix_E=zeros(n_delta,n_delta);
% matrix_v_k=zeros(N_x,N_y);
matrix_E_k=zeros(N_x,N_y);
sum_e_k=0;
opt_delta_diff_g=zeros(num_g,1);



for qq=1:num_g

g_real_space=g_real_space_min+(qq-1)*step_g;


for pp=1:10000





for mm=1:n_delta
    delta_x=delta_x_min+(mm-1)*step_delta_x;
    for nn=1:1%n_delta
    delta_y=0+(nn-1)*step_delta_y;
  E=0;
sum_1_by_E_l=0;

for i=1:N_x


    for j=1:N_y

        l_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        l_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_l=-2*t*[cos(l_x*a)+cos(l_y*a)]-mu;
        R_l=(e_l^2+(abs(delta_x))^2)^(1/2);
        sum_1_by_E_l=sum_1_by_E_l+(1/R_l);


    end

end

temp=0;

for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
         e_k=-2*t*[cos(k_x*a)+cos(k_y*a)]-mu;
%         e_k=calc_e(k_x,k_y,a,t,mu);                                         %epsilon_k
        delta_k=delta_x;
%         delta_k=calc_delta(delta_x,delta_y,k_x,k_y);                        %delta_k
        v_k=(1/2*(1-e_k/((e_k^2+(abs(delta_k))^2)^(1/2))))^(1/2);
%         v_k=calc_mod_v(k_x,k_y,e_k,delta_k);                                %v_k
%         matrix_v_k(i,j)=v_k;
        u_k=(1/2*(1+e_k/((e_k^2+(abs(delta_k))^2)^(1/2))))^(1/2);
%         u_k=calc_mod_u(k_x,k_y,e_k,delta_k);                                %u_k
         R_k=(e_k^2+(abs(delta_k))^2)^(1/2);

         E=E+2*(e_k)*(v_k^2)-(g_real_space/(2*N_x*N_y))*((delta_x^2)/(R_k))*(sum_1_by_E_l);
%             E=E+T45;
        temp=temp+(1-(e_k/R_k));
        if(mm==n_delta & nn==1)

            matrix_E_k(i,j)=2*(e_k)*(v_k^2)-(g_real_space/(2*N_x*N_y))*((delta_x^2)/(R_k))*(sum_1_by_E_l);
%             sum_e_k=sum_e_k+2*e_k*(v_k^2);

        end
            
    end
end

%matrix_E(mm,nn)=abs(E);
matrix_E(mm,nn)=(E);%-(g_real_space/(2*N_x*N_y))*(temp^2);
    end
end
%  [delta_x_axis,delta_y_axis]=meshgrid((0+step_delta_x):step_delta_x:delta_x_max);
%   surf(delta_x_axis,delta_y_axis,matrix_E)
matrix_E_1=matrix_E(:,1);

[min_E,opt_delta_index]=min(matrix_E_1);

opt_delta=delta_x_min+(opt_delta_index-1)*step_delta_x;

if(abs(opt_delta-opt_delta_prev)<tol)

    opt_delta_diff_g(qq,1)=opt_delta;
    break;

end

delta_x_min=opt_delta-step_delta_x;
delta_x_max=opt_delta+step_delta_x;
opt_delta_prev=opt_delta;

end



end


opt_delta_diff_g_dimensionless=opt_delta_diff_g./(8*t);
plot(g_real_space_x_axis_dimensionless,opt_delta_diff_g_dimensionless)
legend('BCS solution','Variational solution')
xlabel('Interaction strength/Bandwidth')
ylabel('BCS solution or optimal order parameter/Bandwidth')
print -deps BCS_sol_and_variational_sol_without_sc_terms_comparison1



%dimensionfull figures

figure

g_real_space_x_axis=g_real_space_min:step_g:g_real_space_max-step_g;
scatter(g_real_space_x_axis,BCS_sol_diff_g);
hold on
plot(g_real_space_x_axis,opt_delta_diff_g);


one_upon_g_real_space_x_axis=1./g_real_space_x_axis;
one_upon_g_real_space_x_axis_dimensionless=(8*t)*one_upon_g_real_space_x_axis;

log_BCS_sol_diff_g=log(BCS_sol_diff_g);
log_opt_delta_diff_g=log(opt_delta_diff_g);
log_opt_delta_diff_g_dimensionless=log(opt_delta_diff_g/(8*t));
log_BCS_sol_diff_g_dimensionless=log(BCS_sol_diff_g/(8*t));

figure
scatter(fliplr(one_upon_g_real_space_x_axis),transpose(fliplr(transpose(log_BCS_sol_diff_g))));
hold on
scatter(fliplr(one_upon_g_real_space_x_axis),transpose(fliplr(transpose(log_opt_delta_diff_g))),'x');

%dimensionless final figure
figure
scatter(fliplr(one_upon_g_real_space_x_axis_dimensionless),transpose(fliplr(transpose(log_BCS_sol_diff_g_dimensionless))));
hold on
scatter(fliplr(one_upon_g_real_space_x_axis_dimensionless),transpose(fliplr(transpose(log_opt_delta_diff_g_dimensionless))),'x');


hold on





% delta_x_axis=(0+step_delta_x):step_delta_x:delta_x_max;
% figure
% plot(delta_x_axis,matrix_E_1)
% title(['free electron band between -0.4 to +0.4  and chemical potential at ',num2str(mu)]);
% xlabel('delta');
% ylabel('Mean field energy functional value');
% [minima,index_minima]=min(matrix_E_1);
% delta_at_minima=index_minima*step_delta_x;
% figure
% surf(matrix_E_k)
% % sum_e_k
% 
% % extrema=find(imregionalmax(matrix_E));

