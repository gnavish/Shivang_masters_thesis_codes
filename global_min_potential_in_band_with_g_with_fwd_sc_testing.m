a=1;
t=0.1;
%mu=-4*t-0.02;
mu=-0.1;
N_x=400;
N_y=400;
%g_real_space=0.26;
%g_real_space=0.09;
g_real_space_min=0.06;
g_real_space_max=0.12;
num_g=20;
step_g=(g_real_space_max-g_real_space_min)/num_g;

%g=g_real_space/((N_x*N_y)^0.5);                                                                        %interaction strength
%pi=3.14;
k_0=1/pi;
delta_x_max=0.03;
delta_x_min=0.000001;
tol=10^(-5);
opt_delta_prev=10000;

delta_y_max=0.0005;
n_delta=100;
step_delta_x=(delta_x_max-(0))/n_delta;
step_delta_y=(delta_y_max-(0))/n_delta;
matrix_E=zeros(n_delta,n_delta);
matrix_E_diff_way=zeros(n_delta,n_delta);
% matrix_v_k=zeros(N_x,N_y);
matrix_E_k=zeros(N_x,N_y);
sum_e_k=0;
sum_1_by_E_l_diff_deltas=zeros(n_delta,1);
sum_nonint_piece=zeros(n_delta,1);
sum_int_piece=zeros(n_delta,1);
v_k_one_delta=zeros(N_x,N_y);
v_k_one_delta_diff_way=zeros(N_x,N_y);
diff_v_k_calc_diff_ways=zeros(N_x,N_y);
e_k_times_v_k_2_one_delta=zeros(N_x,N_y);
digits(40);
opt_delta_diff_g=zeros(num_g,1);

for qq=1:num_g

g_real_space=g_real_space_min+(qq-1)*step_g;
opt_delta_prev=10000;

for pp=1:10000


for mm=1:n_delta
    delta_x=delta_x_min+(mm-1)*step_delta_x;
    for nn=1:1%n_delta
    delta_y=0+(nn-1)*step_delta_y;
  E=0;
  E_diff_way=0;
sum_1_by_E_l=0;
sum_v_k_2=0;

for i=1:N_x


    for j=1:N_y

        % l_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        % l_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
         l_x=(0)+(2*pi/(N_x*a))*(i-1);
        l_y=(0)+(2*pi/(N_y*a))*(j-1);
        e_l=-2*t*(cos(l_x*a)+cos(l_y*a))-mu;
        R_l=((e_l^2+(abs(delta_x))^2)^(1/2));
        sum_1_by_E_l=(sum_1_by_E_l+(1/R_l));


    end

end
sum_1_by_E_l_diff_deltas(mm,1)=(sum_1_by_E_l);


for i=1:N_x
    for j=1:N_y
        
        % k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        % k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        k_x=(0)+(2*pi/(N_x*a))*(i-1);
        k_y=(0)+(2*pi/(N_y*a))*(j-1);
         e_k=(-2*t*(cos(k_x*a)+cos(k_y*a))-mu);
%         e_k=calc_e(k_x,k_y,a,t,mu);                                         %epsilon_k
        delta_k=delta_x;
%         delta_k=calc_delta(delta_x,delta_y,k_x,k_y);                        %delta_k

%         temp1=vpa(((e_k^2+(abs(delta_k))^2)^0.5));
%         temp2=vpa(e_k/temp1);
%         v_k=vpa((0.5*(1-temp2))^0.5);
         v_k=((1/2*(1-e_k/((e_k^2+(abs(delta_k))^2)^(1/2))))^(1/2));
        v_k_diff_way=0.5*delta_k*(1/e_k);
%         v_k=calc_mod_v(k_x,k_y,e_k,delta_k);                                %v_k
%         matrix_v_k(i,j)=v_k;
        u_k=(1/2*(1+e_k/((e_k^2+(abs(delta_k))^2)^(1/2))))^(1/2);
%         u_k=calc_mod_u(k_x,k_y,e_k,delta_k);                                %u_k
         R_k=(e_k^2+(abs(delta_k))^2)^(1/2);

         E=(E+2*(e_k)*(v_k^2));%-(g_real_space/(2*N_x*N_y))*((delta_x^2)/(R_k))*(sum_1_by_E_l);
         sum_v_k_2=sum_v_k_2+(v_k^2);
%             E=E+T45;
        E_diff_way=E_diff_way+2*e_k*(v_k_diff_way^2);

        
        if(mm==n_delta & nn==1)

%               matrix_E_k(i,j)=2*(e_k)*(v_k^2);%-(g_real_space/(2*N_x*N_y))*((delta_x^2)/(R_k))*(sum_1_by_E_l);
%             sum_e_k=sum_e_k+2*e_k*(v_k^2);

        end

        if(mm==657 & nn==1)

%             v_k_one_delta(i,j)=v_k;
%             v_k_one_delta_diff_way(i,j)=0.5*(delta_k)*(1/e_k);
%             e_k_times_v_k_2_one_delta(i,j)=e_k*(v_k^2);

        end
            
    end
end
sum_nonint_piece(mm,1)=E;
sum_int_piece(mm,1)=-(g_real_space/(2*N_x*N_y))*(sum_1_by_E_l^2)*(delta_x^2);
%matrix_E(mm,nn)=abs(E);
matrix_E(mm,nn)=(E-(g_real_space/(2*N_x*N_y))*(sum_1_by_E_l^2)*(delta_x^2)-(2*g_real_space/(N_x*N_y))*(sum_v_k_2^2));
matrix_E_diff_way(mm,nn)=E_diff_way-(g_real_space/(2*N_x*N_y))*(sum_1_by_E_l^2)*(delta_x^2);

    end
end
%  [delta_x_axis,delta_y_axis]=meshgrid((0+step_delta_x):step_delta_x:delta_x_max);
%   surf(delta_x_axis,delta_y_axis,matrix_E)
matrix_E_1=matrix_E(:,1);
matrix_E_1_diff_way=matrix_E_diff_way(:,1);

[min_E,opt_delta_index]=min(matrix_E_1);

opt_delta=delta_x_min+(opt_delta_index-1)*step_delta_x;
if(abs(opt_delta-opt_delta_prev)<tol)

    opt_delta_diff_g(qq,1)=opt_delta;
    break;

end

delta_x_min=opt_delta-2*step_delta_x;
delta_x_max=opt_delta+7*step_delta_x;
opt_delta_prev=opt_delta;

end



end


g_real_space_x_axis=g_real_space_min:step_g:g_real_space_max-step_g;
g_real_space_x_axis_dimensionless=g_real_space_x_axis/(8*t);
opt_delta_diff_g_dimensionless=opt_delta_diff_g/(8*t);
% plot(g_real_space_x_axis_dimensionless,opt_delta_diff_g_dimensionless)
% xlabel('Interaction strength/Bandwidth')
% ylabel('Variational solution/Bandwidth')
% print -deps Variational_sol_with_fwd_sc_terms_for_thesis

%dimensionfull figure
% figure
% scatter(g_real_space_x_axis,opt_delta_diff_g)
% figure
% one_upon_g_real_space_x_axis=1./g_real_space_x_axis;
% log_opt_delta_diff_g=log(opt_delta_diff_g);
% scatter(fliplr(one_upon_g_real_space_x_axis),transpose(fliplr(transpose(log_opt_delta_diff_g))));


%dimensionless figure
one_upon_g_real_space_x_axis_dimensionless=(8*t)*one_upon_g_real_space_x_axis;
log_opt_delta_diff_g_dimensionless_wfsc=log(opt_delta_diff_g/(8*t));
scatter(fliplr(one_upon_g_real_space_x_axis_dimensionless),transpose(fliplr(transpose(log_opt_delta_diff_g_dimensionless_wfsc))));
xlabel('Bandwidth/|g|');
ylabel('log(\Delta/Bandwidth)');
legend('BCS solution','Variational solution without forward scattering terms','Variational solution with forward scattering terms');

txt={'Lattice: 400 \times 400'}
text(6.5,-4.5,txt)
txt={'\nu=0.38'}
text(6.5,-4.8,txt)




% delta_x_axis=(0+step_delta_x):step_delta_x:delta_x_max;
% figure
% plot(delta_x_axis,matrix_E_1)
% title(['free electron band between -0.4 to +0.4  and chemical potential at ',num2str(mu)]);
% xlabel('delta');
% ylabel('Mean field energy functional value');
% [minima,index_minima]=min(matrix_E_1);
% delta_at_minima=index_minima*step_delta_x;
% 
% 
% 
% 
% 
% 
% % figure
% % surf(matrix_E_k)
% % plot(delta_x_axis,sum_1_by_E_l_diff_deltas)
% % figure
% % plot(delta_x_axis,sum_nonint_piece);
% % figure
% % plot(delta_x_axis,sum_int_piece);
% % figure
% % surf(v_k_one_delta);
% % figure
% % surf(v_k_one_delta_diff_way);
% % diff_v_k_calc_diff_ways=v_k_one_delta_diff_way-v_k_one_delta;
% % figure
% % surf(diff_v_k_calc_diff_ways)
% % figure
% % surf(e_k_times_v_k_2_one_delta)
% % figure
% % plot(delta_x_axis,matrix_E_1_diff_way);
% % sum_e_k
% 
% % extrema=find(imregionalmax(matrix_E));
% 
