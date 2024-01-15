

a=1;
t=0.1;
mu=-4*t-(0.02);
%mu=-4*t-0.36;
N_x=80;
N_y=80;
%mod_g_real_space=12.2;                                                      %interaction strength
mod_g_real_space=11.6;
% pi=3.14;
% k_0=1/pi;
delta_x_max=0.45;
delta_y_max=0.45;
% delta_x_max=0.2;
% delta_y_max=0.2;
delta_x_min=0;
delta_y_min=0;
n_delta=40;
step_delta_x=(delta_x_max-(delta_x_min))/n_delta;
step_delta_y=(delta_y_max-(delta_y_min))/n_delta;
matrix_E=zeros(n_delta,n_delta);
gamma=1;
digits(20);
N_particles_array=zeros(n_delta,n_delta);
n_iterations=20;
MFE_iterations=zeros(n_iterations,1);
M1=zeros(1000,20);
M2=zeros(N_x,N_y);
M3=zeros(N_x,N_y);
M4=zeros(N_x,N_y);
a1=(1/1.5)^2;
b=0.5;
c=(1.5)^0.5;
B=1.5;
sigma12=1.3;
sigma22=(sigma12/2);
o=1;
T1_imag=0;
tol=10^(-5);

V_k_space_1d=zeros(N_x,N_y);
num_sum_p=9;


V_p_r=zeros(N_x,N_y);


x_pos_array=(1:1:N_x)';
y_pos_array=1:1:N_y;

% for ii=-(num_sum_p-1)/2:(num_sum_p-1)/2
%     for jj=-(num_sum_p-1)/2:(num_sum_p-1)/2
% 
%         V_p_r=V_p_r - exp( -( (x_pos_array-1+ii*N_x).^2 + (y_pos_array-1+jj*N_y).^2 )/sigma22 );
% 
%     end
% 
% end


for ii=-(num_sum_p-1)/2:(num_sum_p-1)/2
    for jj=-(num_sum_p-1)/2:(num_sum_p-1)/2

        V_p_r=V_p_r+(B-1)*exp( -( (x_pos_array-1+ii*N_x).^2 + (y_pos_array-1+jj*N_y).^2 )/sigma12 ) - B*exp( -( (x_pos_array-1+ii*N_x).^2 + (y_pos_array-1+jj*N_y).^2 )/sigma22 );

    end

end
V_p_r=(1/(num_sum_p^0.5))*V_p_r;
V_p_r(1,1)=0;
surf(V_p_r*mod_g_real_space);
xlabel('Delta_x');
ylabel('Delta_y');
zlabel('potential');






V_k_space_1d=V_k_space_1d*(mod_g_real_space)/((N_x*N_y)^0.5);
sum_V_k=0;


sum_V_k=sum(sum(V_k_space_1d));


















tic
k_x=(0:(2*pi/N_x):(2*pi/N_x)*(N_x-1))';
    k_y=(0:(2*pi/N_y):(2*pi/N_y)*(N_y-1))';
    k_array=[ repelem(k_x,numel(k_y)) repmat(k_y,numel(k_x),1) ];

    X_array=(1:1:N_x)';
    Y_array=(1:1:N_y)';
    R_array=([ repelem(X_array,numel(Y_array)) repmat(Y_array,numel(X_array),1) ])';
    R_array_t=(R_array)';
    V_r=zeros(N_x*N_y,1);
    for ii=-(num_sum_p-1)/2:(num_sum_p-1)/2
        for jj=-(num_sum_p-1)/2:(num_sum_p-1)/2
    
            V_r=V_r+(B-1)*exp( -( (R_array_t(:,1)-1+ii*N_x).^2 + (R_array_t(:,2)-1+jj*N_y).^2 )/sigma12 ) - B*exp( -( (R_array_t(:,1)-1+ii*N_x).^2 + (R_array_t(:,2)-1+jj*N_y).^2)/sigma22 );
    
        end
    end
      V_r(1,1)=0;
      V_r=(1/(num_sum_p^0.5))*V_r;
    ft_matrix=zeros(N_x*N_y,N_x*N_y);
    ft_matrix=exp( 1i*k_array(:,1)*(R_array(1,:)-1)+1i*k_array(:,2)*(R_array(2,:)-1) );
    V_k_1d=ft_matrix*V_r;
    V_k_1d_t=(V_k_1d)';
    sum_V_k=sum(V_k_1d);

    e_k=-2*t*( cos(k_array(:,1)) +cos(k_array(:,2)) )-mu;
    l_array=(k_array)';
    V_k_l=zeros(N_x*N_y,N_x*N_y);
    
    for ii=1:N_x
        temp=zeros(N_x*N_y,1);
        temp=circshift(V_k_1d,(ii-1)*N_y);
%         if(ii==1)
% 
%             plot(real(temp)*mod_g_real_space/((N_x*N_y)^0.5));
%             hold on
% 
%         end
%         if(ii~=1)
% 
%             temp(1:(ii-1)*N_y,1)=flip(temp(1:(ii-1)*N_y,1));
% 
%         end

            for jj=1:N_y
                
                for zz=1:N_x
    
                    V_k_l((((zz-1)*N_y+1):(zz*N_y)),(ii-1)*N_y+jj)=circshift( temp(((zz-1)*N_y+1):(zz*N_y),1), jj-1);
                    
                
                end
            end

    end
  




%     V_k_l1=zeros(N_x*N_y,N_x*N_y);
%     for x=1:N_x
% 
%         for y=1:N_y
%     
%  %           V_k_space_1d=V_k_space_1d+V_p_r(x,y)*exp(1i*(q_x*(x-1)+q_y*(y-1)));
%             V_k_l1=V_k_l1+V_p_r(x,y)*exp(1i*( (mod(k_array(:,1)-l_array(1,:),2*pi))*(x-1) + (mod(k_array(:,2)-l_array(2,:),2*pi))*(y-1) ));
% 
%         end
% 
%     end
% 
%    
% 

toc

V_k_l=V_k_l*(mod_g_real_space)/((N_x*N_y)^0.5);
%  V_k_l1=V_k_l1*(mod_g_real_space)/((N_x*N_y)^0.5);
%  V_k_l=V_k_l1;
% plot(real(V_k_l1(:,1)))
tic
for mm=1:n_delta
    delta_x=delta_x_min+(mm-1)*step_delta_x;

    for nn=1:n_delta%n_delta
    delta_y=delta_y_min+(nn-1)*step_delta_y;

    E=0;
    T1=0;
    T2=0;
    T3_temp=0;
    T3=0;
    T4=0;
    T11=0;
    T21=0;
    T31=0;
    N_particles=0;
    count=0;
   
    
    delta_k=delta_x*sin(k_array(:,1))-1i*delta_y*sin(k_array(:,2));
    R_k=(e_k.^2+(abs(delta_k)).^2).^0.5;
    K1=1-e_k./R_k;
    N_particles_array(mm,nn)=0.5*sum(K1);
    K2=1+e_k./R_k;
    K3=delta_k./R_k;
    E=0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3 + ((K1)')*V_k_l*K2 + V_k_l(1,1)*(sum(sum(K1)))^2);
    matrix_E(mm,nn)=E;
   



    end

end
toc
delta_x_axis=delta_x_min:step_delta_x:(delta_x_max-step_delta_x);
delta_y_axis=delta_y_min:step_delta_y:(delta_y_max-step_delta_y);
%delta_xy_axis=sqrt(delta_x_axis)
delta_xy_axis=(2^0.5)*delta_x_axis;
figure
surf(delta_x_axis,delta_y_axis,real(matrix_E));
title(['free electron band between -4 to +4  and chemical potential at ',num2str(mu)]);
xlabel('delta_x');
ylabel('delta_y');
zlabel('mean-field energy');

matrix_E_xslice=real(matrix_E(:,1));
matrix_E_yslice=real(matrix_E(1,:));
matrix_E_diagslice=real(diag(matrix_E));

delta_x_axis_dimensionless=delta_x_axis/(8*t);
delta_xy_axis_dimensionless=delta_xy_axis/(8*t);
figure
scatter(delta_x_axis_dimensionless,(matrix_E_xslice/(8*t)))
hold on
scatter(delta_x_axis_dimensionless,(matrix_E_yslice/(8*t)),'x')
hold on
scatter(delta_xy_axis_dimensionless(1:floor(n_delta/(2^0.5))+1),(matrix_E_diagslice(1:floor(n_delta/(2^0.5))+1)/(8*t)),'+')

%scatter(delta_xy_axis_dimensionless,(matrix_E_diagslice/(8*t)),'+')
legend('\Delta_y=0 slice','\Delta_x=0 slice','\Delta_x=\Delta_y slice')
xlabel('|\Delta_x+\iota\Delta_y|/Bandwidth')
ylabel('Mean-field energy/Bandwidth')







