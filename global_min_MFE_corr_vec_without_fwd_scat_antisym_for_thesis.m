

a=1;
t=0.1;
mu=-4*t-(0.02);
N_x=80;
N_y=80;
mod_g_real_space=14;                                                      %interaction strength
% pi=3.14;
% k_0=1/pi;
delta_x_max=0.3;
delta_y_max=0.3;
delta_x_min=0;
delta_y_min=0;
n_delta=40;
step_delta_x=(delta_x_max-(delta_x_min))/n_delta;
step_delta_y=(delta_y_max-(delta_y_min))/n_delta;
matrix_E=zeros(n_delta,n_delta);
gamma=1;
digits(20);
N_particles_array=zeros(n_delta,n_delta);
n_iterations=80;
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
%surf(V_p_r*mod_g_real_space);
%xlabel('Delta_x');
%ylabel('Delta_y');
%zlabel('potential');






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
    E=0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3);
    matrix_E(mm,nn)=E;
   



    end

end
toc
delta_x_axis=delta_x_min:step_delta_x:(delta_x_max-step_delta_x);
delta_y_axis=delta_y_min:step_delta_y:(delta_y_max-step_delta_y);
%figure
%surf(delta_x_axis,delta_y_axis,real(matrix_E));
%title(['free electron band between -4 to +4  and chemical potential at ',num2str(mu)]);
%xlabel('delta_x');
%ylabel('delta_y');
%zlabel('mean-field energy');
%figure
%surf(delta_x_axis,delta_y_axis,N_particles_array);







sorted_matrix_E=sort(sort(real(matrix_E)),2);
min_E=sorted_matrix_E(1,1);
MFE_iterations(1)=min_E;
[min_x,min_y]=find(real(matrix_E)==min_E);
min_x=delta_x_min+(min_x-1)*step_delta_x
min_y=delta_y_min+(min_y-1)*step_delta_y

k_x=(0:(2*pi/N_x):(2*pi/N_x)*(N_x-1))';
k_y=(0:(2*pi/N_y):(2*pi/N_y)*(N_y-1))';
k_array=[ repelem(k_x,numel(k_y)) repmat(k_y,numel(k_x),1) ];
delta_k_matrix=min_x*sin(k_array(:,1))-1i*(min_y)*sin(k_array(:,2));
k_x_axis=0:(2*pi/(N_x)):(2*pi-(2*pi/(N_x)));
k_y_axis=0:(2*pi/(N_y)):(2*pi-(2*pi/(N_y)));

% figure
% abs_delta_k_matrix=reshape(abs(delta_k_matrix),N_x,N_y);
% temp1=circshift( abs_delta_k_matrix, N_y-floor( N_y/2+0.1 ), 2 );
% temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
% surf((k_y_axis),(k_x_axis),temp2);

delta_k_matrix_shifted=reshape(delta_k_matrix,N_x,N_y);
delta_k_matrix_shifted=circshift( delta_k_matrix_shifted, N_y-floor( N_y/2+0.1 ), 2 );
delta_k_matrix_shifted=circshift( delta_k_matrix_shifted, N_x-floor( N_x/2+0.1 ), 1 );
%figure
%surf(k_y_axis, k_x_axis, abs(delta_k_matrix_shifted))
%xlabel('k_x');
%ylabel('k_y');
%zlabel('mod(Delta_k)');
%figure
%surf(k_y_axis, k_x_axis, angle(delta_k_matrix_shifted))
%xlabel('k_x');
%ylabel('k_y');
%zlabel('phase(Delta_k)');
% figure
% surf(k_y_axis, k_x_axis, abs(delta_k_matrix_shifted-rot90(delta_k_matrix_shifted)))
% 




phase_delta_k=reshape(atan(imag(delta_k_matrix)./real(delta_k_matrix)),N_x,N_y);
temp1=circshift( phase_delta_k, N_y-floor( N_y/2+0.1 ), 2 );
temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );

%figure
%surf(k_y_axis,k_x_axis,reshape(phase_delta_k,N_x,N_y));



for uu=2:n_iterations


    %antisymmetrization of delta_k
    delta_k_matrix_reshaped_minusky=zeros(N_x,N_y);
    delta_k_matrix_reshaped_minusky_minuskx=zeros(N_x,N_y);


    delta_k_matrix_reshaped=reshape(delta_k_matrix,N_x,N_y);
    delta_k_matrix_reshaped_minusky(:,1)=delta_k_matrix_reshaped(:,1);
    delta_k_matrix_reshaped_minusky(:,2:N_y)=flip(delta_k_matrix_reshaped(:,2:N_y),2);
    delta_k_matrix_reshaped_minusky_minuskx(1,:)=delta_k_matrix_reshaped_minusky(1,:);
    delta_k_matrix_reshaped_minusky_minuskx(2:N_x,:)=flip(delta_k_matrix_reshaped_minusky(2:N_x,:),1);
    delta_k_matrix_reshaped_antisymmetrized=0.5*(delta_k_matrix_reshaped-delta_k_matrix_reshaped_minusky_minuskx);
    delta_k_matrix_antisymmetrized=reshape(delta_k_matrix_reshaped_antisymmetrized,N_x*N_y,1);
    %delta_k_matrix_anmtisymmetrised is the antisymmetrized version of
    %delta_k_matrix. Now we replace delta_k_matrix by this antisymmetrised
    %matrix.
    delta_k_matrix=delta_k_matrix_antisymmetrized;

    



    delta_k_updated_matrix=zeros(N_x,N_y);
    e_k=e_k+sum_V_k/((N_x*N_y)^0.5);
    R_k=(((e_k.^2)+(abs(delta_k_matrix)).^2).^0.5);
    K3=(delta_k_matrix./R_k);
    delta_k_updated_matrix=(1/((N_x*N_y)^0.5))*V_k_l*conj(K3);
if(uu==20 || uu==190 || uu==111)
    %figure
    temp=reshape(abs(delta_k_updated_matrix),N_x,N_y);
    abs_delta_k_updated_matrix=temp;
    temp1=circshift( temp, N_y-floor( N_y/2+0.1 ), 2 );
    temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
   % surf((k_y_axis),(k_x_axis),temp2);



%     surf(k_y_axis,k_x_axis,temp);
    %xlabel('k_x');
    %ylabel('k_y');
   % figure
    temp=reshape(atan(imag(delta_k_updated_matrix)./real(delta_k_updated_matrix)),N_x,N_y);
    phase_delta_k_updated_matrix=temp;
    temp1=circshift( temp, N_y-floor( N_y/2+0.1 ), 2 );
    temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
    %surf((k_y_axis),(k_x_axis),temp2);



%     surf(k_y_axis,k_x_axis,reshape(atan(imag(delta_k_updated_matrix)./real(delta_k_updated_matrix)),N_x,N_y));
    %xlabel('k_x');
    %ylabel('k_y');
end

    K1=(1-e_k./R_k);
    K2=(1+e_k./R_k);
    E=(0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3 ));
    MFE_iterations(uu)=E;

    delta_k_matrix=delta_k_updated_matrix;

end

MFE_iterations_dimensionless=MFE_iterations/(8*t);
no_of_iterations_x_axis=1:(n_iterations-1);
%plot(no_of_iterations_x_axis,real(MFE_iterations_dimensionless(2:n_iterations)));
scatter(no_of_iterations_x_axis,real(MFE_iterations_dimensionless(2:n_iterations)));

xlabel('no of iterations');
ylabel('Mean-field energy/Bandwidth');
hold on
% constant_line_global_min=ones((n_iterations-1),1);

%----------------------------------------------------------------
%starting close to delta_x=delta_y




a=1;
t=0.1;
mu=-4*t-(0.02);
N_x=80;
N_y=80;
mod_g_real_space=14;                                                      %interaction strength
% pi=3.14;
% k_0=1/pi;
delta_x_max=0.3;
delta_y_max=0.3;
delta_x_min=0;
delta_y_min=0;
n_delta=40;
step_delta_x=(delta_x_max-(delta_x_min))/n_delta;
step_delta_y=(delta_y_max-(delta_y_min))/n_delta;
matrix_E=zeros(n_delta,n_delta);
gamma=1;
digits(20);
N_particles_array=zeros(n_delta,n_delta);
n_iterations=80;
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
%surf(V_p_r*mod_g_real_space);
%xlabel('Delta_x');
%ylabel('Delta_y');
%zlabel('potential');






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
    E=0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3);
    matrix_E(mm,nn)=E;
   



    end

end
toc
delta_x_axis=delta_x_min:step_delta_x:(delta_x_max-step_delta_x);
delta_y_axis=delta_y_min:step_delta_y:(delta_y_max-step_delta_y);
%figure
%surf(delta_x_axis,delta_y_axis,real(matrix_E));
%title(['free electron band between -4 to +4  and chemical potential at ',num2str(mu)]);
%xlabel('delta_x');
%ylabel('delta_y');
%zlabel('mean-field energy');
%figure
%surf(delta_x_axis,delta_y_axis,N_particles_array);







sorted_matrix_E=sort(sort(real(matrix_E)),2);
min_E=sorted_matrix_E(1,1);
MFE_iterations(1)=min_E;
[min_x,min_y]=find(real(matrix_E)==min_E);
min_x=delta_x_min+(min_x-1)*step_delta_x
min_y=delta_y_min+(min_y-1)*step_delta_y

k_x=(0:(2*pi/N_x):(2*pi/N_x)*(N_x-1))';
k_y=(0:(2*pi/N_y):(2*pi/N_y)*(N_y-1))';
k_array=[ repelem(k_x,numel(k_y)) repmat(k_y,numel(k_x),1) ];
delta_k_matrix=min_x*sin(k_array(:,1))-1i*(min_y+0.05)*sin(k_array(:,2));
k_x_axis=0:(2*pi/(N_x)):(2*pi-(2*pi/(N_x)));
k_y_axis=0:(2*pi/(N_y)):(2*pi-(2*pi/(N_y)));

% figure
% abs_delta_k_matrix=reshape(abs(delta_k_matrix),N_x,N_y);
% temp1=circshift( abs_delta_k_matrix, N_y-floor( N_y/2+0.1 ), 2 );
% temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
% surf((k_y_axis),(k_x_axis),temp2);

delta_k_matrix_shifted=reshape(delta_k_matrix,N_x,N_y);
delta_k_matrix_shifted=circshift( delta_k_matrix_shifted, N_y-floor( N_y/2+0.1 ), 2 );
delta_k_matrix_shifted=circshift( delta_k_matrix_shifted, N_x-floor( N_x/2+0.1 ), 1 );
%figure
%surf(k_y_axis, k_x_axis, abs(delta_k_matrix_shifted))
%xlabel('k_x');
%ylabel('k_y');
%zlabel('mod(Delta_k)');
%figure
%surf(k_y_axis, k_x_axis, angle(delta_k_matrix_shifted))
%xlabel('k_x');
%ylabel('k_y');
%zlabel('phase(Delta_k)');
% figure
% surf(k_y_axis, k_x_axis, abs(delta_k_matrix_shifted-rot90(delta_k_matrix_shifted)))
% 




phase_delta_k=reshape(atan(imag(delta_k_matrix)./real(delta_k_matrix)),N_x,N_y);
temp1=circshift( phase_delta_k, N_y-floor( N_y/2+0.1 ), 2 );
temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );

%figure
%surf(k_y_axis,k_x_axis,reshape(phase_delta_k,N_x,N_y));



for uu=2:n_iterations


    %antisymmetrization of delta_k
    delta_k_matrix_reshaped_minusky=zeros(N_x,N_y);
    delta_k_matrix_reshaped_minusky_minuskx=zeros(N_x,N_y);


    delta_k_matrix_reshaped=reshape(delta_k_matrix,N_x,N_y);
    delta_k_matrix_reshaped_minusky(:,1)=delta_k_matrix_reshaped(:,1);
    delta_k_matrix_reshaped_minusky(:,2:N_y)=flip(delta_k_matrix_reshaped(:,2:N_y),2);
    delta_k_matrix_reshaped_minusky_minuskx(1,:)=delta_k_matrix_reshaped_minusky(1,:);
    delta_k_matrix_reshaped_minusky_minuskx(2:N_x,:)=flip(delta_k_matrix_reshaped_minusky(2:N_x,:),1);
    delta_k_matrix_reshaped_antisymmetrized=0.5*(delta_k_matrix_reshaped-delta_k_matrix_reshaped_minusky_minuskx);
    delta_k_matrix_antisymmetrized=reshape(delta_k_matrix_reshaped_antisymmetrized,N_x*N_y,1);
    %delta_k_matrix_anmtisymmetrised is the antisymmetrized version of
    %delta_k_matrix. Now we replace delta_k_matrix by this antisymmetrised
    %matrix.
    delta_k_matrix=delta_k_matrix_antisymmetrized;

    



    delta_k_updated_matrix=zeros(N_x,N_y);
    e_k=e_k+sum_V_k/((N_x*N_y)^0.5);
    R_k=(((e_k.^2)+(abs(delta_k_matrix)).^2).^0.5);
    K3=(delta_k_matrix./R_k);
    delta_k_updated_matrix=(1/((N_x*N_y)^0.5))*V_k_l*conj(K3);
if(uu==20 || uu==190 || uu==111)
    %figure
    temp=reshape(abs(delta_k_updated_matrix),N_x,N_y);
    abs_delta_k_updated_matrix=temp;
    temp1=circshift( temp, N_y-floor( N_y/2+0.1 ), 2 );
    temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
   % surf((k_y_axis),(k_x_axis),temp2);



%     surf(k_y_axis,k_x_axis,temp);
    %xlabel('k_x');
    %ylabel('k_y');
   % figure
    temp=reshape(atan(imag(delta_k_updated_matrix)./real(delta_k_updated_matrix)),N_x,N_y);
    phase_delta_k_updated_matrix=temp;
    temp1=circshift( temp, N_y-floor( N_y/2+0.1 ), 2 );
    temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
    %surf((k_y_axis),(k_x_axis),temp2);



%     surf(k_y_axis,k_x_axis,reshape(atan(imag(delta_k_updated_matrix)./real(delta_k_updated_matrix)),N_x,N_y));
    %xlabel('k_x');
    %ylabel('k_y');
end

    K1=(1-e_k./R_k);
    K2=(1+e_k./R_k);
    E=(0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3 ));
    MFE_iterations(uu)=E;

    delta_k_matrix=delta_k_updated_matrix;

end

MFE_iterations_dimensionless=MFE_iterations/(8*t);
no_of_iterations_x_axis=1:(n_iterations-1);
%plot(no_of_iterations_x_axis,real(MFE_iterations_dimensionless(2:n_iterations)));
scatter(no_of_iterations_x_axis,real(MFE_iterations_dimensionless(2:n_iterations)));

xlabel('no of iterations');
ylabel('Mean-field energy/Bandwidth');
hold on




%---------------------------------------------------------------------------
%starting on delta_y=0





a=1;
t=0.1;
mu=-4*t-(0.02);
N_x=80;
N_y=80;
mod_g_real_space=14;                                                      %interaction strength
% pi=3.14;
% k_0=1/pi;
delta_x_max=0.3;
delta_y_max=0.3;
delta_x_min=0;
delta_y_min=0;
n_delta=40;
step_delta_x=(delta_x_max-(delta_x_min))/n_delta;
step_delta_y=(delta_y_max-(delta_y_min))/n_delta;
matrix_E=zeros(n_delta,n_delta);
gamma=1;
digits(20);
N_particles_array=zeros(n_delta,n_delta);
n_iterations=80;
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
%surf(V_p_r*mod_g_real_space);
%xlabel('Delta_x');
%ylabel('Delta_y');
%zlabel('potential');






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
    E=0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3);
    matrix_E(mm,nn)=E;
   



    end

end
toc
delta_x_axis=delta_x_min:step_delta_x:(delta_x_max-step_delta_x);
delta_y_axis=delta_y_min:step_delta_y:(delta_y_max-step_delta_y);
%figure
%surf(delta_x_axis,delta_y_axis,real(matrix_E));
%title(['free electron band between -4 to +4  and chemical potential at ',num2str(mu)]);
%xlabel('delta_x');
%ylabel('delta_y');
%zlabel('mean-field energy');
%figure
%surf(delta_x_axis,delta_y_axis,N_particles_array);







sorted_matrix_E=sort(sort(real(matrix_E)),2);
min_E=sorted_matrix_E(1,1);
MFE_iterations(1)=min_E;
[min_x,min_y]=find(real(matrix_E)==min_E);
min_x=delta_x_min+(min_x-1)*step_delta_x
min_y=delta_y_min+(min_y-1)*step_delta_y

k_x=(0:(2*pi/N_x):(2*pi/N_x)*(N_x-1))';
k_y=(0:(2*pi/N_y):(2*pi/N_y)*(N_y-1))';
k_array=[ repelem(k_x,numel(k_y)) repmat(k_y,numel(k_x),1) ];
delta_k_matrix=min_x*sin(k_array(:,1))-1i*(0)*sin(k_array(:,2));
k_x_axis=0:(2*pi/(N_x)):(2*pi-(2*pi/(N_x)));
k_y_axis=0:(2*pi/(N_y)):(2*pi-(2*pi/(N_y)));

% figure
% abs_delta_k_matrix=reshape(abs(delta_k_matrix),N_x,N_y);
% temp1=circshift( abs_delta_k_matrix, N_y-floor( N_y/2+0.1 ), 2 );
% temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
% surf((k_y_axis),(k_x_axis),temp2);

delta_k_matrix_shifted=reshape(delta_k_matrix,N_x,N_y);
delta_k_matrix_shifted=circshift( delta_k_matrix_shifted, N_y-floor( N_y/2+0.1 ), 2 );
delta_k_matrix_shifted=circshift( delta_k_matrix_shifted, N_x-floor( N_x/2+0.1 ), 1 );
%figure
%surf(k_y_axis, k_x_axis, abs(delta_k_matrix_shifted))
%xlabel('k_x');
%ylabel('k_y');
%zlabel('mod(Delta_k)');
%figure
%surf(k_y_axis, k_x_axis, angle(delta_k_matrix_shifted))
%xlabel('k_x');
%ylabel('k_y');
%zlabel('phase(Delta_k)');
% figure
% surf(k_y_axis, k_x_axis, abs(delta_k_matrix_shifted-rot90(delta_k_matrix_shifted)))
% 




phase_delta_k=reshape(atan(imag(delta_k_matrix)./real(delta_k_matrix)),N_x,N_y);
temp1=circshift( phase_delta_k, N_y-floor( N_y/2+0.1 ), 2 );
temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );

%figure
%surf(k_y_axis,k_x_axis,reshape(phase_delta_k,N_x,N_y));



for uu=2:n_iterations


    %antisymmetrization of delta_k
    delta_k_matrix_reshaped_minusky=zeros(N_x,N_y);
    delta_k_matrix_reshaped_minusky_minuskx=zeros(N_x,N_y);


    delta_k_matrix_reshaped=reshape(delta_k_matrix,N_x,N_y);
    delta_k_matrix_reshaped_minusky(:,1)=delta_k_matrix_reshaped(:,1);
    delta_k_matrix_reshaped_minusky(:,2:N_y)=flip(delta_k_matrix_reshaped(:,2:N_y),2);
    delta_k_matrix_reshaped_minusky_minuskx(1,:)=delta_k_matrix_reshaped_minusky(1,:);
    delta_k_matrix_reshaped_minusky_minuskx(2:N_x,:)=flip(delta_k_matrix_reshaped_minusky(2:N_x,:),1);
    delta_k_matrix_reshaped_antisymmetrized=0.5*(delta_k_matrix_reshaped-delta_k_matrix_reshaped_minusky_minuskx);
    delta_k_matrix_antisymmetrized=reshape(delta_k_matrix_reshaped_antisymmetrized,N_x*N_y,1);
    %delta_k_matrix_anmtisymmetrised is the antisymmetrized version of
    %delta_k_matrix. Now we replace delta_k_matrix by this antisymmetrised
    %matrix.
    delta_k_matrix=delta_k_matrix_antisymmetrized;

    



    delta_k_updated_matrix=zeros(N_x,N_y);
    e_k=e_k+sum_V_k/((N_x*N_y)^0.5);
    R_k=(((e_k.^2)+(abs(delta_k_matrix)).^2).^0.5);
    K3=(delta_k_matrix./R_k);
    delta_k_updated_matrix=(1/((N_x*N_y)^0.5))*V_k_l*conj(K3);
if(uu==20 || uu==190 || uu==111)
    %figure
    temp=reshape(abs(delta_k_updated_matrix),N_x,N_y);
    abs_delta_k_updated_matrix=temp;
    temp1=circshift( temp, N_y-floor( N_y/2+0.1 ), 2 );
    temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
   % surf((k_y_axis),(k_x_axis),temp2);



%     surf(k_y_axis,k_x_axis,temp);
    %xlabel('k_x');
    %ylabel('k_y');
   % figure
    temp=reshape(atan(imag(delta_k_updated_matrix)./real(delta_k_updated_matrix)),N_x,N_y);
    phase_delta_k_updated_matrix=temp;
    temp1=circshift( temp, N_y-floor( N_y/2+0.1 ), 2 );
    temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
    %surf((k_y_axis),(k_x_axis),temp2);



%     surf(k_y_axis,k_x_axis,reshape(atan(imag(delta_k_updated_matrix)./real(delta_k_updated_matrix)),N_x,N_y));
    %xlabel('k_x');
    %ylabel('k_y');
end

    K1=(1-e_k./R_k);
    K2=(1+e_k./R_k);
    E=(0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3 ));
    MFE_iterations(uu)=E;

    delta_k_matrix=delta_k_updated_matrix;

end

MFE_iterations_dimensionless=MFE_iterations/(8*t);
no_of_iterations_x_axis=1:(n_iterations-1);
%plot(no_of_iterations_x_axis,real(MFE_iterations_dimensionless(2:n_iterations)));
scatter(no_of_iterations_x_axis,real(MFE_iterations_dimensionless(2:n_iterations)));

xlabel('no of iterations');
ylabel('Mean-field energy/Bandwidth');
hold on


%---------------------------------------------------------------------------
%starting close to delta_y=0




a=1;
t=0.1;
mu=-4*t-(0.02);
N_x=80;
N_y=80;
mod_g_real_space=14;                                                      %interaction strength
% pi=3.14;
% k_0=1/pi;
delta_x_max=0.3;
delta_y_max=0.3;
delta_x_min=0;
delta_y_min=0;
n_delta=40;
step_delta_x=(delta_x_max-(delta_x_min))/n_delta;
step_delta_y=(delta_y_max-(delta_y_min))/n_delta;
matrix_E=zeros(n_delta,n_delta);
gamma=1;
digits(20);
N_particles_array=zeros(n_delta,n_delta);
n_iterations=80;
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
%surf(V_p_r*mod_g_real_space);
%xlabel('Delta_x');
%ylabel('Delta_y');
%zlabel('potential');






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
    E=0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3);
    matrix_E(mm,nn)=E;
   



    end

end
toc
delta_x_axis=delta_x_min:step_delta_x:(delta_x_max-step_delta_x);
delta_y_axis=delta_y_min:step_delta_y:(delta_y_max-step_delta_y);
%figure
%surf(delta_x_axis,delta_y_axis,real(matrix_E));
%title(['free electron band between -4 to +4  and chemical potential at ',num2str(mu)]);
%xlabel('delta_x');
%ylabel('delta_y');
%zlabel('mean-field energy');
%figure
%surf(delta_x_axis,delta_y_axis,N_particles_array);







sorted_matrix_E=sort(sort(real(matrix_E)),2);
min_E=sorted_matrix_E(1,1);
MFE_iterations(1)=min_E;
[min_x,min_y]=find(real(matrix_E)==min_E);
min_x=delta_x_min+(min_x-1)*step_delta_x
min_y=delta_y_min+(min_y-1)*step_delta_y

k_x=(0:(2*pi/N_x):(2*pi/N_x)*(N_x-1))';
k_y=(0:(2*pi/N_y):(2*pi/N_y)*(N_y-1))';
k_array=[ repelem(k_x,numel(k_y)) repmat(k_y,numel(k_x),1) ];
delta_k_matrix=min_x*sin(k_array(:,1))-1i*(0+0.05)*sin(k_array(:,2));
k_x_axis=0:(2*pi/(N_x)):(2*pi-(2*pi/(N_x)));
k_y_axis=0:(2*pi/(N_y)):(2*pi-(2*pi/(N_y)));

% figure
% abs_delta_k_matrix=reshape(abs(delta_k_matrix),N_x,N_y);
% temp1=circshift( abs_delta_k_matrix, N_y-floor( N_y/2+0.1 ), 2 );
% temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
% surf((k_y_axis),(k_x_axis),temp2);

delta_k_matrix_shifted=reshape(delta_k_matrix,N_x,N_y);
delta_k_matrix_shifted=circshift( delta_k_matrix_shifted, N_y-floor( N_y/2+0.1 ), 2 );
delta_k_matrix_shifted=circshift( delta_k_matrix_shifted, N_x-floor( N_x/2+0.1 ), 1 );
%figure
%surf(k_y_axis, k_x_axis, abs(delta_k_matrix_shifted))
%xlabel('k_x');
%ylabel('k_y');
%zlabel('mod(Delta_k)');
%figure
%surf(k_y_axis, k_x_axis, angle(delta_k_matrix_shifted))
%xlabel('k_x');
%ylabel('k_y');
%zlabel('phase(Delta_k)');
% figure
% surf(k_y_axis, k_x_axis, abs(delta_k_matrix_shifted-rot90(delta_k_matrix_shifted)))
% 




phase_delta_k=reshape(atan(imag(delta_k_matrix)./real(delta_k_matrix)),N_x,N_y);
temp1=circshift( phase_delta_k, N_y-floor( N_y/2+0.1 ), 2 );
temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );

%figure
%surf(k_y_axis,k_x_axis,reshape(phase_delta_k,N_x,N_y));



for uu=2:n_iterations


    %antisymmetrization of delta_k
    delta_k_matrix_reshaped_minusky=zeros(N_x,N_y);
    delta_k_matrix_reshaped_minusky_minuskx=zeros(N_x,N_y);


    delta_k_matrix_reshaped=reshape(delta_k_matrix,N_x,N_y);
    delta_k_matrix_reshaped_minusky(:,1)=delta_k_matrix_reshaped(:,1);
    delta_k_matrix_reshaped_minusky(:,2:N_y)=flip(delta_k_matrix_reshaped(:,2:N_y),2);
    delta_k_matrix_reshaped_minusky_minuskx(1,:)=delta_k_matrix_reshaped_minusky(1,:);
    delta_k_matrix_reshaped_minusky_minuskx(2:N_x,:)=flip(delta_k_matrix_reshaped_minusky(2:N_x,:),1);
    delta_k_matrix_reshaped_antisymmetrized=0.5*(delta_k_matrix_reshaped-delta_k_matrix_reshaped_minusky_minuskx);
    delta_k_matrix_antisymmetrized=reshape(delta_k_matrix_reshaped_antisymmetrized,N_x*N_y,1);
    %delta_k_matrix_anmtisymmetrised is the antisymmetrized version of
    %delta_k_matrix. Now we replace delta_k_matrix by this antisymmetrised
    %matrix.
    delta_k_matrix=delta_k_matrix_antisymmetrized;

    



    delta_k_updated_matrix=zeros(N_x,N_y);
    e_k=e_k+sum_V_k/((N_x*N_y)^0.5);
    R_k=(((e_k.^2)+(abs(delta_k_matrix)).^2).^0.5);
    K3=(delta_k_matrix./R_k);
    delta_k_updated_matrix=(1/((N_x*N_y)^0.5))*V_k_l*conj(K3);
if(uu==20 || uu==190 || uu==111)
    %figure
    temp=reshape(abs(delta_k_updated_matrix),N_x,N_y);
    abs_delta_k_updated_matrix=temp;
    temp1=circshift( temp, N_y-floor( N_y/2+0.1 ), 2 );
    temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
   % surf((k_y_axis),(k_x_axis),temp2);



%     surf(k_y_axis,k_x_axis,temp);
    %xlabel('k_x');
    %ylabel('k_y');
   % figure
    temp=reshape(atan(imag(delta_k_updated_matrix)./real(delta_k_updated_matrix)),N_x,N_y);
    phase_delta_k_updated_matrix=temp;
    temp1=circshift( temp, N_y-floor( N_y/2+0.1 ), 2 );
    temp2=circshift( temp1, N_x-floor( N_x/2+0.1 ), 1 );
    %surf((k_y_axis),(k_x_axis),temp2);



%     surf(k_y_axis,k_x_axis,reshape(atan(imag(delta_k_updated_matrix)./real(delta_k_updated_matrix)),N_x,N_y));
    %xlabel('k_x');
    %ylabel('k_y');
end

    K1=(1-e_k./R_k);
    K2=(1+e_k./R_k);
    E=(0.5*((e_k)')*K1+0.25*(1/((N_x*N_y)^0.5))*( ((K3)')*V_k_l*K3 ));
    MFE_iterations(uu)=E;

    delta_k_matrix=delta_k_updated_matrix;

end

MFE_iterations_dimensionless=MFE_iterations/(8*t);
no_of_iterations_x_axis=1:(n_iterations-1);
scatter(no_of_iterations_x_axis,real(MFE_iterations_dimensionless(2:n_iterations)));
hold on

constant_line_global_min=ones((n_iterations-1),1)*MFE_iterations_dimensionless(1);
plot(no_of_iterations_x_axis,constant_line_global_min,'LineStyle','--')

xlabel('no of iterations');
ylabel('Mean-field energy/Bandwidth');

legend('starting at global optimum on \Delta_x=\Delta_y','starting close to global optimum','starting on \Delta_y=0 axis','starting close to \Delta_y=0 axis','variational global optimum')
txt={'Lattice: 80 \times 80'}

text(4,-10.5,txt)
txt={'g_{int}/Bandwidth=0.53'}

text(4,-11.3,txt)
txt={'E_{gap}/Bandwidth=0.025'}

text(4,-11.9,txt)


































