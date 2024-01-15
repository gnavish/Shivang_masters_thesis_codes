

%critical int strength using variational code


a=1;
t=0.1;

mu_min=-4*t-0.80;
mu_max=-4*t-0.02;
n_mu=20;
step_mu=(mu_max-mu_min)/n_mu;


%mu=-4*t-(0.05);
%mu=-4*t-0.36;
N_x=13;
N_y=13;
%mod_g_real_space=14.1;                                                      %interaction strength


mod_g_real_space_min=11;
mod_g_real_space_max=45;



% pi=3.14;
% k_0=1/pi;
delta_x_max=0.03;
delta_y_max=0.03;
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
tol1=0.001;                   %tolerance for accuracy of critical int strength
critical_int_strength=0;
critical_int_strength_nearest_neighbour=0;
critical_int_strength_diff_mu=zeros(n_mu,1);
critical_int_strength_diff_mu_nearest_neighbour=zeros(n_mu,1);

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
% surf(V_p_r*mod_g_real_space);
% xlabel('Delta_x');
% ylabel('Delta_y');
% zlabel('potential');






%V_k_space_1d=V_k_space_1d*(mod_g_real_space)/((N_x*N_y)^0.5);
sum_V_k=0;


sum_V_k=sum(sum(V_k_space_1d));


















%tic
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
      figure
      surf(reshape(V_r,N_x,N_y))
    ft_matrix=zeros(N_x*N_y,N_x*N_y);
    ft_matrix=exp( 1i*k_array(:,1)*(R_array(1,:)-1)+1i*k_array(:,2)*(R_array(2,:)-1) );
    V_k_1d=ft_matrix*V_r;
    V_k_1d_t=(V_k_1d)';
    sum_V_k=sum(V_k_1d);

    %e_k=-2*t*( cos(k_array(:,1)) +cos(k_array(:,2)) )-mu;
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

%toc
V_k_l_original=V_k_l;   %storing result of fourier transform without rescaling by int strength

for pp=1:n_mu

    mu=mu_min+(pp-1)*step_mu;
    e_k=-2*t*( cos(k_array(:,1)) +cos(k_array(:,2)) )-mu;

    mod_g_real_space_min=11;
    mod_g_real_space_max=45;



for uu=1:1000

mod_g_real_space_mid=(mod_g_real_space_max+mod_g_real_space_min)/2;
mod_g_real_space=mod_g_real_space_mid;




V_k_l=V_k_l_original*(mod_g_real_space)/((N_x*N_y)^0.5);
%  V_k_l1=V_k_l1*(mod_g_real_space)/((N_x*N_y)^0.5);
%  V_k_l=V_k_l1;
% plot(real(V_k_l1(:,1)))





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




%toc
% delta_x_axis=delta_x_min:step_delta_x:(delta_x_max-step_delta_x);
% delta_y_axis=delta_y_min:step_delta_y:(delta_y_max-step_delta_y);
% figure
% surf(delta_x_axis,delta_y_axis,real(matrix_E));
% title(['free electron band between -4 to +4  and chemical potential at ',num2str(mu)]);
% xlabel('delta_x');
% ylabel('delta_y');
% zlabel('mean-field energy');
% figure
% surf(delta_x_axis,delta_y_axis,N_particles_array);







sorted_matrix_E=sort(sort(real(matrix_E)),2);
min_E=sorted_matrix_E(1,1);
MFE_iterations(1)=min_E;
[min_x,min_y]=find(real(matrix_E)==min_E);
min_x=delta_x_min+(min_x-1)*step_delta_x;
min_y=delta_y_min+(min_y-1)*step_delta_y;


if(abs(mod_g_real_space_max-mod_g_real_space_min)<tol1)

    critical_int_strength=mod_g_real_space_mid;
    critical_int_strength_nearest_neighbour=V_r(2)*critical_int_strength;
    critical_int_strength_diff_mu(pp)=critical_int_strength;
    critical_int_strength_diff_mu_nearest_neighbour(pp)=critical_int_strength_nearest_neighbour;
    break;

end



if(size(min_x)==1)

    if(min_x==0 && min_y==0)

        mod_g_real_space_min=mod_g_real_space_mid;


    else

        mod_g_real_space_max=mod_g_real_space_mid;

    end


end


if(size(min_x,1)==2)

    if(min_x(1)==0 && min_x(2)==0 && min_y(1)==0 && min_y(2)==0)

        mod_g_real_space_min=mod_g_real_space_mid;


    else

        mod_g_real_space_max=mod_g_real_space_mid;


    end



end










end


end






%--------------------------------------------------------------------------------------
%critical int strength using ED



tic




a=1;
t=0.1;
%mu=-4*t-(0.02);

%mu_min=-4*t-0.07;
%mu_max=-4*t-0.02;
%n_mu=3;
%step_mu=(mu_max-mu_min)/n_mu;


%mu=-4*t-0.05;
%mu=-2*t-0.001;
%N_x=13;
%N_y=13;
%mod_g_real_space=11.6;
%mod_g_real_space=13.7;                                                      %interaction strength
% pi=3.14;
% k_0=1/pi;
delta_x_max=0.05;
delta_y_max=0.05;
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
%tol1=0.01;                           %tolerance for critical interaction strength
critical_int_strength=0;
critical_int_strength_nearest_neighbour=0;
critical_int_strength_diff_mu_ED=zeros(n_mu,1);
critical_int_strength_diff_mu_nearest_neighbour_ED=zeros(n_mu,1);

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
% V_p_r=(1/(num_sum_p^0.5))*V_p_r*mod_g_real_space;
V_p_r(1,1)=0;
surf(V_p_r);
xlabel('x');
ylabel('y');
zlabel('potential');
% plot(V_p_r)
% xlabel('x')
% ylabel('potential')




% %defining potential
%  X_array=(1:1:N_x)';
%     Y_array=(1:1:N_y)';
%     R_array=([ repelem(X_array,numel(Y_array)) repmat(Y_array,numel(X_array),1) ])';
%     R_array_t=(R_array)';
%     V_r=zeros(N_x*N_y,1);
%     for ii=-(num_sum_p-1)/2:(num_sum_p-1)/2
%         for jj=-(num_sum_p-1)/2:(num_sum_p-1)/2
% 
%             V_r=V_r+(B-1)*exp( -( (R_array_t(:,1)-1+ii*N_x).^2 + (R_array_t(:,2)-1+jj*N_y).^2 )/sigma12 ) - B*exp( -( (R_array_t(:,1)-1+ii*N_x).^2 + (R_array_t(:,2)-1+jj*N_y).^2)/sigma22 );
% 
%         end
%     end
%       V_r(1,1)=0;
%       V_r=(1/(num_sum_p^0.5))*V_r;
%       figure
%       surf(reshape(V_r,N_x,N_y))
% 
% 
% 
% 


% %potential for 3X1 lattice
% V_p_r(2)=-0.14;
% V_p_r(3)=V_p_r(2);
% V_p_r(4)=V_p_r(2);
% V_p_r(5)=V_p_r(2);
% plot(V_p_r)










%position basis construction
position_basis=zeros(N_x*N_y,2);
for ii=1:N_x

    for jj=1:N_y

        position_basis((ii-1)*N_y+jj,1)=ii;
        position_basis((ii-1)*N_y+jj,2)=jj;

    end

end

%two particle Hilbert space construction
hilbert_basis=zeros((N_x*N_y)*(N_x*N_y-1)/2,4);
temp=0;
for ii=1:(N_x*N_y-1)

    for jj=(ii+1):(N_x*N_y)

        temp=temp+1;
        hilbert_basis(temp,1)=position_basis(ii,1);
        hilbert_basis(temp,2)=position_basis(ii,2);
        hilbert_basis(temp,3)=position_basis(jj,1);
        hilbert_basis(temp,4)=position_basis(jj,2);

    end

end


V_p_r_original=V_p_r;


for pp=1:n_mu
    mu=mu_min+(pp-1)*step_mu;

    mod_g_real_space_min=11;
    mod_g_real_space_max=45;




for uu=1:1000

    mod_g_real_space_mid=(mod_g_real_space_max+mod_g_real_space_min)/2;
    mod_g_real_space=mod_g_real_space_mid;

    V_p_r=(1/(num_sum_p^0.5))*V_p_r_original*mod_g_real_space;








%populating the Hamiltonian matrix
%H=sparse((N_x*N_y)*(N_x*N_y-1)/2,(N_x*N_y)*(N_x*N_y-1)/2);
H=zeros((N_x*N_y)*(N_x*N_y-1)/2,(N_x*N_y)*(N_x*N_y-1)/2);

for ii=1:(N_x*N_y)*(N_x*N_y-1)/2

    for jj=1:(N_x*N_y)*(N_x*N_y-1)/2

        %bra
        r1_x=hilbert_basis(ii,1);
        r1_y=hilbert_basis(ii,2);
        r2_x=hilbert_basis(ii,3);
        r2_y=hilbert_basis(ii,4);
        %ket
        r1_prime_x=hilbert_basis(jj,1);
        r1_prime_y=hilbert_basis(jj,2);
        r2_prime_x=hilbert_basis(jj,3);
        r2_prime_y=hilbert_basis(jj,4);

        % r1_x==r1_prime_x 
        % r1_y==r1_prime_y 
        %  ((mod(r2_x-r2_prime_x,N_x))^2+(mod(r2_y-r2_prime_y,N_y))^2)
        %  ((mod(r2_prime_x-r2_x,N_x))^2+(mod(r2_prime_y-r2_y,N_y))^2)

        if(r1_x==r1_prime_x && r1_y==r1_prime_y && ( ((mod(r2_x-r2_prime_x,N_x))^2+(mod(r2_y-r2_prime_y,N_y))^2)==1 || ((mod(r2_prime_x-r2_x,N_x))^2+(mod(r2_prime_y-r2_y,N_y))^2)==1 ) )

            H(ii,jj)=H(ii,jj)-t;


        end

         if(r1_x==r2_prime_x && r1_y==r2_prime_y && ( ((mod(r2_x-r1_prime_x,N_x))^2+(mod(r2_y-r1_prime_y,N_y))^2)==1 || ((mod(r1_prime_x-r2_x,N_x))^2+(mod(r1_prime_y-r2_y,N_y))^2)==1 ) )

            H(ii,jj)=H(ii,jj)+t;

        end

         if(r2_x==r1_prime_x && r2_y==r1_prime_y && ( ((mod(r1_x-r2_prime_x,N_x))^2+(mod(r1_y-r2_prime_y,N_y))^2)==1 || ((mod(r2_prime_x-r1_x,N_x))^2+(mod(r2_prime_y-r1_y,N_y))^2)==1 ) )

            H(ii,jj)=H(ii,jj)+t;


        end

         if(r2_x==r2_prime_x && r2_y==r2_prime_y && ( ((mod(r1_x-r1_prime_x,N_x))^2+(mod(r1_y-r1_prime_y,N_y))^2)==1 || ((mod(r1_prime_x-r1_x,N_x))^2+(mod(r1_prime_y-r1_y,N_y))^2)==1 ) )

            H(ii,jj)=H(ii,jj)-t;


         end

         if(r1_x==r2_prime_x && r1_y==r2_prime_y && r2_x==r1_prime_x && r2_y==r1_prime_y)

            H(ii,jj)=H(ii,jj)+2*mu;

         end

        if(r1_x==r1_prime_x && r1_y==r1_prime_y && r2_x==r2_prime_x && r2_y==r2_prime_y)

            H(ii,jj)=H(ii,jj)-2*mu;

        end


% mod(r1_prime_x-r2_prime_x,N_x)+1
% mod(r1_prime_y-r2_prime_y,N_y)+1
% mod(r2_prime_x-r1_prime_x,N_x)+1
% mod(r2_prime_y-r1_prime_y,N_y)+1 

        if(r1_x==r2_prime_x && r1_y==r2_prime_y && r2_x==r1_prime_x && r2_y==r1_prime_y)

            H(ii,jj)=H(ii,jj)-V_p_r( mod(r1_prime_x-r2_prime_x,N_x)+1,mod(r1_prime_y-r2_prime_y,N_y)+1 );
            H(ii,jj)=H(ii,jj)-V_p_r( mod(r2_prime_x-r1_prime_x,N_x)+1,mod(r2_prime_y-r1_prime_y,N_y)+1 );


        end

        if(r1_x==r1_prime_x && r1_y==r1_prime_y && r2_x==r2_prime_x && r2_y==r2_prime_y)

            H(ii,jj)=H(ii,jj)+V_p_r( mod(r1_prime_x-r2_prime_x,N_x)+1  , mod(r1_prime_y-r2_prime_y,N_y)+1  );
            H(ii,jj)=H(ii,jj)+V_p_r( mod(r2_prime_x-r1_prime_x,N_x)+1  , mod(r2_prime_y-r1_prime_y,N_y)+1  );


        end




    end

end

%E=eigs(H,50,'smallestabs');
E=sort(eig(H));
ground_state_energy=E(1);

if(abs(mod_g_real_space_max-mod_g_real_space_min)<tol1)

    critical_int_strength=mod_g_real_space_mid;
    critical_int_strength_nearest_neighbour=V_p_r(1,2);
    critical_int_strength_diff_mu_ED(pp)=critical_int_strength;
    critical_int_strength_diff_mu_nearest_neighbour_ED(pp)=critical_int_strength_nearest_neighbour;
    break;

end

if(ground_state_energy>0)

    mod_g_real_space_min=mod_g_real_space_mid;

else

    mod_g_real_space_max=mod_g_real_space_mid;

end






end






end


toc


gap_x_axis=(-4*t-mu_max+step_mu):step_mu:(-4*t-mu_min);
gap_x_axis_dimensionless=gap_x_axis/(8*t);
critical_int_strength_diff_mu_nearest_neighbour_ED_rev=transpose(fliplr(transpose(critical_int_strength_diff_mu_nearest_neighbour_ED)));
critical_int_strength_diff_mu_nearest_neighbour_rev=transpose(fliplr(transpose(critical_int_strength_diff_mu_nearest_neighbour)));
critical_int_strength_diff_mu_nearest_neighbour_rev_dimless=critical_int_strength_diff_mu_nearest_neighbour_rev/(8*t);
critical_int_strength_diff_mu_nearest_neighbour_ED_rev_dimless=critical_int_strength_diff_mu_nearest_neighbour_ED_rev/(8*t);
figure
plot(gap_x_axis_dimensionless,abs(critical_int_strength_diff_mu_nearest_neighbour_rev_dimless))
hold on
plot(gap_x_axis_dimensionless,abs(critical_int_strength_diff_mu_nearest_neighbour_ED_rev_dimless))

xlabel('E_{gap}/Bandwidth')
ylabel('g_{c}/Bandwidth')
legend('variational result','ED result')



%final plot
%if you have saved variables, just copy the below lines in command window
%to get the final figure
figure
scatter(gap_x_axis_dimensionless,abs(critical_int_strength_diff_mu_nearest_neighbour_rev_dimless))
hold on
scatter(gap_x_axis_dimensionless,abs(critical_int_strength_diff_mu_nearest_neighbour_ED_rev_dimless))

xlabel('E_{gap}/Bandwidth')
ylabel('g_{c}/Bandwidth')
legend('variational result','ED result in the two-fermion sector')
hold on

txt={'Lattice dimensions: 13 \times 13'};

t=text(0.1,1,txt);
axes('Position',[0.2,0.7,0.2,0.2])
box on
diff_between_plots=abs(critical_int_strength_diff_mu_nearest_neighbour_rev_dimless-critical_int_strength_diff_mu_nearest_neighbour_ED_rev_dimless);
scatter(gap_x_axis_dimensionless,diff_between_plots);
xlabel('E_{gap}/Bandwidth')
ylabel('Diff b/w g_{c}/Bandwidth')
% txt={'Lattice dimensions:',num2str(N_x),'\times',num2str(N_y)};
% 
% t=text(0.2,0.6,txt);
% 
% 
