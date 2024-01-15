tic
mu_below_bott_band_max=-0.4;
mu_below_bott_band_min=-0.01;
n_steps=20;
g_real_space_to_reach_mu_diff_mu=zeros(n_steps,1);
g_real_space_to_reach_mu_diff_mu_dimensionless=zeros(n_steps,1);
mu_below_bott_band_step=(mu_below_bott_band_max-mu_below_bott_band_min)/n_steps;
N_x=6;
N_y=6;
% mu=-0.00;
h=0.1;
%pi=3.14;
g_real_space_init=-1;
% lowest_eigenvalue=0;
brill_zone_points=zeros(2*N_x*N_y,3);
basis_set=zeros((N_x*N_y)^2,6);
sorted_eigenvalues=zeros((N_x*N_y)^2,1);
tol=10^(-9);
tol1=10^(-4);



for q=1:2 

    for i=1:N_x
    
        for j=1:N_y
            
            k_x=((2*pi)/(N_x))*(i-1);
            k_y=((2*pi)/(N_y))*(j-1);
            brill_zone_points(N_x*N_y*(q-1)+N_y*(i-1)+j,1)=k_x;
            brill_zone_points(N_x*N_y*(q-1)+N_y*(i-1)+j,2)=k_y;
            brill_zone_points(N_x*N_y*(q-1)+N_y*(i-1)+j,3)=q;
    
        end
    
    end

end

x=1;

for i=1:(N_x*N_y)

    for j=(N_x*N_y+i):(2*N_x*N_y)
        
        
        basis_set(x,1)=brill_zone_points(i,1);
        basis_set(x,2)=brill_zone_points(i,2);
        basis_set(x,3)=brill_zone_points(i,3);
        basis_set(x,4)=brill_zone_points(j,1);
        basis_set(x,5)=brill_zone_points(j,2);
        basis_set(x,6)=brill_zone_points(j,3);
        x=x+1;

    end

end

for i=(N_x*N_y+1):(2*N_x*N_y)

    for j=(i-(N_x*N_y)):(N_x*N_y)

        if((i-j)~=(N_x*N_y))

            basis_set(x,1)=brill_zone_points(i,1);
            basis_set(x,2)=brill_zone_points(i,2);
            basis_set(x,3)=brill_zone_points(i,3);
            basis_set(x,4)=brill_zone_points(j,1);
            basis_set(x,5)=brill_zone_points(j,2);
            basis_set(x,6)=brill_zone_points(j,3);
            x=x+1;

        end

    end

end


for nn=1:n_steps

    mu_below_bott_band=mu_below_bott_band_min+mu_below_bott_band_step*nn;







    

g_real_space_prev=0;
g_real_space_curr=g_real_space_init;
g_real_space_next=0;


lowest_eigenvalue_prev=0;
%calculating lowest_eigenvalue_curr


H=sparse((N_x*N_y)^2,(N_x*N_y)^2);
for i=1:((N_x*N_y)^2)

    for j=1:((N_x*N_y)^2)

        k1_x=basis_set(i,1);
        k1_y=basis_set(i,2);
        q=basis_set(i,3);
        k2_x=basis_set(i,4);
        k2_y=basis_set(i,5);
        r=basis_set(i,6);
        k1_prime_x=basis_set(j,1);
        k1_prime_y=basis_set(j,2);
        s=basis_set(j,3);
        k2_prime_x=basis_set(j,4);
        k2_prime_y=basis_set(j,5);
        t=basis_set(j,6);
        e_k1=-2*h*(cos(k1_x)+cos(k1_y)-2);
        e_k2=-2*h*(cos(k2_x)+cos(k2_y)-2);

        if(abs(k1_x-k1_prime_x)<tol && abs(k1_y-k1_prime_y)<tol && abs(k2_x-k2_prime_x)<tol && abs(k2_y-k2_prime_y)<tol && q==s && r==t)
                                                
            H(i,j)=H(i,j)+e_k1+e_k2;

        
        end

        if(abs(k1_x-k2_prime_x)<tol && abs(k1_y-k2_prime_y)<tol && abs(k1_prime_x-k2_x)<tol && abs(k1_prime_y-k2_y)<tol && q==t && r==s)

            H(i,j)=H(i,j)-(e_k1+e_k2);
            

        end
        
        sum_k_x=k2_prime_x+k1_prime_x-k1_x-k2_x;
        sum_k_y=k2_prime_y+k1_prime_y-k1_y-k2_y;

        if(abs(sum_k_x-2*pi)<tol  || abs(sum_k_x)<tol)

            sum_k_x=0;

        end

        if(abs(sum_k_y-2*pi)<tol  || abs(sum_k_y)<tol)

            sum_k_y=0;

        end

        if(abs(mod(sum_k_x,2*pi))<tol && abs(mod(sum_k_y,2*pi))<tol && q==t && r==s)

            H(i,j)=H(i,j)-2*g_real_space_curr/(N_x*N_y);

        end

        if(abs(mod(sum_k_x,2*pi))<tol && abs(mod(sum_k_y,2*pi))<tol && r==t && q==s)

            H(i,j)=H(i,j)+2*g_real_space_curr/(N_x*N_y);
            

        end


   end

end

lowest_eigenvalue_curr=eigs(H,1,-130);







for uu=1:1000

    if(abs(lowest_eigenvalue_curr-2*mu_below_bott_band)<tol1)

        break;

   

    else

        g_real_space_next=(g_real_space_prev+g_real_space_curr)/2;

        %calculating lowest_eigenvalue_next


        H=sparse((N_x*N_y)^2,(N_x*N_y)^2);
        for i=1:((N_x*N_y)^2)
        
            for j=1:((N_x*N_y)^2)
        
                k1_x=basis_set(i,1);
                k1_y=basis_set(i,2);
                q=basis_set(i,3);
                k2_x=basis_set(i,4);
                k2_y=basis_set(i,5);
                r=basis_set(i,6);
                k1_prime_x=basis_set(j,1);
                k1_prime_y=basis_set(j,2);
                s=basis_set(j,3);
                k2_prime_x=basis_set(j,4);
                k2_prime_y=basis_set(j,5);
                t=basis_set(j,6);
                e_k1=-2*h*(cos(k1_x)+cos(k1_y)-2);
                e_k2=-2*h*(cos(k2_x)+cos(k2_y)-2);
        
                if(abs(k1_x-k1_prime_x)<tol && abs(k1_y-k1_prime_y)<tol && abs(k2_x-k2_prime_x)<tol && abs(k2_y-k2_prime_y)<tol && q==s && r==t)
                                                        
                    H(i,j)=H(i,j)+e_k1+e_k2;
        
                
                end
        
                if(abs(k1_x-k2_prime_x)<tol && abs(k1_y-k2_prime_y)<tol && abs(k1_prime_x-k2_x)<tol && abs(k1_prime_y-k2_y)<tol && q==t && r==s)
        
                    H(i,j)=H(i,j)-(e_k1+e_k2);
                    
        
                end
                
                sum_k_x=k2_prime_x+k1_prime_x-k1_x-k2_x;
                sum_k_y=k2_prime_y+k1_prime_y-k1_y-k2_y;
        
                if(abs(sum_k_x-2*pi)<tol  || abs(sum_k_x)<tol)
        
                    sum_k_x=0;
        
                end
        
                if(abs(sum_k_y-2*pi)<tol  || abs(sum_k_y)<tol)
        
                    sum_k_y=0;
        
                end
        
                if(abs(mod(sum_k_x,2*pi))<tol && abs(mod(sum_k_y,2*pi))<tol && q==t && r==s)
        
                    H(i,j)=H(i,j)-2*g_real_space_next/(N_x*N_y);
        
                end
        
                if(abs(mod(sum_k_x,2*pi))<tol && abs(mod(sum_k_y,2*pi))<tol && r==t && q==s)
        
                    H(i,j)=H(i,j)+2*g_real_space_next/(N_x*N_y);
                    
        
                end
        
        
           end
        
        end
        
        
        lowest_eigenvalue_next=eigs(H,1,-130);

        if(lowest_eigenvalue_next>(2*mu_below_bott_band))

%             g_real_space_prev=g_real_space_curr
%             lowest_eigenvalue_prev=lowest_eigenvalue_curr
            g_real_space_prev=g_real_space_next;
            lowest_eigenvalue_prev=lowest_eigenvalue_next;

            

        end
% lowest_eigenvalue_prev
        if(lowest_eigenvalue_next<(2*mu_below_bott_band) )
            
           
            g_real_space_curr=g_real_space_next;
            lowest_eigenvalue_curr=lowest_eigenvalue_next;
            

        end

        








    end

end

g_real_space_to_reach_mu_diff_mu(nn,1)=g_real_space_curr;
%g_real_space_to_reach_mu_diff_mu_dimensionless(nn,1)=g_real_space_to_reach_mu_diff_mu(nn,1)/(mu_below_bott_band);

end


mu_array=(mu_below_bott_band_min+mu_below_bott_band_step):mu_below_bott_band_step:mu_below_bott_band_max;
%mu_array_rev=fliplr(mu_array);
mu_array=abs(mu_array);
g_real_space_to_reach_mu_diff_mu_rev=fliplr(g_real_space_to_reach_mu_diff_mu);
%plot(mu_array,g_real_space_to_reach_mu_diff_mu);
plot(mu_array,abs(g_real_space_to_reach_mu_diff_mu_rev))
title('int strength required for bound state to reach 2*mu v/s mu')
xlabel('mu below bottom of band')
ylabel('int strength required for bound state to reach 2*mu')
figure
%plot(mu_array,g_real_space_to_reach_mu_diff_mu_dimensionless)
mu_array_dimensionless=mu_array/(8*h);
scatter(mu_array_dimensionless,g_real_space_to_reach_mu_diff_mu_dimensionless);

hold on








% mu_below_bott_band_max=-0.08;
% mu_below_bott_band_min=-0.01;
% n_steps=40;
% g_real_space_to_reach_mu_diff_mu=zeros(n_steps,1);
% g_real_space_to_reach_mu_diff_mu_dimensionless=zeros(n_steps,1);
% mu_below_bott_band_step=(mu_below_bott_band_max-mu_below_bott_band_min)/n_steps;
% N_x=6;
% N_y=6;
% % mu=-0.00;
% h=0.1;
% %pi=3.14;
toc

%-----------------------------------------------------------------------------------------
%variational solution
tic
a=1;
%t=0.1;
t=h;
%mu=-4*t-0.02;
% mu_max=-4*t-0.01;
% mu_min=-4*t-0.08;
mu_max=-4*t+mu_below_bott_band_min;
mu_min=-4*t+mu_below_bott_band_max;
%n_mu=40;
n_mu=n_steps;
step_mu=(mu_max-mu_min)/n_mu;


%mu=-0.2;
% N_x=6;
% N_y=6;
%g_real_space=0.26;
g_real_space=0.2979;
g_real_space_min=10^(-3);
g_real_space_max=5;
minimum_delta_at_g_min=0;
minimum_delta_at_g_max=0;
minimum_delta_at_g_mid=0;




%g_real_space=0.05;
g=g_real_space/((N_x*N_y)^0.5);                                                                        %interaction strength
%pi=3.14;
k_0=1/pi;
delta_x_max=0.03;
delta_y_max=0.03;
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
tol=10^(-4);
critical_g_diff_mu=zeros(n_mu,1);
critical_g_diff_mu_dimensionless=zeros(n_mu,1);


for qq=1:n_mu
    mu=mu_min+(step_mu)*(qq-1);
    g_real_space_min=10^(-3);
g_real_space_max=7;
minimum_delta_at_g_min=0;
minimum_delta_at_g_max=0;
minimum_delta_at_g_mid=0;

for pp=1:10000

    g_real_space=g_real_space_min;


parfor mm=1:n_delta
    delta_x=0+(mm-1)*step_delta_x;
    for nn=1:1%n_delta
    delta_y=0+(nn-1)*step_delta_y;
  E=0;
  E_diff_way=0;
sum_1_by_E_l=0;
sum_v_k_2=0;

for i=1:N_x


    for j=1:N_y

        l_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        l_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_l=-2*t*(cos(l_x*a)+cos(l_y*a))-mu;
        R_l=((e_l^2+(abs(delta_x))^2)^(1/2));
        sum_1_by_E_l=(sum_1_by_E_l+(1/R_l));


    end

end
sum_1_by_E_l_diff_deltas(mm,1)=(sum_1_by_E_l);


for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
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
[val,minimum_delta_at_g_min]=min(matrix_E_1_diff_way);
minimum_delta_at_g_min=step_delta_x*(minimum_delta_at_g_min-1);




%for g_real_space_max


    g_real_space=g_real_space_max;


parfor mm=1:n_delta
    delta_x=0+(mm-1)*step_delta_x;
    for nn=1:1%n_delta
    delta_y=0+(nn-1)*step_delta_y;
  E=0;
  E_diff_way=0;
sum_1_by_E_l=0;
sum_v_k_2=0;

for i=1:N_x


    for j=1:N_y

        l_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        l_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_l=-2*t*(cos(l_x*a)+cos(l_y*a))-mu;
        R_l=((e_l^2+(abs(delta_x))^2)^(1/2));
        sum_1_by_E_l=(sum_1_by_E_l+(1/R_l));


    end

end
sum_1_by_E_l_diff_deltas(mm,1)=(sum_1_by_E_l);


for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
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
[val,minimum_delta_at_g_max]=min(matrix_E_1_diff_way);
minimum_delta_at_g_max=step_delta_x*(minimum_delta_at_g_max-1);





%for g_real_space_mid

    g_real_space=(g_real_space_min+g_real_space_max)/2;
    g_real_space_mid=(g_real_space_min+g_real_space_max)/2;


parfor mm=1:n_delta
    delta_x=0+(mm-1)*step_delta_x;
    for nn=1:1%n_delta
    delta_y=0+(nn-1)*step_delta_y;
  E=0;
  E_diff_way=0;
sum_1_by_E_l=0;
sum_v_k_2=0;

for i=1:N_x


    for j=1:N_y

        l_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        l_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_l=-2*t*(cos(l_x*a)+cos(l_y*a))-mu;
        R_l=((e_l^2+(abs(delta_x))^2)^(1/2));
        sum_1_by_E_l=(sum_1_by_E_l+(1/R_l));


    end

end
sum_1_by_E_l_diff_deltas(mm,1)=(sum_1_by_E_l);


for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
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
[val,minimum_delta_at_g_mid]=min(matrix_E_1_diff_way);
minimum_delta_at_g_mid=step_delta_x*(minimum_delta_at_g_mid-1);



if(abs(g_real_space_max-g_real_space_min)<tol)

    break;

end

if(minimum_delta_at_g_mid<=minimum_delta_at_g_max && minimum_delta_at_g_mid~=minimum_delta_at_g_min)

    g_real_space_max=g_real_space_mid;

end

if(minimum_delta_at_g_mid==minimum_delta_at_g_min)

    g_real_space_min=g_real_space_mid;

end







end
critical_g_diff_mu(qq,1)=(g_real_space_max+g_real_space_min)/2;
critical_g_diff_mu_dimensionless(qq,1)=critical_g_diff_mu(qq,1)/(-4*t-mu);


end

gap_x_axis=(-4*t-(mu_max-step_mu)):step_mu:(-4*t-mu_min);
gap_x_axis_dimensionless=gap_x_axis/(8*t);
%critical_g_diff_mu_dimensionless_rev=fliplr(critical_g_diff_mu_dimensionless);
%figure
critical_g_diff_mu_dimensionless_rev=transpose(fliplr(transpose(critical_g_diff_mu_dimensionless)));
plot(gap_x_axis_dimensionless,critical_g_diff_mu_dimensionless_rev);
ylabel('g_{c}/E_{gap}')
xlabel('E_{gap}/Bandwidth')
legend('Exact diagonalization prediction','Variational solution')
%print -deps critical_int_strength_vs_gap_spinful_fermions_for_thesis

















delta_x_axis=(0+step_delta_x):step_delta_x:delta_x_max;
% figure
% plot(delta_x_axis,matrix_E_1)
% title(['free electron band between -0.4 to +0.4  and chemical potential at ',num2str(mu)]);
% xlabel('delta');
% ylabel('Mean field energy functional value');
% [minima,index_minima]=min(matrix_E_1);
% delta_at_minima=index_minima*step_delta_x;
% figure
% surf(matrix_E_k)
% plot(delta_x_axis,sum_1_by_E_l_diff_deltas)
% figure
% plot(delta_x_axis,sum_nonint_piece);
% figure
% plot(delta_x_axis,sum_int_piece);
% figure
% surf(v_k_one_delta);
% figure
% surf(v_k_one_delta_diff_way);
% diff_v_k_calc_diff_ways=v_k_one_delta_diff_way-v_k_one_delta;
% figure
% surf(diff_v_k_calc_diff_ways)
% figure
% surf(e_k_times_v_k_2_one_delta)
% figure
% plot(delta_x_axis,matrix_E_1_diff_way);
% sum_e_k

% extrema=find(imregionalmax(matrix_E));


%figure
min_gap=-4*t-mu_max;
max_gap=-4*t-mu_min;
gap_x_axis=(min_gap+step_mu):step_mu:max_gap;
gap_x_axis_dimensionless=gap_x_axis/(8*t);


%g_real_space_to_reach_mu_diff_mu_dimless and critical_g_diff_mu_dimless
%are the two relevant variables and are reversed in order. We reverse on of
%them below for plotting.


g_real_space_to_reach_mu_diff_mu_dimless=g_real_space_to_reach_mu_diff_mu/(8*t);
critical_g_diff_mu_dimless=critical_g_diff_mu/(8*t);

% figure
% %plot(gap_x_axis,abs(g_real_space_to_reach_mu_diff_mu))
% %hold on
% scatter(gap_x_axis,transpose(fliplr(transpose(critical_g_diff_mu))))
% 
% %figure
% %plot(gap_x_axis,transpose(fliplr(transpose(critical_g_diff_mu_dimless))));
% hold on
% scatter(gap_x_axis,abs(g_real_space_to_reach_mu_diff_mu_dimless));

%final figure
figure
scatter(gap_x_axis_dimensionless,transpose(fliplr(transpose(critical_g_diff_mu_dimless))))
hold on
scatter(gap_x_axis_dimensionless,abs(g_real_space_to_reach_mu_diff_mu_dimless),'x');

%scatter(gap_x_axis_dimensionless,transpose(fliplr(transpose(critical_g_diff_mu_dimless))),'d')
% xlabel('E_{gap}/Bandwidth');
% ylabel('g_{c}/Bandwidth')
% legend('Variational result','ED result')
hold on

toc



















%----------------------------------------------------------------------------
%variational solution in thermodynamic limit


tic
a=1;
%t=0.1;
t=h;
%mu=-4*t-0.02;
% mu_max=-4*t-0.01;
% mu_min=-4*t-0.08;
mu_max=-4*t+mu_below_bott_band_min;
mu_min=-4*t+mu_below_bott_band_max;
%n_mu=40;
n_mu=n_steps;
step_mu=(mu_max-mu_min)/n_mu;


%mu=-0.2;
 N_x=200;
 N_y=200;
%g_real_space=0.26;
g_real_space=0.2979;
g_real_space_min=10^(-3);
g_real_space_max=5;
minimum_delta_at_g_min=0;
minimum_delta_at_g_max=0;
minimum_delta_at_g_mid=0;




%g_real_space=0.05;
g=g_real_space/((N_x*N_y)^0.5);                                                                        %interaction strength
%pi=3.14;
k_0=1/pi;
delta_x_max=0.03;
delta_y_max=0.03;
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
tol=10^(-4);
critical_g_diff_mu_th_limit=zeros(n_mu,1);
critical_g_diff_mu_dimensionless=zeros(n_mu,1);


for qq=1:n_mu
    mu=mu_min+(step_mu)*(qq-1);
    g_real_space_min=10^(-3);
g_real_space_max=7;
minimum_delta_at_g_min=0;
minimum_delta_at_g_max=0;
minimum_delta_at_g_mid=0;

for pp=1:10000

    g_real_space=g_real_space_min;


parfor mm=1:n_delta
    delta_x=0+(mm-1)*step_delta_x;
    for nn=1:1%n_delta
    delta_y=0+(nn-1)*step_delta_y;
  E=0;
  E_diff_way=0;
sum_1_by_E_l=0;
sum_v_k_2=0;

for i=1:N_x


    for j=1:N_y

        l_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        l_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_l=-2*t*(cos(l_x*a)+cos(l_y*a))-mu;
        R_l=((e_l^2+(abs(delta_x))^2)^(1/2));
        sum_1_by_E_l=(sum_1_by_E_l+(1/R_l));


    end

end
sum_1_by_E_l_diff_deltas(mm,1)=(sum_1_by_E_l);


for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
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
[val,minimum_delta_at_g_min]=min(matrix_E_1_diff_way);
minimum_delta_at_g_min=step_delta_x*(minimum_delta_at_g_min-1);




%for g_real_space_max


    g_real_space=g_real_space_max;


parfor mm=1:n_delta
    delta_x=0+(mm-1)*step_delta_x;
    for nn=1:1%n_delta
    delta_y=0+(nn-1)*step_delta_y;
  E=0;
  E_diff_way=0;
sum_1_by_E_l=0;
sum_v_k_2=0;

for i=1:N_x


    for j=1:N_y

        l_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        l_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_l=-2*t*(cos(l_x*a)+cos(l_y*a))-mu;
        R_l=((e_l^2+(abs(delta_x))^2)^(1/2));
        sum_1_by_E_l=(sum_1_by_E_l+(1/R_l));


    end

end
sum_1_by_E_l_diff_deltas(mm,1)=(sum_1_by_E_l);


for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
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
[val,minimum_delta_at_g_max]=min(matrix_E_1_diff_way);
minimum_delta_at_g_max=step_delta_x*(minimum_delta_at_g_max-1);





%for g_real_space_mid

    g_real_space=(g_real_space_min+g_real_space_max)/2;
    g_real_space_mid=(g_real_space_min+g_real_space_max)/2;


parfor mm=1:n_delta
    delta_x=0+(mm-1)*step_delta_x;
    for nn=1:1%n_delta
    delta_y=0+(nn-1)*step_delta_y;
  E=0;
  E_diff_way=0;
sum_1_by_E_l=0;
sum_v_k_2=0;

for i=1:N_x


    for j=1:N_y

        l_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        l_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
        e_l=-2*t*(cos(l_x*a)+cos(l_y*a))-mu;
        R_l=((e_l^2+(abs(delta_x))^2)^(1/2));
        sum_1_by_E_l=(sum_1_by_E_l+(1/R_l));


    end

end
sum_1_by_E_l_diff_deltas(mm,1)=(sum_1_by_E_l);


for i=1:N_x
    for j=1:N_y
        
        k_x=(-pi/a)+(2*pi/(N_x*a))*(i-1);
        k_y=(-pi/a)+(2*pi/(N_y*a))*(j-1);
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
[val,minimum_delta_at_g_mid]=min(matrix_E_1_diff_way);
minimum_delta_at_g_mid=step_delta_x*(minimum_delta_at_g_mid-1);



if(abs(g_real_space_max-g_real_space_min)<tol)

    break;

end

if(minimum_delta_at_g_mid<=minimum_delta_at_g_max && minimum_delta_at_g_mid~=minimum_delta_at_g_min)

    g_real_space_max=g_real_space_mid;

end

if(minimum_delta_at_g_mid==minimum_delta_at_g_min)

    g_real_space_min=g_real_space_mid;

end







end
critical_g_diff_mu_th_limit(qq,1)=(g_real_space_max+g_real_space_min)/2;
%critical_g_diff_mu_dimensionless(qq,1)=critical_g_diff_mu(qq,1)/(-4*t-mu);


end


critical_g_diff_mu_th_limit_dimless=critical_g_diff_mu_th_limit/(8*t);


scatter(gap_x_axis_dimensionless,transpose(fliplr(transpose(critical_g_diff_mu_th_limit_dimless))))
xlabel('E_{gap}/Bandwidth');
ylabel('g_{c}/Bandwidth')
legend('Variational result for 6 \times 6 system','ED result for 6 \times 6 system','Variational result for 200 \times 200 system')
axes('Position',[0.25,0.6,0.2,0.2])
box on
scatter(gap_x_axis_dimensionless,transpose(fliplr(transpose(critical_g_diff_mu_th_limit_dimless-critical_g_diff_mu_dimless))))
xlabel('E_{gap}/Bandwidth')
ylabel({'Diff in g_{c}/Bandwidth b/w';' 200 \times 200 and 6 \times 6 system'})






figure
scatter(gap_x_axis_dimensionless,transpose(fliplr(transpose(critical_g_diff_mu_th_limit_dimless-critical_g_diff_mu_dimless))))
xlabel('E_{gap}/Bandwidth')
ylabel('Diff in g_{c} b/w 200 \times 200 and 6 \times 6 system')

toc


















