tic
clear
eps=1;gamma_y=1;kappa=-2;
% eps=0.3;gamma_y=1;kappa=-1.9718;
h=1/50;delta_t=0.0005;
time_step=1000;

%create grid
x=-10:h:10;
y=x;

M=length(x);

%initial condition
[X,Y]=meshgrid(x,y);  
exp_mat=exp(-(X.^2+Y.^2)/(2*eps));
Psi=1/sqrt(pi*eps)*exp_mat;
% Psi(:,:,time_step+1)=0;%extend timestamp



value=zeros(1,time_step+1);
value(1,1)=(abs(Psi(501,501))).^2;
% value(1,1)=(abs(Psi(501,501,1))).^2;
% sigma_x=zeros(1,time_step+1);
% sigma_x(1,1)=sqrt(inner_product((X-inner_product(Psi(:,:,1),x)).^2,x));
%TSSP iteration
i=sqrt(-1);
for step=1:time_step
%     tempt=Psi(:,:,step);%Psi^n
    tempt=Psi;
    %ode n-->*
    psi_star=tempt.*exp(-(X.^2+Y.^2+kappa*(abs(tempt)).^2)*i*delta_t/(2*eps));%psi^*
    %heat equation *--->**
    psi_star_star=solve_heat(psi_star,delta_t,eps);
    
    %ode **--->n+1
    Psi=psi_star_star.*exp(-(X.^2+Y.^2+kappa*(abs(psi_star_star)).^2)*i*delta_t/(2*eps));
%     sigma_x(1,step+1)=sqrt(inner_product((X-inner_product(Psi(:,:,step+1),x)).^2,x));
    value(1,step+1)=(abs(Psi(501,501))).^2;%value at (0,0)
    
end
meshc(abs(Psi).^2)
toc