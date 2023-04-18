%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Solving 2D Allen-Cahn Eq using pseudo-spectral with Implicit/Explicit
%u_t= epsilon(u_{xx}+u_{yy}) 
%where u-u^3 is treated explicitly and epsilon(u_{xx} + u_{yy}) is treated implicitly
%BC = Periodic
%IC=v=sin(2*pi*x)+0.001*cos(16*pi*x);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%Grid
eps=1;gamma_y=1;kappa=-2;delta_t=0.005;
% eps=0.3;gamma_y=1;kappa=-1.9718;delta_t=0.0005;


%time steps
time_step=100;


N = 1001; h = 20/(N-1); epsilon=eps*1i/2;
x=linspace(-10,10,N);
dt = delta_t;

% x and y meshgrid
y=x';
[xx,yy]=meshgrid(x,y);


% initial conditions
exp_mat=exp(-(xx.^2+yy.^2)/(2*eps));
v=1/sqrt(pi*eps)*exp_mat;


% (ik) and (ik)^2 vectors in x and y direction
kx=(1i*[0:N/2-1 0 -N/2+1:0]);
ky=(1i*[0:N/2-1 0 -N/2+1:0]');
k2x=kx.^2;
k2y=ky.^2;

[kxx,kyy]=meshgrid(k2x,k2y);
     
value=zeros(1,time_step+1);
value(1,1)=(abs(v(501,501))).^2;   

for n = 1:time_step
    v_nl=v.^3;  %calculates nonlinear term in real space
    %FFT for linear and nonlinear term
    v_nl = 2/eps*fft2(v_nl);
    v_v=(-1i)/(2*eps)*fft2(v.*(xx.^2+yy.^2));
    v_hat=fft2(v);
    vnew=(v_hat*(1/dt)+v_nl+v_v)./ ...
       (-(kxx+kyy)*epsilon+1/dt); %Implicit/Explicit timestepping
    %converts to real space in x-direction
    v=ifft2(vnew);   
    value(1,n+1)=(abs(v(501,501))).^2;
    %Plots each timestep
    meshc(abs(v).^2); title(['Time ',num2str(n)]); axis([0 N 0 N 0 0.4]); 
    xlabel x; ylabel y; zlabel \Psi^2;
    view(43,22); drawnow;  
end
