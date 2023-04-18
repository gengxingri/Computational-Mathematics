function [v] = solve_heat(psi,dt,eps)
N = 1001; 

% dt = .012
epsilon=eps*1i/2;




%(ik) and (ik)^2 vectors in x and y direction
kx=(1i*[0:N/2-1 0 -N/2+1:0]);
ky=(1i*[0:N/2-1 0 -N/2+1:0]');
k2x=kx.^2;
k2y=ky.^2;

[kxx,kyy]=meshgrid(k2x,k2y);
v_hat=fft2(psi);
vnew=(v_hat*(1/dt))./ ...
       (-(kxx+kyy)*epsilon+1/dt); %Implicit/Explicit timestepping
    %converts to real space in x-direction
v=ifft2(vnew);   
% mesh(abs(v))
end
