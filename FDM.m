clear
tic
eps=1;gamma_y=1;kappa=-2;
h=1/50;delta_t=0.0005;
time_step=80000;

%create grid
x=-10:h:10;
y=x;

%initial condition
[X,Y]=meshgrid(x,y);
exp_mat=exp(-(X.^2+Y.^2)/(2*eps));
Psi=1/sqrt(pi*eps)*exp_mat;
% Psi(:,:,time_step+1)=0;%extend timestamp


value=zeros(1,time_step+1);
value(1,1)=(abs(Psi(501,501))).^2;
%explicit iteration
i=sqrt(-1);
for step=1:time_step
%     tempt=Psi(:,:,step);%current matrix
    tempt=Psi;
    %four point arount i,j,at current time level
    psi1=zeros(size(tempt));
    psi1(:,1:end-1)=tempt(:,1:end-1);
    psi1(:,end)=0;
    
    psi2=zeros(size(tempt));
    psi2(:,2:end)=tempt(:,2:end);
    psi2(:,1)=0;
    
    psi3=zeros(size(tempt));
    psi3(1:end-1,:)=tempt(1:end-1,:);
    psi3(end,:)=0;
    
    psi4=zeros(size(tempt));
    psi4(2:end,:)=tempt(2:end,:);
    psi4(1,:)=0;
    %iteration
    %Psi(:,:,step+1)=((2+0.5*(X.^2+Y.^2)-2*tempt.^2).*tempt-0.5*(psi1+psi2+psi3+psi4))*delta_t/i+tempt;
    Psi=((2+0.5*(X.^2+Y.^2)-2*tempt.^2).*tempt-0.5*(psi1+psi2+psi3+psi4))*delta_t/i+tempt;
    value(1,step+1)=(abs(Psi(501,501))).^2;%value at (0,0)
    %     meshc(abs(Psi(:,:,step+1)).^2);title(['Time ',num2str(step)]);
%     view(43,22); drawnow;  
end
toc