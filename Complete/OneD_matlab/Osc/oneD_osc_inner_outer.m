%This plots the diagram showing the regions

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\Complete\Functions')

%Note: Test is different values of T, causes different bifurcation
%locations

A=2;Omega=100;
mu=linspace(0,10/Omega,200);
N=length(mu);
x=zeros(N,1);
C=sqrt(abs(A)/(2*(1-2/pi)));
boundary=2*abs(A)/Omega;
bif=4*abs(A)/(pi*Omega);

t=linspace(0,1,N);
for i=1:N
    outer=@(T)(1-sqrt(1+mu(i))-A*cos(Omega*T)/Omega);
    innerC1=@(T)(-mu(i)/2-A*cos(Omega*T)/Omega);
    innerC2=@(T)(-C*sqrt((mu(i)-4*abs(A)/(pi*Omega))/Omega)-A*cos(Omega*T)/Omega);

    if mu(i)>= boundary
        x(i)=innerC1(t(i));
    elseif mu(i)>= bif
        x(i)=innerC2(t(i));
    else
        x(i)=2;
    end

end

critvec=boundary*ones(N,1);
critrange=linspace(-8/Omega,6/Omega,N);

bifvec=bif*ones(N,1);

lower=@(mu)(1-sqrt(1+mu));
xlower=lower(mu);


middle=@(mu)(1-sqrt(1-mu));
xmid=middle(mu);

figure(1)
plot(mu,x,'k--')
hold on
plot(mu,xlower,'r','linewidth',2)
plot(mu,xmid,'k-.')
plot(critvec,critrange,'b')
plot(bifvec,critrange,'g')
xlabel('\mu');ylabel('x')
axis([0,.1,-.08,.06])
