%Replicates some code from Zhu paper


% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\Complete\Functions')

a1=linspace(-.065,.45,10);
mu=.005;
lambda=.8;
A=1;

x1=zeros(1,length(a1));

for i=1:length(a1)
   x1(i)=-mu^(1/3)*airy(1,a1(i)*mu^(-2/3))/airy(0,a1(i)*mu^(-2/3)) ;
end

ode=@(t,a2)(-mu);
[t,a2]=RK2(ode,.45,0,103,.1);


x2=zeros(1,length(a2));

for i=1:length(t)
    x2(i)=mu^lambda*(-A*cos(mu^(-lambda)*t(i)))-mu^(1/3)*(airy(1,(a2(i)-mu^(2*lambda)*A^2/2)/mu^(2/3))/airy(0,(a2(i)-mu^(2*lambda)*A^2/2)/mu^(2/3)));
end

ode=@(t,a3)(-mu);
[t,a3]=RK2(ode,.45,0,88,.1);

x3=zeros(1,length(a3));

for i=1:length(t)
    x3(i)=sqrt(a3(i))+mu/(4*a3(i))-mu^lambda*(A*cos(mu^(-lambda)*t(i)));
end

figure(1)
plot(a1,x1,'r*')
hold on
plot(a3,x3,'g')
hold on
plot(a2,x2,'--b')
xlabel('a');
ylabel('x');
legend('(3)','(16)','(22)','location','northwest')

a1=linspace(-.065,.45,10);
mu=.005;
lambda=.2;
A=1;

x1=zeros(1,length(a1));

for i=1:length(a1)
   x1(i)=-mu^(1/3)*airy(1,a1(i)*mu^(-2/3))/airy(0,a1(i)*mu^(-2/3)) ;
end

ode=@(t,a2)(-mu);
[t,a2]=RK2(ode,.45,0,90,.1);


x2=zeros(1,length(a2));

for i=1:length(t)
x2(i)=mu^lambda*(-A*cos(mu^(-lambda)*t(i)))-mu^(1/3)*(airy(1,(a2(i)-mu^(2*lambda)*A^2/2)/mu^(2/3))/airy(0,(a2(i)-mu^(2*lambda)*A^2/2)/mu^(2/3)));
end

ode=@(t,a3)(-mu);
[t,a3]=RK2(ode,.45,0,88,.1);

x3=zeros(1,length(a3));

for i=1:length(t)
    x3(i)=sqrt(a3(i))+mu/(4*a3(i))-mu^lambda*(A*cos(mu^(-lambda)*t(i)));
end

figure(2)
plot(a1,x1,'r*')
hold on
plot(a3,x3,'g')
hold on
plot(a2,x2,'b--')

xlabel('a');
ylabel('x');
legend('(3)','(16)','(22)','location','northwest')

