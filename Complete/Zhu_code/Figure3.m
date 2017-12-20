% Replicates some code from Zhu paper

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\Complete\Functions')

mu=.05;A=5;

lambda=0:.01:1;
a1=zeros(1,length(lambda));
h=.01;
steadysta1=zeros(length(lambda),1);
a0=.2;init1=.04;
w=zeros(length(lambda),1);

for i=1:length(lambda)
    a1(i)=mu^(2/3)*(-2.33810)+A^2*mu^(2*lambda(i))/2;
    w(i)=mu^(-lambda(i));
    xode=@(t,x,a)(a-x^2+A*sin(w(i)*t));
    aode=@(t,x,a)(-mu);
    [~,x,a]=RK2sys(xode,aode,init1,a0,0,200,h);
    I=find(x<=-10);
    steadysta1(i)=a(I(1));
end

b1=mu^(2/3)*(-2.33810)*ones(1,length(lambda));

mu=.01;A=5;

lambda=0:.01:1;
a2=zeros(1,length(lambda));
h=.01;
steadysta2=zeros(length(lambda),1);
a0=.1;init1=.01;
w=zeros(length(lambda),1);

for i=1:length(lambda)
    a2(i)=mu^(2/3)*(-2.33810)+A^2*mu^(2*lambda(i))/2;
    w(i)=mu^(-lambda(i));
    xode=@(t,x,a)(a-x^2+A*sin(w(i)*t));
    aode=@(t,x,a)(-mu);
    [~,x,a]=RK2sys(xode,aode,init1,a0,0,200,h);
    I=find(x<=-10);
    steadysta2(i)=a(I(1));
end

b2=mu^(2/3)*(-2.33810)*ones(1,length(lambda));

mu=.001;A=5;

lambda=0:.01:1;
a3=zeros(1,length(lambda));
h=.01;
steadysta3=zeros(length(lambda),1);
a0=.2;init1=.04;
w=zeros(length(lambda),1);

for i=1:length(lambda)
    a3(i)=mu^(2/3)*(-2.33810)+A^2*mu^(2*lambda(i))/2;
    w(i)=mu^(-lambda(i));
    xode=@(t,x,a)(a-x^2+A*sin(w(i)*t));
    aode=@(t,x,a)(-mu);
    [~,x,a]=RK2sys(xode,aode,init1,a0,0,500,h);
    I=find(x<=-10);
    steadysta3(i)=a(I(1));
end

b3=mu^(2/3)*(-2.33810)*ones(1,length(lambda));

figure(1)
p1=plot(a1(60:101),lambda(60:101),'r')
hold on
plot(steadysta1(60:101),lambda(60:101),'r-.')
plot(b1,lambda,'r--')

p2=plot(a2(50:101),lambda(50:101),'b')
plot(steadysta2(50:101),lambda(50:101),'b-.')
plot(b2,lambda,'b--')

p3=plot(a3(30:101),lambda(30:101),'k')
plot(steadysta3(30:101),lambda(30:101),'k-.')
plot(b3,lambda,'k--')
hold off
xlabel('a');
ylabel('\lambda');
legend([p1 p2 p3],'\mu=.05','\mu=.01','\mu=.001')

mu=.05;A=1;

lambda=0:.01:1;
a1=zeros(1,length(lambda));
h=.01;
steadysta1=zeros(length(lambda),1);
a0=-.2;init1=-.1;
w=zeros(length(lambda),1);

for i=1:length(lambda)
    a1(i)=mu^(2/3)*(-2.33810)+A^2*mu^(2*lambda(i))/2;
    w(i)=mu^(-lambda(i));
    xode=@(t,x,a)(a-x^2+A*sin(w(i)*t));
    aode=@(t,x,a)(-mu);
    [~,x,a]=RK2sys(xode,aode,init1,a0,0,200,h);
    I=find(x<=-10);
    steadysta1(i)=a(I(1));
end

mu=.01;A=1;

lambda=0:.01:1;
a2=zeros(1,length(lambda));
h=.01;
steadysta2=zeros(length(lambda),1);
a0=.13;init1=.1;
w=zeros(length(lambda),1);

for i=1:length(lambda)
    a2(i)=mu^(2/3)*(-2.33810)+A^2*mu^(2*lambda(i))/2;
    w(i)=mu^(-lambda(i));
    xode=@(t,x,a)(a-x^2+A*sin(w(i)*t));
    aode=@(t,x,a)(-mu);
    [~,x,a]=RK2sys(xode,aode,init1,a0,0,200,h);
    I=find(x<=-10);
    steadysta2(i)=a(I(1));
end

mu=.001;A=1;

lambda=0:.01:1;
a3=zeros(1,length(lambda));
b3=zeros(1,length(lambda));
h=.01;
steadysta3=zeros(length(lambda),1);
a0=.13;init1=.1;
w=zeros(length(lambda),1);

for i=1:length(lambda)
    a3(i)=mu^(2/3)*(-2.33810)+A^2*mu^(2*lambda(i))/2;
    b3(i)=A^2*mu^(2*lambda(i))/2;
    w(i)=mu^(-lambda(i));
    xode=@(t,x,a)(a-x^2+A*sin(w(i)*t));
    aode=@(t,x,a)(-mu);
    [~,x,a]=RK2sys(xode,aode,init1,a0,0,200,h);
    I=find(x<=-10);
    steadysta3(i)=a(I(1));
end


figure(2)
p1=plot(a1(40:101),lambda(40:101),'r')
hold on
plot(steadysta1(20:101),lambda(20:101),'r-.')
p2=plot(a2(20:101),lambda(20:101),'b')
plot(steadysta2(20:101),lambda(20:101),'b-.')
p3=plot(a3(10:101),lambda(10:101),'k')
plot(steadysta3(10:101),lambda(10:101),'k-.')
plot(b3(10:101),lambda(10:101),'k--')

hold off
xlabel('a');
ylabel('\lambda');
legend([p1,p2,p3], '\mu=.05','\mu=.01','\mu=.001')
