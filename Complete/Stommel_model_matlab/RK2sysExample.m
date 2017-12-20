%This is a script to run RK2 for a 2x2 system
 

%For bottom

%Set up step size and parameter grid
h=.1;
n=100;

%Initial values and set up
check=false;
while check==false
init1=2.9;init2=2.9;t1=0;t2=100;L=length(t1:h:t2);
eta1=3;eta3=.3;
eta2bot=linspace(0,2,n);
steadystbot=zeros(1,n);
ode1=@(t,y1,y2)(eta1-y1*(1+abs(y1-y2)));
ode2=@(t,y1,y2)(eta2bot(n)-y2*(eta3+abs(y1-y2)));

%Numerical accuracy evaluation
TOL=.1; mult=2;
[~,y1,y2]=RK2sys(ode1,ode2,init1,init2,t1,t2,h);
[~,k1,k2]=RK2sys(ode1,ode2,init1,init2,t1,t2,h/mult);
nodetest1=zeros(1,L);
nodetest1(1)=y1(1);
nodetest1(L)=y1(L);
nodetest2=zeros(1,L);
nodetest2(1)=y2(1);
nodetest2(L)=y2(L);

for i=2:L-1
    nodetest1(i)=k1(mult*(i-1));
    nodetest2(i)=k2(mult*(i-1));
end
if [abs(y1-nodetest1)<TOL ; abs(y2-nodetest2)<TOL]
    check=true;
    disp('The step size that is accurate enough for bottom is')
    disp(h)
else
    disp(h);disp('is not small enough.')
    h=h/mult;
end

end

for i=1:length(eta2bot)
ode1=@(t,y1,y2)(eta1-y1*(1+abs(y1-y2)));
ode2=@(t,y1,y2)(eta2bot(i)-y2*(eta3+abs(y1-y2)));
[~,y1,y2]=RK2sys(ode1,ode2,init1,init2,t1,t2,h);
psi=y1-y2;
steadystbot(i)=psi(n);
end
 
%For top
 
%Set up step size and parameter grid
h=.1;
n=100;

%Initial values and set up
check=false;
while check==false
init1=1;init2=1;t1=0;t2=100;L=length(t1:h:t2);
eta1=3;eta3=.3;
eta2top=linspace(0,2,n);
steadystup=zeros(1,n);
ode1=@(t,y1,y2)(eta1-y1*(1+abs(y1-y2)));
ode2=@(t,y1,y2)(eta2top(n)-y2*(eta3+abs(y1-y2)));

%Numerical accuracy evaluation
TOL=.1; mult=2;
[~,y1,y2]=RK2sys(ode1,ode2,init1,init2,t1,t2,h);
[~,k1,k2]=RK2sys(ode1,ode2,init1,init2,t1,t2,h/mult);
nodetest1=zeros(1,L);
nodetest1(1)=y1(1);
nodetest1(L)=y1(L);
nodetest2=zeros(1,L);
nodetest2(1)=y2(1);
nodetest2(L)=y2(L);

for i=2:L-1
    nodetest1(i)=k1(mult*(i-1));
    nodetest2(i)=k2(mult*(i-1));
end
if [abs(y1-nodetest1)<TOL ; abs(y2-nodetest2)<TOL]
    check=true;
    disp('The step size that is accurate enough for top is')
    disp(h)
else
    disp(h);disp('is not small enough.')
    h=h/mult;
end

end

for i=1:length(eta2top)
ode1=@(t,y1,y2)(eta1-y1*(1+abs(y1-y2)));
ode2=@(t,y1,y2)(eta2top(i)-y2*(eta3+abs(y1-y2)));
[t,y1,y2]=RK2sys(ode1,ode2,init1,init2,t1,t2,h);
psi=y1-y2;
steadystup(i)=psi(n);
end

figure(1)
plot(eta2bot,steadystbot,'r')
hold on
plot(eta2top,steadystup,'b')
xlabel('eta2')
ylabel('steady state')
hold off
eta1=3;eta3=.3; eta2=1;
figure(2)
syms t s
curve=ezplot(eta1-t*(1+abs(t-s))==eta2-s*(eta3+abs(t-s)));
set(curve,'Color','red');