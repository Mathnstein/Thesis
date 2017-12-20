%This is a script to run RK2, pulls final values of an RK2 scheme for an
%ODE with a certain parameter and then takes the next parameter. Plots
%what this looks like.

%Set up step size and parameter grid
h=1;
n=1000;


%Initial values and set up
check=false;
while check==false
init=0;t1=0;t2=100;L=length(t1:h:t2);
a=linspace(0,16,n);
steadyst=zeros(1,n);
ode=@(t,y)(a(n)-y^2);

%Numerical accuracy evaluation
TOL=.75; mult=2;
[~,y]=RK2(ode,init,t1,t2,h);
[~,k]=RK2(ode,init,t1,t2,h/mult);
nodetest=zeros(1,L);
nodetest(1)=y(1);
nodetest(L)=y(L);
for i=2:L-1
    nodetest(i)=k(mult*(i-1));
end
if abs(y-nodetest)<TOL
    check=true;
    disp('The step size that is accurate enough is')
    disp(h)
else
    disp(h)
    disp('is not small enough.')
    h=h/mult;
end

end

%Create vector of steady states
for i=1:length(a)
ode=@(t,y)(a(i)-y^2);
[~,y]=RK2(ode,init,t1,t2,h);
steadyst(i)=y(n);
end

%Plot the steady states against the parameter values.
plot(a,steadyst)
xlabel('a')
ylabel('steady state')
