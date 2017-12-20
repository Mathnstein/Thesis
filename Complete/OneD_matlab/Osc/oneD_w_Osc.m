%Full analysis of One D w/ Osc.
%(1)A dynamics plot and a (2)comparison plot against frequency will be produced

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\Complete\Functions')
 
%Analysis

%Set up step size and parameter grid
h=.01;
A=1;
invOmega=10;
start=1;stop=-1;
yInit=1-sqrt(1+start);

 
%Fixed mu
muper=fliplr(linspace(stop,start,1500));
N=length(muper);
xfinal=zeros(1,N);
for i=1:N
    t1=0;
    t2=30;
    xinit=1-sqrt(1+muper(i));
    xDE=@(t,x)(-muper(i)+2*abs(x)-x*abs(x)+A*sin(invOmega*t));
    [t,x]=RK2(xDE,xinit,t1,t2,h);
    xfinal(i)=x(t2/h-N+i);
end
 
 
%Search for bifurcation
n=3000;
range=linspace(0,3,n);

%Create range of candidate values
potentialval=fliplr(linspace(0,A,n));

%Create a logic loop to find bifurcation
truth=zeros(1,n+1);
truth(1)=1;

for i=1:n
    if truth(i)==1
        fun=@(x)(abs(-A*cos(x)+potentialval(i)));
        equation=@(m)(-m+(1/pi)*integral(fun,0,2*pi));
        [m,root]=root_algorithm(equation,range);
        %Condition to continue or stop at bif
         if isa(root,'double')
            truth(i+1)=1;
         end
    else
    end
    
end

bif=m/invOmega;
 
tipactual=muper(find(xfinal>1,1));
tipactualper=tipactual*ones(1,100);
dif=abs(tipactual-bif);
fprintf('The difference in predicted tip for fixed parameter is %f\n',dif)
 
 
mulower=linspace(0,1.5,150);
lower=@(mu)(1-sqrt(1+mu));
xlower=lower(mulower);
 
muupper=linspace(-1,1,200);
upper=@(mu)(1+sqrt(1-mu));
xupper=upper(muupper);
 
mumid=linspace(0,1,100);
middle=@(mu)(1-sqrt(1-mu));
xmid=middle(mumid);
 
 
%(1) Dynamics plot
figure(1)
%plot(mu,y,'r-.')
plot(muper,xfinal,'k--','linewidth',2)
hold on
%plot(tipxdel,tipy,'r')
xlim([stop,start])
plot(mulower,xlower,'r','linewidth',2)
hold on
plot(mumid,xmid,'k:')
plot(muupper,xupper,'r','linewidth',2)
xlabel('\mu')
ylabel('x')
 
%{

%(2) Comparison plot
invOmegavec=linspace(.001,.3,20);
bifvec=zeros(1,length(invOmegavec));
bifactualvec=zeros(1,length(invOmegavec));
M=length(invOmegavec);
 
for i=1:length(invOmegavec)
    invOmega=invOmegavec(i);
    %Make sure to have a time step that will capture highly oscillatory
    %pieces
    h=2*pi*invOmega/(20);
    start=1;stop=-1;
    %Have a parameter grid fine enough to see the different bifurcation
    %locations for highly oscillatory cases
    muper=fliplr(linspace(stop,start,1000));
    N=length(muper);
    yfinal=zeros(1,N);
     
    %Numerical Solution
    for j=1:N
        yInit=1-sqrt(1+muper(j));
        yDE=@(t,y)(-muper(j)+2*abs(y)-y*abs(y)+A*sin(t/invOmega));
        [t,y]=RK2(yDE,yInit,0,15,h);
        yfinal(j)=y(end);
    end
 
    bifvec(i)=4*abs(A)*invOmega/(pi);
    cutoff=1;
    bifactualvec(i)=muper(find(yfinal>cutoff,1));
end
 
figure(2)
plot(invOmegavec,bifvec,'k','linewidth',2);
hold on
plot(invOmegavec,bifactualvec,'r*');
xlabel('\Omega^{-1}');ylabel('\mu-\mu_{ns}')

%}
