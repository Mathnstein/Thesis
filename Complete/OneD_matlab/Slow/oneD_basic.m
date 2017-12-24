%Bifurcation diagram for slow pass w/o osc

criteria=.7;

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\OneD_Basic')

%Initial values and set up
h=.01;
eps=.01;
start=1.5;stop=-1;
time=abs(stop-start)/eps;
muInit=start;
yInit=1-sqrt(1+muInit);

yDE=@(t,y,mu)(-mu+2*abs(y)-y*abs(y));
muDE=@(t,y,mu)(-eps);

%Numerical Solution
[~,y,mu]=RK2sys(yDE,muDE,yInit,muInit,0,time,h);

%Comparison

tipdelayed=eps*log(eps)/2;
tipxdel=tipdelayed*ones(1,100);
tipy=linspace(min(y),max(y),100);

tipactual=mu(find(y>criteria,1)-1);
dif=abs(tipdelayed-tipactual);
fprintf('The difference in predicted tip for eps=%f is %f\n',eps,dif)

%Figures

figure(1)

%plot(tipxdel,tipy,'b-.')
syms mu1 x
z=@(mu,y)(-mu+2*abs(y)-y*abs(y));
h1=ezplot(z,[-1,1.5,-.6,2.5]);
set(h1,'linestyle',':','color','k')
hold on
h2=ezplot(z,[-1,1.5,-.6,0]);
set(h2,'color','r','linewidth', 2)
h3=ezplot(z,[-1,1.5,1,2.5]);
set(h3,'color','r','linewidth', 2)
plot(mu,y,'k--','linewidth',2)
xlabel('\mu');ylabel('x')
title('')
%legend('Ramped Solution','Equilibrium')

print('-f1','slow_bif_diagram','-djpeg');

fileID = fopen('slow_bif_diagram_information.txt','w');
fprintf(fileID,'tipping actual =%f\n',tipactual);
fprintf(fileID,'tipping estimate =%f\n',tipdelayed);
fprintf(fileID,'epsilon=%f\n',eps);
fprintf(fileID,'Criteria=%f\', criteria);
fclose(fileID);

%Comparison plot
epsvec=linspace(.001,.1,20);
tipdelayedvec=zeros(1,length(epsvec));
tipactualvec=zeros(1,length(epsvec));

for i=1:length(epsvec)
    h=.01;
    eps=epsvec(i);
    start=1.5;stop=-1;
    time=abs(stop-start)/eps;
    muInit=start;
    yInit=1-sqrt(1+muInit);

    yDE=@(t,y,mu)(-mu+2*abs(y)-y*abs(y));
    muDE=@(t,y,mu)(-eps);

    %Numerical Solution
    [~,y,mu]=RK2sys(yDE,muDE,yInit,muInit,0,time,h);

    tipdelayedvec(i)=eps*log(eps);
    tipactualvec(i)=mu(find(y>criteria,1));
    
end
figure(2)
plot(epsvec,tipdelayedvec,'k','linewidth',2);
hold on
plot(epsvec,tipactualvec,'r*');
xlabel('\epsilon');ylabel('\mu-\mu_{ns}')


print('-f2','slow_epscomp','-djpeg');

