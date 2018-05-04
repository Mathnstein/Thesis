%Bifurcation diagram for slow pass w/o osc

criteria=.5;

% Load function folder
addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\OneD_Basic')

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

tipdelayed=eps*log(eps);
tipxdel=tipdelayed*ones(1,100);
tipy=linspace(min(y),max(y),100);

tipactual=mu(find(y>criteria,1)-1);
truetipvec=tipactual*ones(1,100);

%Figures
close(figure(1))
figure(1); % Full picture

syms mu1 x
z=@(mu,y)(-mu+2*abs(y)-y*abs(y));
h1=fimplicit(z,[-1,1.5,-.6,2.5]);
set(h1,'linestyle',':','color','k')
hold on
h2=fimplicit(z,[-1,1.5,-.6,0]);
set(h2,'color','r','linewidth', 2)
h3=fimplicit(z,[-1,1.5,1,2.5]);
set(h3,'color','r','linewidth', 2)
plot(mu,y,'k--','linewidth',2)
set(gca,'fontsize',14)
xlabel('\mu','fontsize',20)
ylabel('x','fontsize',20)
title('')

% Remove vectorized warning, not important for speed
[~,war]=lastwarn();
warning('off',war);
print('-f1','slow_bif_diagram','-djpeg');

%Zoom
plot(tipxdel, tipy, 'b--','linewidth',2)
plot(truetipvec, tipy, 'k','linewidth',2)
xlim([-.05 .1])
ylim([-.05 .2])

print('-f1','slow_bif_diagram_zoom','-djpeg');

%Comparison plot
epsvec=linspace(.001,.2,20);
tipdelayedvec=zeros(1,length(epsvec));
tipactualvec=zeros(1,length(epsvec));

for i=1:length(epsvec)
    h=.01;
    eps=epsvec(i);
    start=1;stop=-1;
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

close(figure(2))
figure(2)
plot(epsvec,tipdelayedvec,'k','linewidth',2);
hold on
plot(epsvec,tipactualvec,'r*');
set(gca,'fontsize',14)
xlabel('\epsilon','fontsize',20)
ylabel('\mu','fontsize',20)


print('-f2','slow_epscomp','-djpeg');