% Plots for oneD_osc

criteria =.5;

% Load function folder
addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\OneD_w_Osc')
 
% Time Series

%Set up sep size and parameter grid
h=.01;
A=2;
Omega=10;
start=1;stop=-1;
yInit=1-sqrt(1+start);
bif=4*abs(A)/(pi*Omega);
region1 = 2*A/Omega;

bifvec = bif*ones(1,100);
region1vec = region1*ones(1,100);
yvec = linspace(-1,3,100);

 
%Fixed mu
muper=[region1+4/Omega (bif+region1)/2 bif-1/Omega ];
N=length(muper);
close(figure(1))
labs ={};
c = {};
for i=1:N
    t1=0;
    t2=20;
    xinit=1-sqrt(1+muper(i));
    xDE=@(t,x)(-muper(i)+2*abs(x)-x*abs(x)+A*sin(Omega*t));
    [t,x]=RK2(xDE,xinit,t1,t2,h);
    figure(1)
    hold on
    P1 = plot(t,x,'linewidth',2);
    set(gca,'Xdir','reverse')
    c{i} = get(P1,'Color');
    if i == 1
        xbefore = x(500:end);
        labs{i}='\mu in Case I';
    elseif i == 2
        xat = x(500:end);
        labs{i}='\mu in Case II';
    elseif i == 3
        xafter = x(500:end);
        labs{i}='\mu after bifurcation';
    end
end
set(gca,'fontsize',18)
xlabel('Time','fontsize',32)
ylabel('x','fontsize',32)

print('-f1','osc_timeseries','-djpeg')

% Phase plot
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
close(figure(2))
figure(2)

xlim([stop,start])
plot(mulower,xlower,'r','linewidth',2)
hold on
plot(mumid,xmid,'k:')
plot(muupper,xupper,'r','linewidth',2)
set(gca,'fontsize',18)
xlabel('\mu','fontsize',32)
ylabel('x','fontsize',32)

% x=0 axis
plot(linspace(-1,1,100),zeros(1,100),'k')

% Parameter ranges
plot(bifvec, yvec, 'b--','linewidth',2)
plot(region1vec, yvec, 'g--','linewidth',2)
xlim([0 .6])
ylim([-.5 .4])

print('-f2','osc_cases','-djpeg')

%{
% Phase plane time series
plot(muper(1)*ones(1,length(xbefore)),xbefore,'color',c{1},'linewidth',2)
plot(muper(2)*ones(1,length(xat)),xat,'color',c{2},'linewidth',2)
plot(muper(3)*ones(1,length(xafter)),xafter,'color',c{3},'linewidth',2)
xlim([-.5 1])
ylim([-.5 2.2])
print('-f2','osc_bif_diagram','-djpeg')

xlim([0 .85])
ylim([-.5 .4])
print('-f2','osc_bif_diagram_zoom','-djpeg')
%}
