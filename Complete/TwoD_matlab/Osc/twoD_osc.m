% Load function folder
addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\TwoD_Osc')

% Zoom Yes == 1, No == 0

zoom=1;

eta1 = 4; eta3 = 3/8;
nsbif = eta1*eta3;
A = 2; B = 2; Omega = 10;
h = min(.01, 2*pi/(10*Omega));

caseIboundary = abs(A)*(eta1*(1-eta3)+eta3)/Omega+eta1*eta3;
mosc = 2*eta1*(1-eta3)/(pi*abs(A))-pi*abs(A)*eta3^2/(4*eta1*(1-eta3));
muosc = mosc/Omega+nsbif;

% Locate the smooth bifurcation
derivative = [-2 -(eta3+4) -2*(eta3+1) eta1-eta3*(eta1+1)];
r = roots(derivative);
Vsmooth = r(imag(r)==0);
smootheta2 = Vcurve(Vsmooth,eta1,eta3);

% Time series
close(figure(1))
close(figure(2))
eta2vec = [caseIboundary+2/Omega (caseIboundary+muosc)/2 eta1*eta3+1/Omega];
N = length(eta2vec);
labs={};
c={};
for i = 1:N
    vDEosc = @(t,V,T)((eta1-eta2vec(i))-V*abs(V)-T+eta3*(T-V)+A*sin(Omega*t));
    tDEosc = @(t,V,T)(eta1-T*(1+abs(V))+B*sin(Omega*t));
    
    vInit = -eta2vec(i)/eta1+eta3;
    tInit = eta1/(1-vInit);
    
    [t,Vosc,Tosc] = RK2sys(vDEosc,tDEosc,vInit,tInit,0,20,h);
    figure(1)
    hold on
    P1 = plot(t,Vosc,'linewidth',2);
    set(gca,'Xdir','reverse')
    set(gca,'fontsize',18)
    c{i} = get(P1,'Color');
    if i == 1
        Vbefore = Vosc(500:end);
        labs{i}='\eta_2 in Case I';
    elseif i == 2
        Vat = Vosc(500:end);
        labs{i}='\eta_2 in Case II';
    elseif i == 3
        Vafter = Vosc(500:end);
        labs{i}='\eta_2 after bifurcation';
    end
    xlabel('Time','fontsize',32)
    ylabel('V','fontsize',32)
    ylim([-.6 1.3])
    print('-f1','osc_Vtimeseries','-djpeg')
    
    figure(2)
    hold on
    P2 = plot(t,Tosc,'linewidth',2);
    set(gca,'Xdir','reverse')
    set(gca,'fontsize',18)
    c{i} = get(P2,'Color');
     if i == 1
        Tbefore = Tosc(500:end);
        labs{i}='\eta_2 in Case I';
    elseif i == 2
        Tat = Tosc(500:end);
        labs{i}='\eta_2 in Case II';
    elseif i == 3
        Tafter = Tosc(500:end);
        labs{i}='\eta_2 after bifurcation';
    end
    xlabel('Time','fontsize',32)
    ylabel('T','fontsize',32)
    print('-f2','osc_Ttimeseries','-djpeg')
end

yvec = linspace(-1,2,300);
caseIvec = caseIboundary* ones(1,300);
bifpredvec = muosc* ones(1,300);


% Phase plot time series
close(figure(3))
figure(3)
z = @(eta2,V)(eta1-eta2-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V));
h1 = fimplicit(z,[nsbif-1,nsbif+1,-.5,1.5]);
set(h1,'linestyle','-.','color','k')
hold on
h2 = fimplicit(z,[nsbif-1,nsbif+1,-.5,0]);
set(h2,'color','r','linewidth',2)
hold on
h3 = fimplicit(z,[nsbif-1,nsbif+1,Vsmooth,1.5]);
set(h3,'color','r','linewidth',2)

% Remove vectorized warning, not important for speed
[~,war] = lastwarn();
warning('off',war);
set(gca,'fontsize',18)
xlabel('\eta_2','fontsize',32)
ylabel('V','fontsize',32)
title('')

% x=0 axis
plot(linspace(nsbif-1,nsbif+1,100),zeros(1,100),'k')

% Parameter ranges
plot(bifpredvec, yvec, 'b--','linewidth',2)
plot(caseIvec, yvec, 'g--','linewidth',2)
xlim([nsbif 2.3])
ylim([-.5 .3])

%print('-f3','osc_cases','-djpeg')

% Phase plane time series
plot(eta2vec(1)*ones(1,length(Vbefore)),Vbefore,'color',c{1},'linewidth',2)
plot(eta2vec(2)*ones(1,length(Vat)),Vat,'color',c{2},'linewidth',2)
plot(eta2vec(3)*ones(1,length(Vafter)),Vafter,'color',c{3},'linewidth',2)
xlim([nsbif-.5 nsbif+1])
ylim([-.6 1.2])
print('-f3','osc_bif_diagram','-djpeg')

if zoom==1
    xlim([nsbif 2.3])
    ylim([-.5 .3])
    print('-f3','osc_bif_diagram_zoom','-djpeg');
end

% T phase plot
m = 300;
Vlow = linspace(-1,0,m);
Vmid = linspace(0,Vsmooth,m);
Vup = linspace(Vsmooth,1.5,m);


func=@(x)(eta1./(1+abs(x)));

Tlow = func(Vlow);
Tmid = func(Vmid);
Tup = func(Vup);

close(figure(4))
figure(4)
plot(Vlow,Tlow,'r','linewidth',2)
hold on
plot(Vmid,Tmid, 'k-.')
plot(Vup,Tup,'r','linewidth',2)
plot(Vbefore,Tbefore,'color',c{1},'linewidth',2)
plot(Vat,Tat,'color',c{2},'linewidth',2)
plot(Vafter,Tafter,'color',c{3},'linewidth',2)
set(gca,'fontsize',18)
xlabel('V','fontsize',32)
ylabel('T','fontsize',32)

print('-f4','osc_bif_Tplot','-djpeg');

if zoom == 1
    xlim([-.5 .3])
    ylim([2.7 eta1])
    print('-f4','osc_bif_Tplot_zoom','-djpeg');
end


% Calculates the V curve
function eta2 = Vcurve(V,eta1,eta3)
    T0 = eta1/(1+abs(V));
    eta2 = eta1+eta3*(T0-V)-T0-V*abs(V);
end