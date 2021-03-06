% Load function folder
addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\TwoD_SlowOsc')

% Determine whether or not to zoom 1 Yes, 0 No
zoom = 0;

% Start parallel pooling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Code Begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start on NS stable branch
% Works for n3<1

% Hyperparameters
eta1 = 4; eta3 = .375;
eps = .01;
lambda = .8;
Omega = eps^(-lambda);
A=2;
B=2;

% Locate the smooth bifurcation
derivative = [-2 -(eta3+4) -2*(eta3+1) eta1-eta3*(eta1+1)];
r = roots(derivative);
Vsmooth = r(real(r)>=0&imag(r)==0);
smootheta2 = Vcurve(Vsmooth,eta1,eta3);

criteria = Vsmooth;

%-------------------------------------------------
% System 
start = eta1*eta3-1; stop=eta1*eta3+1;
time = (stop-start)/(eps);

% Non-dim Stommel Equations
vDE = @(t,V,T,eta2)((eta1-eta2)-V*abs(V)-T+eta3*(T-V)+A*sin(Omega*t));
tDE = @(t,V,T,eta2)(eta1-T*(1+abs(V))+B*sin(Omega*t));
eta2DE = @(t,V,T,eta2)(-eps);

eta2Init = stop;
vInit = -stop/eta1+eta3;
tInit = eta1/(1-vInit);
h = min(.01,2*pi/(10*Omega));

% Solve equations for numerical soln, step size h
[~,Vnum,Tnum,eta2num] = RK2sys3(vDE,tDE,eta2DE,vInit,tInit,eta2Init,0,time,h);

% Ignore transient behavior
Vnum = Vnum(100:end);
Tnum = Tnum(100:end);
eta2num = eta2num(100:end);

Min = min(Vnum);
Max = max(Vnum);

n=100;
% Tipping
Mat = [eta3 1-eta3;eta1 1];
tipslow = eta1*eta3-eps*log(eps)/min(eig(Mat));
tippredlargevec = tipslow*ones(1,n);

scaledcoef = eps^(lambda-1)*A;
c0 = 2*eta1*(1-eta3)*scaledcoef/pi;
c1 = eta1*(1-eta3)/(pi*scaledcoef);
tippredRED = eta1*eta3+eps*(c0-eta3^2/(4*c1)-c1^(-1/3)*2.3381);


yvec = linspace(Min,Max,n);
tipactual=eta2num(find(Vnum>criteria,1));
tipactualvec=tipactual*ones(1,n);
tippredREDvec = tippredRED*ones(1,n);
tippredslowvec = tipslow*ones(1,n);


%--------------------Plot the dynamics plot--------------------------

close(figure(1))
figure(1)
z = @(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h1 = fimplicit(z,[eta1*eta3-1,eta1*eta3+1,-.5,1.5]);
set(h1,'linestyle','-.','color','k')
hold on
z = @(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h2 = fimplicit(z,[eta1*eta3-1,eta1*eta3+1,-.5,0]);
set(h2,'color','r','linewidth',2)
hold on
z = @(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h3 = fimplicit(z,[eta1*eta3-1,eta1*eta3+1,Vsmooth,1.5]);
set(h3,'color','r','linewidth',2)

% Remove vectorized warning, not important for speed
[~,war]=lastwarn();
warning('off',war);

set(gca,'fontsize',18)
xlabel('\eta_2','fontsize',32)
ylabel('V','fontsize',32)
title('')
hold on
plot(eta2num,Vnum,'k--','linewidth',2)
xlim([start,stop])

% Zoom
if zoom ==1
    if lambda>1.5
        print('-f1','slowosc_bif_diagram_large','-djpeg');
        plot(tippredslowvec,yvec,'g--','linewidth',2)
        plot(tipactualvec,yvec,'k','linewidth',2)
        plot(tippredREDvec,yvec,'b--','linewidth',2)
        xlim([tipactual-.05 tipactual+.1])
        ylim([-.05 .2])
        print('-f1','slowosc_bif_diagram_large_zoom','-djpeg');
    elseif lambda>1
        print('-f1','slowosc_bif_diagram_medium','-djpeg');
        plot(tippredslowvec,yvec,'g--','linewidth',2)
        plot(tipactualvec,yvec,'k','linewidth',2)
        plot(tippredREDvec,yvec,'b--','linewidth',2)
        xlim([tipactual-.05 tipactual+.15])
        ylim([-.05 .2])
        print('-f1','slowosc_bif_diagram_medium_zoom','-djpeg');
    else
        print('-f1','slowosc_bif_diagram_small','-djpeg');
        plot(tippredREDvec,yvec,'b--','linewidth',2)
        plot(tipactualvec,yvec,'k','linewidth',2)
        xlim([eta1*eta3 eta1*eta3+.2])
        ylim([-.1 .2])
        print('-f1','slowosc_bif_diagram_small_zoom','-djpeg');
    end
end

% T plot
m=500;
Vlow=linspace(-1,0,m);
Vmid=linspace(0,Vsmooth,m);
Vup = linspace(Vsmooth,1.5,m);
Tlow=zeros(m,1);
Tmid = zeros(m,1);
Tup = zeros(m,1);

func=@(V)(eta1/(1+abs(V)));
for i=1:m
    Tlow(i) = func(Vlow(i));
    Tmid(i) = func(Vmid(i));
    Tup(i) = func(Vup(i));
end

close(figure(2))
figure(2)
plot(Vlow,Tlow,'r','linewidth',2)
hold on
plot(Vmid,Tmid, 'k-.')
plot(Vup,Tup,'r','linewidth',2)
plot(Vnum,Tnum,'k--','linewidth',2)
set(gca,'fontsize',18)
xlabel('V','fontsize',32)
ylabel('T','fontsize',32)
xlim([-.75 Max])
ylim([1.5 eta1])

% Zoom
if zoom ==1
    if lambda>1.5
        print('-f2','slowosc_Tplot_large','-djpeg');
        xlim([-.05 .1])
        ylim([3.8 4])    
        print('-f2','slowosc_Tplot_large_zoom','-djpeg');
    elseif lambda >1
        print('-f2','slowosc_Tplot_medium','-djpeg');
        xlim([-.05 .1])
        ylim([3.8 4])    
        print('-f2','slowosc_Tplot_medium_zoom','-djpeg');
    else
        print('-f2','slowosc_Tplot_small','-djpeg');
        xlim([-.07 .12])
        ylim([3.75 4])
        print('-f2','slowosc_Tplot_small_zoom','-djpeg');
    end
end

fprintf('The tipping is eta2= %f\n',tipactual)
fprintf('The reduced estimated tipping is eta2 = %f\n',tippredRED)


% Calculates the V curve
function eta2 = Vcurve(V,eta1,eta3)
    T0 = eta1/(1+abs(V));
    eta2 = eta1+eta3*(T0-V)-T0-V*abs(V);
end
