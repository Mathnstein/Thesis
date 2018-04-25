% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\TwoD_Osc')

% Zoom: 0 = No, 1 = Yes
zoom = 1;

% Comparison: 0 = No, 1 = Yes
comp = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Code Begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start on NS stable branch
%Works for n3<1
%Parameters
n = 2000;
eta1 = 4; eta3 = 3/8;
nsbif = eta1*eta3;
A = 1; B = 1; Omega = 30;
h = min(.01, 2*pi/(10*Omega));

%Fixed parameter system
start = nsbif-.5; stop = nsbif+1;
eta2 = fliplr(linspace(start,stop,n));
vOscTimeSeries = zeros(1,length(eta2));
tOscTimeSeries = vOscTimeSeries;

%Create equilibrium vector
parfor i = 1:length(eta2)
    %Non-dim Stommel Equations    
    vDEosc = @(t,V,T)((eta1-eta2(i))-V*abs(V)-T+eta3*(T-V)+A*sin(Omega*t));
    tDEosc = @(t,V,T)(eta1-T*(1+abs(V))+B*sin(Omega*t));
    
    vInit = -eta2(i)/eta1+eta3;
    tInit = eta1/(1-vInit);
    
    
    %Solve equations for long term equilibrium, step size h
    [~,Vosc,Tosc] = RK2sys(vDEosc,tDEosc,vInit,tInit,0,50,h);
    
    %Keep steady state
    vOscTimeSeries(i) = Vosc(floor(50/h)-n+i);
    tOscTimeSeries(i) = Tosc(floor(50/h)-n+i);
end

% Locate the smooth bifurcation
derivative = [-2 -(eta3+4) -2*(eta3+1) eta1-eta3*(eta1+1)];
r = roots(derivative);
Vsmooth = r(imag(r)==0);
smootheta2 = Vcurve(Vsmooth,eta1,eta3);

%-------------------------------------------------
%Osc system bif values

criteria = Vsmooth;

caseIboundary = abs(A)*(eta1*(1-eta3)+eta3)/Omega+eta1*eta3;
mosc = 2*eta1*(1-eta3)/(pi*abs(A))-pi*abs(A)*eta3^2/(4*eta1*(1-eta3));
muosc = mosc/Omega+nsbif;


bifactual = eta2(find(vOscTimeSeries>criteria,1));

Min = min(vOscTimeSeries);
Max = max(vOscTimeSeries);

Vaxis = zeros(1,length(eta2));

yvec = linspace(Min,Max,300);
caseIvec = caseIboundary* ones(1,300);
bifpredvec = muosc* ones(1,300);
bifactualvec = bifactual*ones(1,300);


error = abs(bifactual-muosc);

%--------------------Plot the dynamics plot--------------------------
close(figure(1))
figure(1)
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

xlabel('\eta_2')
ylabel('V')
title('')
hold on
plot(eta2,vOscTimeSeries,'k--','linewidth',2)
xlim([start,stop])
%title(['V Equilibria for \epsilon=' num2str(eps)])
xlabel('\eta_2')
ylabel('V')

print('-f1','osc_bif_diagram','-djpeg')

%Zoom
if zoom == 1
    plot(caseIvec,yvec,'g-.','linewidth',2)
    plot(bifpredvec,yvec,'b-.','linewidth',2)
    plot(eta2,Vaxis,'k','linewidth',2)
    xlim([nsbif caseIboundary+2/Omega])
    ylim([-.1 .1])
    print('-f1','osc_cases','-djpeg');
    
    plot(bifpredvec,yvec,'b-.','linewidth',2)
    xlim([nsbif muosc+2/Omega])
    print('-f1','osc_bif_diagram_zoom','-djpeg');
end

% T plot
m = 5*n;
Vlow = linspace(-1,0,m);
Vmid = linspace(0,Vsmooth,m);
Vup = linspace(Vsmooth,1.5,m);


func=@(V)(eta1/(1+abs(V)));

Tlow = func(Vlow);
Tmid = func(Vmid);
Tup = func(Vup);

close(figure(2))
figure(2)
plot(Vlow,Tlow,'r','linewidth',2)
hold on
plot(Vmid,Tmid, 'k-.')
plot(Vup,Tup,'r','linewidth',2)
plot(vOscTimeSeries,tOscTimeSeries,'k--','linewidth',2)
xlabel('V')
ylabel('T')
xlim([Min Max])
ylim([1.7 eta1])

print('-f2','osc_bif_Tplot','-djpeg');

if zoom == 1
    xlim([-.1 .11])
    ylim([eta1-.2 eta1])
    print('-f2','osc_bif_Tplot_zoom','-djpeg');
end

%print('-f2','slow_bif_Tplot_zoom','-djpeg');

fprintf('The error is= %f\n',error)


if comp == 1
    %(2) Comparison plot
    invOmegavec = linspace(.001,.2,20);
    M = length(invOmegavec);
    bifpredvec = zeros(1,M);
    bifactualvec = zeros(1,M);

    parfor i = 1:M
        invOmega = invOmegavec(i);
        %Make sure to have a time step that will capture highly oscillatory
        %pieces
        h = min(.01,2*pi*invOmega/(5));
        start = nsbif-.5; stop = nsbif+1;
        %Have a parameter grid fine enough to see the different bifurcation
        %locations for highly oscillatory cases
        eta2 = fliplr(linspace(start,stop,600));
        N = length(eta2);
        Vfinal = zeros(1,M);
        Tfinal = zeros(1,M);

        %Numerical Solution
        for j = 1:N
            vDEosc = @(t,V,T)((eta1-eta2(j))-V*abs(V)-T+eta3*(T-V)+A*sin(t/invOmega));
            tDEosc = @(t,V,T)(eta1-T*(1+abs(V))+B*sin(t/invOmega));

            vInit = -eta2(j)/eta1+eta3;
            tInit = eta1/(1-vInit);

            %Solve equations for long term equilibrium, step size h
            [~,Vosc,Tosc] = RK2sys(vDEosc,tDEosc,vInit,tInit,0,50,h);

             Vfinal(j) = Vosc(end);
             Tfinal(j) = Tosc(end);
        end

        mosc = 2*eta1*(1-eta3)/(pi*abs(A))-pi*abs(A)*eta3^2/(4*eta1*(1-eta3));
        bifpredvec(i) = mosc*invOmega;
        bifactualvec(i) = eta2(find(Vfinal>criteria,1))-eta1*eta3;

    end

    close(figure(3))
    figure(3)
    plot(invOmegavec,bifpredvec,'k','linewidth',2);
    hold on
    plot(invOmegavec,bifactualvec,'r*');
    xlabel('\Omega^{-1}');ylabel('\eta_2-\eta_{2ns}')

     print('-f3','osc_Omegacomp','-djpeg')
end

% Calculates the V curve
function eta2 = Vcurve(V,eta1,eta3)
    T0 = eta1/(1+abs(V));
    eta2 = eta1+eta3*(T0-V)-T0-V*abs(V);
end