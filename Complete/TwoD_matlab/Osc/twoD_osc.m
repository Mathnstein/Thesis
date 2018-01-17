%This script produces:
%i) A plot containing a prediction on tipping for the 2-D Stommel model    
% along with the ramped numerical solution.
%ii) A plot comparing the actual tipping solution vs the asymptotic
%solution for different choices of epsilon.

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\TwoD_Osc')

criteria = .2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Code Begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start on NS stable branch
%Works for n3<1
%Parameters
n=400;
eta1=4;eta3=.375;
A=1; B=1; Omega=30;
h=min(.01,2*pi/(10*Omega));

%Fixed parameter system
start=eta1*eta3-.5;stop=eta1*eta3+1;
eta2=fliplr(linspace(start,stop,n));
vOscTimeSeries=zeros(1,length(eta2));
tOscTimeSeries=vOscTimeSeries;

%Create equilibrium vector
for i=1:length(eta2)
    %Non-dim Stommel Equations    
    vDEosc=@(t,V,T)((eta1-eta2(i))-V*abs(V)-T+eta3*(T-V)+A*sin(Omega*t));
    tDEosc=@(t,V,T)(eta1-T*(1+abs(V))+B*sin(Omega*t));
    
    vInit=-eta2(i)/eta1+eta3;
    tInit=eta1/(1-vInit);
    
    
    %Solve equations for long term equilibrium, step size h
    [~,Vosc,Tosc] = RK2sys(vDEosc,tDEosc,vInit,tInit,0,100,h);
    
    %Keep steady state
    vOscTimeSeries(i) = Vosc(floor(100/h)-n+i);
    tOscTimeSeries(i) = Tosc(floor(100/h)-n+i);
end

%-------------------------------------------------
%Osc system bif values

mosc = 2*eta1*(1-eta3)/(pi*abs(A))-pi*abs(A)*eta3^2/(4*eta1*(1-eta3));
muosc = mosc/Omega+eta1*eta3;

bifactual = eta2(find(vOscTimeSeries>criteria,1));

Min = min(vOscTimeSeries);
Max = max(vOscTimeSeries);

yvec = linspace(Min,Max,300);
bifpredvec = muosc* ones(1,300);
bifactualvec = bifactual*ones(1,300);

error = abs(bifactual-muosc);


%--------------------Plot the dynamics plot--------------------------

figure(1)
%plot(eta2,vSteady,'r')
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h1=ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,1.5]);
set(h1,'linestyle','-.','color','k')
hold on
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h2=ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,0]);
set(h2,'color','r','linewidth',2)
hold on
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h3=ezplot(z,[eta1*eta3-1,eta1*eta3+1,.426651,1.5]);
set(h3,'color','r','linewidth',2)

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
plot(bifpredvec,yvec,'b')
plot(bifactualvec,yvec,'k--')
xlim([1.5 1.7])
ylim([-.15 .2])

print('-f1','osc_bif_diagram_zoom','-djpeg');

% T plot
V=linspace(-2,2,n);
T=zeros(length(V),1);
func=@(V)(eta1/(1+abs(V)));
for i=1:length(V)
    T(i)=func(V(i));
end
figure(2)
plot(V,T,'r','linewidth',2)
hold on
plot(vOscTimeSeries,tOscTimeSeries,'k','linewidth',2)
xlabel('V')
ylabel('T')
xlim([Min Max])
ylim([1.7 eta1])

print('-f2','osc_bif_Tplot','-djpeg');

fprintf('The error is= %f\n',error)
%{
%(2) Comparison plot
invOmegavec=linspace(.001,.3,20);
M=length(invOmegavec);
bifpredvec=zeros(1,M);
bifactualvec=zeros(1,M);
 
for i=1:M
    invOmega=invOmegavec(i);
    %Make sure to have a time step that will capture highly oscillatory
    %pieces
    h=min(.01,2*pi*invOmega/(5));
    start=eta1*eta3-.5;stop=eta1*eta3+1;
    %Have a parameter grid fine enough to see the different bifurcation
    %locations for highly oscillatory cases
    eta2=fliplr(linspace(start,stop,200));
    N=length(eta2);
    Vfinal=zeros(1,M);
    Tfinal=zeros(1,M);
    
     
    %Numerical Solution
    for j=1:N
        vDEosc=@(t,V,T)((eta1-eta2(j))-V*abs(V)-T+eta3*(T-V)+A*sin(t/invOmega));
        tDEosc=@(t,V,T)(eta1-T*(1+abs(V))+B*sin(t/invOmega));
    
        vInit=-eta2(j)/eta1+eta3;
        tInit=eta1/(1-vInit);
    
    
        %Solve equations for long term equilibrium, step size h
        [~,Vosc,Tosc] = RK2sys(vDEosc,tDEosc,vInit,tInit,0,100,h);
    
         Vfinal(j)=Vosc(end);
         Tfinal(j)=Tosc(end);
    end
 
    mosc = 2*eta1*(1-eta3)/(pi*abs(A))-pi*abs(A)*eta3^2/(4*eta1*(1-eta3));
    bifpredvec(i) = mosc*invOmega+eta1*eta3;
    bifactualvec(i) = eta2(find(Vfinal>criteria,1));
    
end
 
figure(2)
    plot(invOmegavec,bifpredvec,'k','linewidth',2);
    hold on
    plot(invOmegavec,bifactualvec,'r*');
    xlabel('\Omega^{-1}');ylabel('\eta_2-{\eta_2}_{ns}')
    
 print('-f2','osc_Omegacomp','-djpeg')

%}