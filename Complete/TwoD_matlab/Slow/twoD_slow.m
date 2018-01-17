%This script produces:
%i) A plot containing a prediction on tipping for the 2-D Stommel model    
% along with the ramped numerical solution.
%ii) A plot comparing the actual tipping solution vs the asymptotic
%solution for different choices of epsilon.

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\TwoD_Slow')

criteria = .5;
% Previous criteria was -C(1,1), the eigvector component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Code Begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start on NS stable branch
% Works for n3<1

% Parameters
n = 400;
eta1 = 4; eta3 = .375;
h = .01;

%-------------------------------------------------
% Slow system 
eps = .01;
start = eta1*eta3-1;stop=eta1*eta3+1;
time = (stop-start)/(abs(eps));

% Non-dim Stommel Equations
vDE = @(t,V,T,eta2)((eta1-eta2)-V*abs(V)-T+eta3*(T-V));
tDE = @(t,V,T,eta2)(eta1-T*(1+abs(V)));
eta2DE = @(t,V,T,eta2)(-eps);

eta2Init = stop;
vInit = -stop/eta1+eta3;
tInit = eta1/(1-vInit);

% Solve equations for numerical soln, step size h
[~,Vnum,Tnum,eta2num] = RK2sys3(vDE,tDE,eta2DE,vInit,tInit,eta2Init,0,time,h);

% Ignore transient behavior
Vnum = Vnum(100:end);
Tnum = Tnum(100:end);
eta2num = eta2num(100:end);

% Determine the tipping points for each case
A=[eta3 1-eta3;eta1 1];
[C,eigs] = eig(A);
tipvec1 = eta1*eta3-eps*log(eps)/eigs(1,1);
tippredvec = tipvec1*ones(1,n);
yvec = linspace(Min,Max,n);

tipactual=eta2num(find(Vnum>criteria,1));
tipactualvec=tipactual*ones(1,n);

%--------------------Plot the dynamics plot--------------------------

figure(1)
z = @(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h1 = ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,1.5]);
set(h1,'linestyle','-.','color','k')
hold on
z = @(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h2 = ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,0]);
set(h2,'color','r','linewidth',2)
hold on
z = @(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h3 = ezplot(z,[eta1*eta3-1,eta1*eta3+1,.426651,1.5]);
set(h3,'color','r','linewidth',2)

xlabel('\eta_2')
ylabel('V')
title('')
hold on
plot(eta2num,Vnum,'k--','linewidth',2)
plot(tippredvec,yvec,'b')
plot(tipactualvec,yvec,'k--')
xlim([start,stop])
xlabel('\eta_2')
ylabel('V')

%Zoom

plot(tippredvec,yvec,'b')
plot(tipactualvec,yvec,'k--')
xlim([1.44 1.6])
ylim([-.05 .2])

print('-f1','slow_bif_zoom','-djpeg');


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
plot(Vnum,Tnum,'k','linewidth',2)
xlabel('V')
ylabel('T')
xlim([-.75 Max])
ylim([1.5 eta1])


print('-f2','slow_bif_Tplot','-djpeg');

fprintf('Tipping point for %f is %f\nPredicted tip is %f\nError: %f\n',eps,tipactual,tipvec1,abs(tipactual-tipvec1))


%%%%%%%%%%%%%%Comparison plot for different epsilon%%%%%%%%%%%%%%%%%%%%%%%
%3-system for tipping
epsvec = linspace(.001,.2,20);
M = length(epsvec);
tipactualvec = zeros(1,M);
tipvec1 = zeros(1,M);
tipvec2 = zeros(1,M);

for j = 1:M
    eps = epsvec(j);
    time = (stop-start)/(abs(eps));

    % Non-dim Stommel Equations
    vDE = @(t,V,T,eta2)((eta1-eta2)-V*abs(V)-T+eta3*(T-V));
    tDE = @(t,V,T,eta2)(eta1-T*(1+abs(V)));
    eta2DE = @(t,V,T,eta2)(-eps);

    eta2Init = stop;
    vInit = -stop/eta1+eta3;
    tInit = eta1/(1-vInit);
    
    % Solve equations for numerical soln, step size .01
    [~,Vnum,Tnum,eta2num] = RK2sys3(vDE,tDE,eta2DE,vInit,tInit,eta2Init,0,time,h);
    
     tipactual = eta2num(find(Vnum>criteria,1));
     tipactualvec(j) = tipactual-eta1*eta3;
     % Asymptotic Function
     asympfunc1 = @(epsilon)(-epsilon*log(epsilon)/eigs(1,1));
     asympfunc2 = @(epsilon)(-epsilon*log(epsilon)/eigs(2,2));
     % Form Asymptotic line
     tipvec1(j) = asympfunc1(epsvec(j));
     tipvec2(j) = asympfunc2(epsvec(j));
end

%Plot the Asymptotic solution vs TIpping point
figure(3)
plot(epsvec,tipvec1,'k','linewidth',2)
hold on
plot(epsvec,tipactualvec,'r*')
hold off
xlabel('\epsilon')
ylabel('\eta_2-\eta_{2ns}')

print('-f3','slow_epscomp','-djpeg')
