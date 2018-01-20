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

Min = min(Vnum);
Max = max(Vnum);

% Determine the tipping points for each case
A=[eta3 1-eta3;eta1 1];
[C,eigs] = eig(A);
eig1 = eigs(1,1);
eig2 = eigs(2,2);
tippred1 = eta1*eta3-eps*log(eps)/eig1;
tippred1vec = tippred1*ones(1,n);
yvec = linspace(Min,Max,n);

tipactual=eta2num(find(Vnum>criteria,1));
tipactualvec=tipactual*ones(1,n);

%--------------------Plot the dynamics plot--------------------------

close(figure(1))
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
xlim([start,stop])
xlabel('\eta_2')
ylabel('V')

print('-f1','slow_bif_diagram','-djpeg');


%Zoom

plot(tippred1vec,yvec,'b')
plot(tipactualvec,yvec,'k--')
xlim([1.44 1.6])
ylim([-.05 .2])

print('-f1','slow_bif_diagram_zoom','-djpeg');


% Locate the smooth bifurcation
derivative = [-2 -(eta3+4) -2*(eta3+1) eta1-eta3*(eta1+1)];
r = roots(derivative);
Vsmooth = r(real(r)>=0&imag(r)==0);
smootheta2 = Vcurve(Vsmooth,eta1,eta3);

% T plot
m=5*n;
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
xlabel('V')
ylabel('T')
xlim([-.75 Max])
ylim([1.5 eta1])

print('-f2','slow_bif_Tplot','-djpeg');

xlim([-.05 .1])
ylim([3.8 4.01])

print('-f2','slow_bif_Tplot_zoom','-djpeg');

error = abs(tipactual-tippred1);

% Print a report of all values
fileID = fopen('slow_bif_information.txt','w');
fprintf(fileID,'The actual tipping: %f\n Estimate: %f\n', tipactual, tippred1);
fprintf(fileID,'Absolute error: %f\n', error);
fprintf(fileID,'Tipping Criteria: >%f\n', criteria);
fprintf(fileID,'epsilon=%f\n',eps);
fprintf(fileID, 'lambda_1 used = %f, lambda_2 not used = %f\n',eig1,eig2);
fclose(fileID);

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

close(figure(3))
%Plot the Asymptotic solution vs TIpping point
figure(3)
plot(epsvec,tipvec1,'k','linewidth',2)
hold on
plot(epsvec,tipactualvec,'r*')
hold off
xlabel('\epsilon')
ylabel('\eta_2-\eta_{2ns}')

print('-f3','slow_epscomp','-djpeg')

% Calculates the V curve
function eta2 = Vcurve(V,eta1,eta3)
    T0 = eta1/(1+abs(V));
    eta2 = eta1+eta3*(T0-V)-T0-V*abs(V);
end