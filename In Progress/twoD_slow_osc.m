% Full analysis on slow+osc twoD

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\Complete\Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Code Begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start on NS stable branch
%Works for n3<1
%In progress for making work for n3>1

%Parameters
n=100;
eta1=4;eta3=.375;
h=.01;
A=1;
Omega=100;

%Fixed parameter system
start=eta1*eta3-1;stop=eta1*eta3+1;
eta2=linspace(start,stop,n);
vSteady=zeros(1,length(eta2));
tSteady=vSteady;

%Create equilibrium vector
for i=1:length(eta2)
    %Non-dim Stommel Equations
    vDE=@(t,V,T)((eta1-eta2(i))-V*abs(V)-T+eta3*(T-V));
    tDE=@(t,V,T)(eta1-T*(1+abs(V)));
    if eta3<= 1
        vInit=-eta2(i)/eta1+eta3;
        tInit=eta1/(1-vInit);
    else
        vInit=(-eta2(i))/eta1+eta1;
        tInit=eta1/(1+vInit);
    end
    
    
    %Solve equations for long term equilibrium, step size h
    [~,V,T]=RK2sys(vDE,tDE,vInit,tInit,0,100,h);
    
    %Keep steady state
    vSteady(i)=V(end);
    tSteady(i)=T(end);
end
Min=min(vSteady);
Max=max(vSteady);

%-------------------------------------------------
%Ramped system 
eps=.01;
time=(stop-start)/(abs(eps));

%Non-dim Stommel Equations
vDE=@(t,V,T,eta2)((eta1-eta2)-V*abs(V)-T+eta3*(T-V)+A*sin(Omega * t));
tDE=@(t,V,T,eta2)(eta1-T*(1+abs(V)));

if eta3<=1
    eta2DE=@(t,V,T,eta2)(-eps);
else
    eta2DE=@(t,V,T,eta2)(eps);
end

if eta3<=1
    eta2Init=stop;
    vInit=-stop/eta1+eta3;
    tInit=eta1/(1-vInit);
else
    eta2Init=start;
    vInit=start/eta1+eta3;
    tInit=eta1/(1+vInit);
end

%Solve equations for numerical soln, step size h
[~,Vnum,Tnum,eta2num]=RK2sys3(vDE,tDE,eta2DE,vInit,tInit,eta2Init,0,time,h);


%Determine the tipping points for each case
if eta3<=1
    A=[eta3 1-eta3;eta1 1];
    [C,eigs]=eig(A);
    tip1=eta1*eta3-eps*log(eps)/eigs(1,1);
    tipx1=tip1*ones(1,n);
    tipy=linspace(Min,Max,n);

    realtip=eta2num(find(Vnum>-C(1,1),1));
    realtipvec=realtip*ones(1,n);
else 
    A=[-eta3 1-eta3;-eta1 -1];
    [C,eigs]=eig(A);
    tip1=eta1*eta3+eps*log(abs(eps))/eigs(1,1);
    tipx1=tip1*ones(1,n);
    tipy=linspace(Min,Max,n);

    realtip=eta2num(find(Vnum<C(1,1),1));
    realtipvec=realtip*ones(1,n);
end

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
plot(eta2num,Vnum,'k--')
%plot(realtipvec,tipy,'b')
%plot(tipx1,tipy,'b-.')
xlim([start,stop])
%title(['V Equilibria for \epsilon=' num2str(eps)])
xlabel('\eta_2')
ylabel('V')
%legend('Numerical w/ Fixed \eta_2','Ramped System','Tipping Criteria','Tipping Point: \lambda_1')

fprintf('Tipping point for %f is %f\nPredicted tip is %f\nError: %f\n',eps,realtip,tip1,abs(realtip-tip1))

%Plot the T-equilibria against the numerical solution
figure(2)
plot(eta2,tSteady,'r','linewidth',2)
hold on
plot(eta2num,Tnum,'k--','linewidth',2)
xlabel('\eta_2')
ylabel('T')