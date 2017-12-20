% This code will create a slow comparison plot for tip vs estimate

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\Complete\Functions')

%Start on NS stable branch
%Parameters
n=100;
eta1=4;eta3=.375;
start=eta1*eta3-1;stop=eta1*eta3+1;
eta2=linspace(start,stop,n);

A=[eta3 1-eta3;eta1 1];
[C eigs]=eig(A);

%Max slope
epsvec=linspace(.0001,.2,20);
maxslope=zeros(1,length(epsvec));
analysis=zeros(1,length(epsvec));
tip1=zeros(1,length(epsvec));
tip2=zeros(1,length(epsvec));

for j=1:length(epsvec)
eps=epsvec(j);
%If statement to determine if we go from left or right
    if eta3<=1
    else
        eps=-eps;
    end
%Non-dim Stommel Equations
vDE=@(t,V,T,eta2)((eta1-eta2)-V*abs(V)-T+eta3*(T-V));
tDE=@(t,V,T,eta2)(eta1-T*(1+abs(V)));
eta2DE=@(t,V,T,eta2)(-eps);
    
vInit=-stop/eta1+eta3;
tInit=eta1/(1-vInit);
    if eta3<=1
        eta2Init=stop;
    else
        eta2Init=start;
    end
time=(stop-start)/(abs(eps));
%Solve equations for numerical soln, step size .01
[~,Vnum,Tnum,eta2num]=RK2sys3(vDE,tDE,eta2DE,vInit,tInit,eta2Init,0,time,.01);

%Comparison
transient=round(length(Vnum)/3); %Ignore initial transient behavior
differences=zeros(1,length(Vnum)-transient);
    for i=transient:length(Vnum)-1
        differences(i+1-transient)=abs(Vnum(i+1)-Vnum(i));
    end
[val,pos]=max(differences);
%Collect each tipping point
maxslope(j)=eta2num(pos+transient)-eta1*eta3; %Add together the starting transient part of vector
realtip=eta2num(find(Vnum>-C(1,1),1)-1);
analysis(j)=realtip-eta1*eta3;
%Asymptotic Function
asympfunc1=@(epsilon)(-epsilon*log(epsilon)/eigs(1));
asympfunc2=@(epsilon)(-epsilon*log(epsilon)/eigs(2));
%Form Asymptotic line
tip1(j)=asympfunc1(epsvec(j));
tip2(j)=asympfunc2(epsvec(j));
end

figure(1)
plot(epsvec,maxslope,'k*')
hold on
plot(epsvec,analysis,'r*')
plot(epsvec,tip1,'g')

title('Tipping Criteria Comparison')
xlabel('\epsilon')
ylabel('\eta_2-\eta_{2ns}')
legend('Max Slope','Analysis C value','Asympt: \lambda_1')
