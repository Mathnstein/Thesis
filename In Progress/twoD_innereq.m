% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Code Begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start on NS stable branch
% Works for n3<1

% Hyperparameters
n = 500;
eta1 = 4; eta3 = .375;
h = .01;
eps = .01;
lambda = .7;
Omega = eps^(-lambda);
A=1;
B=1;
C=eps^(lambda-1)*A;
criteria = .5/eps;

%-------------------------------------------------
% Inner System 
start = 10; stop=-10;
time = (start-stop);

quadcoef = (eta1/(pi*abs(C)));
constant = 2*eta1*abs(C)/pi;

% Non-dim Stommel Equations
xDE = @(t,x,y,m)(-eta3*x-(1-eta3)*y-m);
yDE = @(t,x,y,m)(-y-quadcoef*x^2-constant);
mDE = @(t,x,y,m)(-1);

a = (1-eta3)*quadcoef;
c = (1-eta3)*constant;

mInit = start;
xInit = (eta3-sqrt(eta3^2+4*a*(mInit-c)))/(2*a);
yInit = -quadcoef*xInit^2-constant;

% Solve equations for numerical soln, step size h
[~,xnum,ynum,mnum] = RK2sys3(xDE,yDE,mDE,xInit,yInit,mInit,0,time,h);



tipactual=mnum(find(xnum>criteria,1));
tipactualvec=tipactual*ones(1,n);
yvec = linspace(-criteria,criteria,n);
fprintf('The tipping point is m=%f\n',tipactual)

%--------------------Plot the dynamics plot--------------------------
% Compares the bifurcation plot in terms of x
close(figure(1))
figure(1)

plot(mnum,xnum,'k--','linewidth',2)
hold on
plot(tipactualvec,yvec,'k--')
ylim([-3 criteria])
xlabel('m')
ylabel('x')

% Here we recreate the outer solution
tvec = linspace(0,time,length(xnum));
eta2num_inner = eps*mnum+eta1*eta3;
Vnum_inner = eps*xnum-eps^(lambda)*A*cos(Omega*tvec);
Tnum_inner = eta1+eps*ynum-eps^(lambda)*B*cos(Omega*tvec);

% Code from twoD_slowosc
start = eta1*eta3-1; stop=eta1*eta3+1;
time = (stop-start)/(abs(eps));

% Non-dim Stommel Equations
vDE = @(t,V,T,eta2)((eta1-eta2)-V*abs(V)-T+eta3*(T-V)+A*sin(Omega*t));
tDE = @(t,V,T,eta2)(eta1-T*(1+abs(V))+B*sin(Omega*t));
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

% End code from twoD slowosc


% Compares the solutions V vs T
close(figure(2))
figure(2)
plot(Vnum,Tnum,'k','linewidth',2)
hold on
plot(Vnum_inner,Tnum_inner,'r','linewidth',2)
xlim([-.75 1.5])
ylim([1.5 4])
xlabel('V')
ylabel('T')

% Compares eta2 vs V

close(figure(3))
figure(3)
plot(eta2num_inner,Vnum_inner,'r','linewidth',2)
hold on
plot(eta2num,Vnum,'k','linewidth',2)
%plot(tipactualvec,yvec,'k--')
ylim([-.2 1])
xlim([1.4 1.6])
xlabel('\eta_2')
ylabel('V')

%{
% Lets try to solve this problem assuming y is in equilibrium

% System assuming y stays in equil

quadcoef = (eta1/(pi*abs(C)));
constant = 2*eta1*abs(C)/pi;

% Non-dim Stommel Equations
xDE = @(t,x,m)(-m-eta3*x-(1-eta3)*(-quadcoef*x^2-constant));
mDE = @(t,x,m)(-1);

a = (1-eta3)*quadcoef;
c = (1-eta3)*constant;

mInit = start;
xInit = (eta3-sqrt(eta3^2+4*a*(mInit-c)))/(2*a);

% Solve equations for numerical soln, step size h
[~,xnum,mnum] = RK2sys(xDE,mDE,xInit,mInit,0,time,h);

tipactualequil=mnum(find(xnum>criteria,1));
%tipactualequilvec=tipactualequil*ones(1,n);
yvec = linspace(-criteria,criteria,n);

fprintf('The equil. tipping point is m=%f\n',tipactualequil)

%--------------------Plot the dynamics plot--------------------------
% Compares the bifurcation plot in terms of x
figure(1)

plot(mnum,xnum,'r--','linewidth',2)
hold on
%plot(tipactualequilvec,yvec,'r--')

xlabel('m')
ylabel('x')

%}

