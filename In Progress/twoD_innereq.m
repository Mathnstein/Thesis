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
lambda = .6;
Omega = eps^(-lambda);
A=1;
C=eps^(lambda-1)*A;
criteria = .5/eps;

%-------------------------------------------------
% System 
start = 10;stop=-5;
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

% Ignore transient behavior
%Vnum = Vnum(100:end);
%Tnum = Tnum(100:end);
%eta2num = eta2num(100:end);

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

% Compares the solutions x vs y
close(figure(2))
figure(2)

plot(xnum,ynum,'k--','linewidth',2)
xlim([-1 10])
xlabel('x')
ylabel('y')

close(figure(3))
figure(3)

plot(mnum,ynum,'k--','linewidth',2)
hold on
plot(tipactualvec,yvec,'k--')
ylim([-criteria 1])
xlabel('m')
ylabel('y')


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
tipactualequilvec=tipactualequil*ones(1,n);
yvec = linspace(-criteria,criteria,n);

fprintf('The equil. tipping point is m=%f\n',tipactualequil)

%--------------------Plot the dynamics plot--------------------------
% Compares the bifurcation plot in terms of x
figure(1)

plot(mnum,xnum,'r--','linewidth',2)
hold on
plot(tipactualequilvec,yvec,'r--')
ylim([-3 criteria])
xlabel('m')
ylabel('x')

figure(3)
plot(tipactualequilvec,yvec,'r--')
