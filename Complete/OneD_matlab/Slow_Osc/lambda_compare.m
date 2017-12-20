% Create an epsilon comparison plot an lambda comparison plot

% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')

%Hyperparameter
criteria=.2;

% This will compare lambda over a single choice in epsilon
n=20;

lambdavec = linspace(.5,1.5,n);

truetipvec = zeros(1,n);
estimatedtipvec = zeros(1,n);
slowtipvec = zeros(1,n);

for i = 1:length(lambdavec)
%Initial values and set up
eps = .01;
lambda = lambdavec(i);
Omega = eps^(-lambda);
A = 1;

mudel = eps * (-2.33811);
muosc = 4 * abs(A)/(pi * Omega);
B = eps^(lambda-1) * abs(A);

%Numerics
h = min(.01, 2 * pi/(10 * Omega));
start = 1.5; stop = -.5;
time = abs(stop - start)/eps;
muInit = start;
yInit = 1-sqrt(1 + muInit);

yDE = @(t, y, mu)(-mu + 2 * abs(y)-y * abs(y) + A * sin(Omega*t));
muDE = @(t, y, mu)(-eps);

%Numerical Solution
[~, y, mu] = RK2sys(yDE, muDE, yInit, muInit, 0, time, h);

%I use the solution passing epsilon as tipping here.
%Comparison
truetipvec(i) = mu(find(y < criteria, 1, 'last'));
estimatedtipvec(i) = muosc + mudel * (pi * B/2)^(1/3);
slowtipvec(i) = eps*log(eps)/2;
end

figure(1)
plot(lambdavec, truetipvec, 'r*')
hold on
plot(lambdavec, estimatedtipvec, 'k')
plot(lambdavec, slowtipvec, 'b')
xlabel('\lambda'); ylabel('\mu')

pos = [.6 .2 .5 .5];
str = {['\epsilon= ' num2str(eps)],...
   ['A= ' num2str(A)], ['Tip. Crit > ' num2str(criteria)]};
annotation('textbox', pos, 'String', str, 'FitBoxToText', 'on');

cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\OneD_slowosc')
print('-f1','oneD_slowosc_lambdacomp','-dpdf')

% This will compare epsilon over a single choice in lambda
n=20;

epsvec = linspace(.005,.01,n);

truetipvec = zeros(1,n);
estimatedtipvec = zeros(1,n);
slowtipvec = zeros(1,n);

for i = 1:length(epsvec)
%Initial values and set up
eps = epsvec(i);
lambda = .8;
Omega = eps^(-lambda);
A = 1;

mudel = eps * (-2.33811);
muosc = 4 * abs(A)/(pi * Omega);
B = eps^(lambda-1) * abs(A);

%Numerics
h = min(.01, 2 * pi/(10 * Omega));
start = 1.5; stop = -.5;
time = abs(stop - start)/eps;
muInit = start;
yInit = 1 - sqrt(1 + muInit);

yDE = @(t, y, mu)(-mu + 2 * abs(y)-y * abs(y) + A * sin(Omega*t));
muDE = @(t, y, mu)(-eps);

%Numerical Solution
[~, y, mu] = RK2sys(yDE, muDE, yInit, muInit, 0, time, h);

%I use the solution passing epsilon as tipping here.
%Comparison
truetipvec(i) = mu(find(y < criteria, 1, 'last'));
estimatedtipvec(i) = muosc + mudel * (pi * B/2)^(1/3);
slowtipvec(i) = eps * log(eps)/2;
end

figure(2)
plot(epsvec, truetipvec, 'r*')
hold on
plot(epsvec, estimatedtipvec, 'k')
plot(epsvec, slowtipvec, 'b')
xlabel('\epsilon'); ylabel('\mu')

pos = [.6 .2 .5 .5];
str = {['\lambda= ' num2str(lambda)],...
   ['A= ' num2str(A)], ['Tip. Crit > ' num2str(criteria)]};
annotation('textbox', pos, 'String', str, 'FitBoxToText', 'on');

cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\OneD_slowosc')
print('-f2','oneD_slowosc_epscomp','-dpdf')

