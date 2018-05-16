%Full analysis of one-dimensional slow+osc

%Load function folder
addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\OneD_slowosc')

%Initial values and set up
eps = .01;
lambda = 1.05;
Omega = eps^(-lambda);
A = 4;
criteria=.5;

mudel = eps * (-2.33811);
muosc = 4 * abs(A)/(pi * Omega);
B = eps^(lambda-1) * abs(A);

%Numerics
h = min(.01, 2 * pi/(10 * Omega));
start = 1; stop = -.5;
time = abs(stop - start)/eps;
muInit = start;
yInit = 1-sqrt(1 + muInit);

yDE = @(t, y, mu)(-mu + 2 * abs(y)-y * abs(y) + A * sin(Omega*t));
muDE = @(t, y, mu)(-eps);

%Numerical Solution
[~, y, mu] = RK2sys(yDE, muDE, yInit, muInit, 0, time, h);

%I use the solution passing epsilon as tipping here.
%Comparison
truetip = mu(find(y < criteria, 1, 'last'));
region1boundary = 2 * abs(A)/Omega;
mixedtipping = muosc + mudel * (pi * B/2)^(1/3);
slowtipping = .5 * eps * log(eps);


%Create region dividers
region1x = region1boundary * ones(1, 1000);
mixedtippingvec = mixedtipping * ones(1,1000);
slowtippingvec = slowtipping * ones(1,1000);
truetipvec = truetip * ones(1, 1000);
yvec = linspace(min(y), max(y), 1000);



%The bifurcation diagram
close(figure(1))
figure(1)

%plots the underlying bifurcation diagram
syms mu1 x
z = @(mu, y)(-mu+2*abs(y)-y*abs(y));
h1 = ezplot(z, [-1, 1.5, -.6, 2.5]);
set(h1, 'linestyle', '-.', 'color', 'k')
hold on
h2 = ezplot(z, [-1, 1.5, -.6, 0]);
set(h2, 'color', 'r', 'linewidth', 2)
h3 = ezplot(z, [-1, 1.5, 1, 2.5]);
set(h3, 'color', 'r', 'linewidth', 2)
plot(mu, y, 'k--', 'linewidth', 1.2)

%Label the graph
set(gca,'fontsize',18)
xlabel('\mu','fontsize',32); ylabel('x','fontsize',32)
title('')

if lambda> 1.5
    print('-f1','slowosc_bif_diagram_large','-djpeg')
elseif lambda>1 && lambda<1.5
    print('-f1', 'slowosc_bif_diagram_medium','-djpeg')
else
    print('-f1','slowosc_bif_diagram_small','-djpeg')
end


%zoom
plot(mixedtippingvec, yvec, 'b--','linewidth',2)
plot(truetipvec, yvec, 'k','linewidth',2)
xlim([mixedtipping-.05 mixedtipping+.1])
ylim([-.1 .3])
if lambda>1.5
    plot(slowtippingvec, yvec, 'g--','linewidth',2)
    print('-f1', 'slowosc_bif_diagram_large_zoom','-djpeg')
elseif lambda>1
    plot(slowtippingvec, yvec, 'g--','linewidth',2)
    print('-f1', 'slowosc_bif_diagram_medium_zoom','-djpeg')
else
    print('-f1', 'slowosc_bif_diagram_small_zoom','-djpeg')
end
