%Full analysis of one-dimensional slow+osc

%Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')

%Initial values and set up
eps = .01;
lambda = 2.5;
Omega = eps^(-lambda);
A = 1;
criteria=.2;

mudel = eps * (-2.33811);
muosc = 4 * abs(A)/(pi * Omega);
B = eps^(lambda-1) * abs(A);

%Numerics
h = min(.01, 2 * pi/(10 * Omega));
start = 1.5; stop = -1;
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
region2boundary = muosc + mudel * (pi * B/2)^(1/3);
slowtip = .5 * eps * log(eps);

%Inner Asymptotic solution

m = zeros(1, length(mu));

%Here I distinguish between the osc Dom and slow Dom
outer = zeros(1, length(mu));
case1 = zeros(1, length(mu));
    
if lambda < 1.5
    case2 = zeros(1, length(mu));
    for i = 1:length(m)
       %Contains the cases from exploding, chose an artificial lower
       %threshold
         if mu(i) > region2boundary+eps
               m(i) = mu(i);
         else
               m(i) = 0;
         end
       outer(i) = 1 - sqrt(1 + m(i)) - eps/(4 * (1 + m(i)));
       case1(i) = -m(i)/2;
       inside = (2/(pi * B))^(1/3) * (m(i)/eps - muosc/eps);
       case2(i) = eps*(pi * B/2)^(2/3) * airy(1, inside)/airy(inside);
    end
    abserror = abs(truetip - region2boundary);
    estimated_tipping = region2boundary;
else
    case3 = zeros(1, length(mu));
    for i = 1:length(m)
        %Contains the cases from exploding, chose an artificial lower
       %threshold
        if mu(i) > eps * log(eps)/2
             m(i) = mu(i);
        else
             m(i) = 0;
        end 
        outer(i) = 1 - sqrt(1 + m(i)) - eps/(4 * (1 + m(i)));
        case1(i) = -m(i)/2;
        case3(i) = eps * exp(-2 * m(i)/eps) + m(i)/2 - 1/4 * eps;
    end
    abserror=abs(truetip - slowtip);
    estimated_tipping = slowtip;
end

%Forms the correct inner by pairing the right cases together;
asymp = zeros(1, length(mu));

for i = 1:length(mu)
    %Builds outer region
    if mu(i) > 1
        asymp(i) = outer(i) - eps^lambda * A * cos(-eps^(-lambda - 1) * mu(i));
    else
        %Builds region 1
        if mu(i) > region1boundary
            asymp(i) = case1(i) - eps^lambda * A * cos(-eps^(-lambda - 1) * mu(i));
        else
            %Builds region 2
            if lambda < 1.5
                asymp(i) = case2(i);
            else
                asymp(i) = case3(i);
            end
        end
    end
end

%Create region dividers
region1x = region1boundary * ones(1, 1000);
region2x = region2boundary * ones(1, 1000);
tipping = estimated_tipping * ones(1,1000);
truetipvec = truetip * ones(1, 1000);
yvec = linspace(min(y), max(y), 1000);


%Numerics vs Asymptotics
figure(3)
plot(mu, y, 'k--')
hold on
plot(mu, asymp, 'r')
plot(region1x, yvec, 'g')
plot(region2x, yvec, 'b')
plot(tipping, yvec, 'r-.')
plot(truetipvec, yvec, 'k--')

pos = [.6 .2 .5 .5];
str = {
   ['\lambda=' num2str(lambda)], ['\epsilon= ' num2str(eps)],...
   ['\Omega= ' num2str(Omega)], ['A= ' num2str(A)],...
   ['Abs. Err.= ' num2str(abserror)], ['Tip. Crit > ' num2str(criteria)]};
annotation('textbox', pos, 'String', str, 'FitBoxToText', 'on');


%The bifurcation diagram
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
xlabel('\mu'); ylabel('x')
title('')

cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\OneD_slowosc')
print('-f1','slowosc_bif_diagram','-djpeg')


%zoom
plot(region1x, yvec, 'g')
plot(region2x, yvec, 'b')
plot(truetipvec, yvec, 'k--')
xlim([-.1 .5])
ylim([-.3 .5])
print('-f1', 'slowosc_bif_diagram_zoom','-djpeg')

%Print a report of frequency and tipping comparison values
fileID = fopen('slowosc_bif_information.txt','w');
fprintf(fileID,'The actual tipping: %f\n Estimate: %f\n', truetip, region2boundary);
fprintf(fileID,'The slow tip: %f\n', slowtip);
fprintf(fileID,'The osc tip: %f\n', muosc);
fprintf(fileID,'Absolute error: %f\n', abserror);
fprintf(fileID,'Tipping Criteria: >%f\n', criteria);
fprintf(fileID,'epsilon=%f\n',eps);
fprintf(fileID,'lambda=%f\n',lambda);
fprintf(fileID,'A=%f\n',A);
fprintf(fileID,'Omega=%f',Omega);
fclose(fileID);

