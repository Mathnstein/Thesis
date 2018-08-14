% Plot comparison plots between NS and smooth
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\Conclusion')

% Slow

eps = linspace(.000001,.3,50);
slow_ns = .5*eps.*log(eps);
slow_smooth = eps.^(2/3)*(-2.3381);

figure(1)
plot(eps,slow_ns,'linewidth',2)
hold on
plot(eps,slow_smooth,'r--','linewidth',2)
set(gca,'fontsize',18)
xlabel('\epsilon','fontsize',32)
ylabel('x','fontsize',32)
xlim([0,.3])
ylim([-1,0])
print('-f1','slowcomp','-djpeg')


% Oscillatory

invOmega = linspace(.000001,.2,50);
A=2;
osc_ns = 4*abs(A)*invOmega/pi;
osc_smooth = A^2*invOmega.^2/2;

figure(2)
plot(invOmega,osc_ns,'linewidth',2)
hold on
plot(invOmega,osc_smooth,'r--','linewidth',2)
set(gca,'fontsize',18)
xlabel('\Omega^{-1}','fontsize',32)
ylabel('x','fontsize',32)
xlim([0,.2])
ylim([0,.5])
print('-f2','osccomp','-djpeg')

% Mixed

% Fix Omega
A=2;
invOmega = .1;
eps = linspace(.000001,.3,50);
mixed_ns = (2*abs(A)*invOmega/pi)^(1/3)*eps.^(2/3)*(-2.3381)+(4*abs(A)*invOmega/pi);
mixed_smooth = eps.^(2/3)*(-2.3381)+(A^2*invOmega^2/2);

figure(3)
plot(eps,mixed_ns,'linewidth',2)
hold on
plot(eps,mixed_smooth,'r--','linewidth',2)
set(gca,'fontsize',18)
xlabel('\epsilon','fontsize',32)
ylabel('x','fontsize',32)
xlim([0,.3])
ylim([-1,.5])
print('-f3','mixedcomp_eps','-djpeg')

% Fix eps
A=2;
invOmega = linspace(.000001,.2,50);
eps = .1;
mixed_ns = (2*abs(A)*invOmega/pi).^(1/3)*eps^(2/3)*(-2.3381)+(4*abs(A)*invOmega/pi);
mixed_smooth = eps^(2/3)*(-2.3381)+(A^2*invOmega.^2/2);

figure(4)
plot(invOmega,mixed_ns,'linewidth',2)
hold on
plot(invOmega,mixed_smooth,'r--','linewidth',2)
set(gca,'fontsize',18)
xlabel('\Omega^{-1}','fontsize',32)
ylabel('x','fontsize',32)
xlim([0,.2])
ylim([-.6,.3])
print('-f4','mixedcomp_omega','-djpeg')