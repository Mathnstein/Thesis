% This script will create a graph of the non-smooth function we integrate
% over

addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\OneD_w_Osc')

n = 200;
A = 2;
v0 = 1.3;
time = linspace(0,2*pi,n);
horiz = v0*ones(1,n);

ylength = linspace(0,v0+A,n);
vertical1 = acos(v0/A)*ones(1,n);
vertical2 = (2*pi-acos(v0/A))*ones(1,n);

func = @(t)(abs(v0-A*cos(t)));
v = func(time);

close(figure(1))
figure(1)
plot(time,v,'r','LineWidth',2)
hold on
plot(time, horiz,'b--','LineWidth',2)
plot(vertical1,ylength,'k-.','LineWidth',2)
plot(vertical2,ylength,'k-.','LineWidth',2)
set(gca,'fontsize',18)
xlabel('T','fontsize',32)
ylabel('|v_0-Acos(T)|','fontsize',32)
xlim([0 2*pi])
ylim([0 v0+A])

print('-f1', 't1t2_graphic', '-djpeg')