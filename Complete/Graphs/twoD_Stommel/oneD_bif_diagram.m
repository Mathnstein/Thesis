%One dimensional bifurcation diagram.

mulower=linspace(0,1.5,150);
lower=@(mu)(1-sqrt(1+mu));
xlower=lower(mulower);

muupper=linspace(-1,1,200);
upper=@(mu)(1+sqrt(1-mu));
xupper=upper(muupper);

mumid=linspace(0,1,100);
middle=@(mu)(1-sqrt(1-mu));
xmid=middle(mumid);

figure(1)
plot(mulower,xlower,'r','linewidth',2)
hold on
plot(mumid,xmid,'k:')
plot(muupper,xupper,'r','linewidth',2)
plot(0,0,'ko','markersize',10)
plot(1,1,'kx','markersize',10)
xlabel('\mu')
ylabel('x')
legend('x_l','x_m','x_u','Non-smooth Bifurcation','Smooth Bifurcation')
