% Phase plot
mulower=linspace(0,1.5,150);
lower=@(mu)(1-sqrt(1+mu));
xlower=lower(mulower);
 
muupper=linspace(-1,1,200);
upper=@(mu)(1+sqrt(1-mu));
xupper=upper(muupper);
 
mumid=linspace(0,1,100);
middle=@(mu)(1-sqrt(1-mu));
xmid=middle(mumid);
 

%(1) Dynamics plot
close(figure(2))
figure(2)

xlim([stop,start])
plot(mulower,xlower,'r','linewidth',2)
hold on
plot(mumid,xmid,'k:')
plot(muupper,xupper,'r','linewidth',2)
set(gca,'fontsize',18)
xlabel('\mu','fontsize',32)
ylabel('x','fontsize',32)
xlim([-1 1.5])
ylim([-1 2.5])