%This will plot a few vector fields along with the equilibrium solutions
%and a couple initial value solutions.

%Initial set up
eta2=.2;
eta1=3;eta3=.3;
ode1=@(t,s)(eta1-t*(1+abs(t-s)));
ode2=@(t,s)(eta2-s*(eta3+abs(t-s)));

figure(1)
%Algebraic roots of nonlinear ode
curvet=ezplot(ode1)
set(curvet,'Color','black')
hold on
curves=ezplot(ode2)
set(curves,'Color','red')
hold on
%Vector Field plot
[T,S]=meshgrid(0:0.2:5,0:0.2:5);
dt=eta1-T.*(1+abs(T-S));
ds=eta2-S.*(eta3+abs(T-S));
quiver(T,S,dt,ds);
axis([0 5 0 5])
hold on
%Solutions with chosen init values from graph
ode1=@(time,t,s)(eta1-t*(1+abs(t-s)));
ode2=@(time,t,s)(eta2-s*(eta3+abs(t-s)));
[time,solt,sols]=RK2sys(ode1,ode2,3,3,0,100,1000);
plot(solt,sols,'b')
title('eta2=.2')
hold off



eta2=1;
eta1=3;eta3=.3;
ode1=@(t,s)(eta1-t*(1+abs(t-s)));
ode2=@(t,s)(eta2-s*(eta3+abs(t-s)));
figure(2)
%Algebraic roots of nonlinear ode
curvet=ezplot(ode1)
set(curvet,'Color','black')
hold on
curves=ezplot(ode2)
set(curves,'Color','red')
hold on
%Vector Field plot
[T,S]=meshgrid(0:0.1:5,0:0.1:5);
dt=eta1-T.*(1+abs(T-S));
ds=eta2-S.*(eta3+abs(T-S));
quiver(T,S,dt,ds);
axis([0 5 0 5])
hold on
%Solutions with chosen init values from graph
ode1=@(time,t,s)(eta1-t*(1+abs(t-s)));
ode2=@(time,t,s)(eta2-s*(eta3+abs(t-s)));
[~,solt,sols]=RK2sys(ode1,ode2,1.25,3,3,100,1000);
plot(solt,sols,'b')
hold on
[~,solt,sols]=RK2sys(ode1,ode2,3,3.5,0,100,1000);
plot(solt,sols,'g')
title('eta2=1')


eta2=1.5;
eta1=3;eta3=.3;
ode1=@(t,s)(eta1-t*(1+abs(t-s)));
ode2=@(t,s)(eta2-s*(eta3+abs(t-s)));
figure(3)
%Algebraic roots of nonlinear ode
curvet=ezplot(ode1)
set(curvet,'Color','black')
hold on
curves=ezplot(ode2)
set(curves,'Color','red')
hold on
%Vector Field plot
[T,S]=meshgrid(0:0.2:5,0:0.2:5);
dt=eta1-T.*(1+abs(T-S));
ds=eta2-S.*(eta3+abs(T-S));
quiver(T,S,dt,ds);
axis([0 5 0 5])
hold on
%Solutions with chosen init values from graph
ode1=@(time,t,s)(eta1-t*(1+abs(t-s)));
ode2=@(time,t,s)(eta2-s*(eta3+abs(t-s)));
[time,solt,sols]=RK2sys(ode1,ode2,3,3,0,100,1000);
plot(solt,sols,'b')
title('eta2=1.5')
hold off