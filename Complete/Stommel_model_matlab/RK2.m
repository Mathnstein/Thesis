%This function implements and graphs the time vector against the approx.
%function value according to an RK2 method (midpoint).

function [t,y] = RK2(ode,init,L,R,h)
t =L:h:R;
y =zeros(1,length(t));
y(1)=init;
for i=1:(length(t)-1)
y(i+1)=y(i)+h*ode(t(i)+h/2,y(i)+h/2*ode(t(i),y(i)));
end




