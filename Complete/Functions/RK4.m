%This is code to implement RK4

function [t,y] = RK4(ode,init,L,R,h)
t =L:h:R;
y =zeros(1,length(t));
y(1)=init;
for i=1:(length(t)-1)                              
    k_1 = ode(t(i),y(i));
    k_2 = ode(t(i)+0.5*h,y(i)+0.5*h*k_1);
    k_3 = ode((t(i)+0.5*h),(y(i)+0.5*h*k_2));
    k_4 = ode((t(i)+h),(y(i)+k_3*h));

    y(i+1) = y(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h; 
end