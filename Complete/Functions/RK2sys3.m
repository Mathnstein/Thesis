%RK2 for systems

function [t,y1,y2,y3] = RK2sys3(ode1,ode2,ode3,init1,init2,init3,L,R,h)
t=L:h:R;
y1=zeros(1,length(t));
y2=zeros(1,length(t));
y3=zeros(1,length(t));
y1(1)=init1;
y2(1)=init2;
y3(1)=init3;
for i=1:(length(t)-1)
y1(i+1)=y1(i)+h*ode1(t(i)+h/2,y1(i)+h/2*ode1(t(i),y1(i),y2(i),y3(i)),y2(i)+h/2*ode2(t(i),y1(i),y2(i),y3(i)),y3(i)+h/2*ode3(t(i),y1(i),y2(i),y3(i)));
y2(i+1)=y2(i)+h*ode2(t(i)+h/2,y1(i)+h/2*ode1(t(i),y1(i),y2(i),y3(i)),y2(i)+h/2*ode2(t(i),y1(i),y2(i),y3(i)),y3(i)+h/2*ode3(t(i),y1(i),y2(i),y3(i)));
y3(i+1)=y3(i)+h*ode3(t(i)+h/2,y1(i)+h/2*ode1(t(i),y1(i),y2(i),y3(i)),y2(i)+h/2*ode2(t(i),y1(i),y2(i),y3(i)),y3(i)+h/2*ode3(t(i),y1(i),y2(i),y3(i)));
end