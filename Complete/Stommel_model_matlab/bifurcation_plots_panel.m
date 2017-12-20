%This script creates the array of plots on the stommel bifurcation

%The parameter values for each subplot
p1_eta1=4;p1_eta3=.375;
p2_eta1=4;p2_eta3=1;
p3_eta1=4;p3_eta3=1.875;


subplot(3,1,1)
eta1=p1_eta1;eta3=p1_eta3;
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h1=ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,2]);
set(h1,'linestyle','-.','color','k')

hold on
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h2=ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,0]);
set(h2,'color','r','linewidth',2)

hold on
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h3=ezplot(z,[eta1*eta3-1,eta1*eta3+1,.426651,2.5]);
set(h3,'color','r','linewidth',2)
xlabel('')

ylabel('V')
title('')
axis([.5 2.5 -.5 1.5])

subplot(3,1,2)
eta1=p2_eta1;eta3=p2_eta3;
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h1=ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,.5]);
set(h1,'color','r','linewidth',2)
title('')
xlabel('')

subplot(3,1,3)
eta1=p3_eta1;eta3=p3_eta3;
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
ezplot(z,[eta1*eta3-.5,eta1*eta3+.5,-1,.5])
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h1=ezplot(z,[eta1*eta3-1,eta1*eta3+1,-1,.5]);
set(h1,'linestyle','-.','color','k')
hold on
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h2=ezplot(z,[eta1*eta3-1,eta1*eta3+1,0,1]);
set(h2,'color','r','linewidth',2)
hold on
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h3=ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.226269,-1]);
set(h3,'color','r','linewidth',2)

title('')
axis([6.5 8.5 -1 .5])
