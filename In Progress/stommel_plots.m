%This script creates the array of plots on the stommel bifurcation

%The parameter values for each subplot
p1_eta1=4;p1_eta3=.375;
p2_eta1=4;p2_eta3=1;
p3_eta1=4;p3_eta3=1.875;


subplot(3,1,1)
eta1=p1_eta1;eta3=p1_eta3;
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,1.5])
title('')
xlabel('')

subplot(3,1,2)
eta1=p2_eta1;eta3=p2_eta3;
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,.5])
title('')
xlabel('')

subplot(3,1,3)
eta1=p3_eta1;eta3=p3_eta3;
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
ezplot(z,[eta1*eta3-.5,eta1*eta3+.5,-1,.5])
title('')
