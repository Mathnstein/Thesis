eta3=.5;
eta1=1/eta3-1;

figure(3)
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h1=ezplot(z,[eta1*eta3-2,eta1*eta3+2]);
set(h1,'linestyle','-.','color','k')

xlabel('\eta_2')
ylabel('V')
title('')

