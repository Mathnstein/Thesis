addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\TwoD_Stommel')


%Standard orientation bif plot
eta3=0.375;
eta1=4;

figure(1)
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

xlabel('\eta_2')
ylabel('V')
title('')
axis([.5 2.5 -.5 1.5])
print('-f1','V_bif','-djpeg');

%Saddle orientation

eta1=4;
eta3=1;
figure(2)
z=@(eta2,V)(eta1-eta2)-V*abs(V)-eta1/(1+abs(V))+eta3*(eta1/(1+abs(V))-V);
h1=ezplot(z,[eta1*eta3-1,eta1*eta3+1,-.5,.5]);
set(h1,'color','r','linewidth',2)
xlabel('\eta_2')
ylabel('V')
title('')

%Reverse orientation

eta1=4;
eta3=1.875;
print('-f2','V_bif_reverse','-djpeg');

figure(3)
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

xlabel('\eta_2')
ylabel('V')
title('')
axis([6.5 8.5 -1 .5])

print('-f3','V_bif_collapse','-djpeg');

%T equilibrium
eta1 = 4; eta3 = 3/8;

% Locate the smooth bifurcation
derivative = [-2 -(eta3+4) -2*(eta3+1) eta1-eta3*(eta1+1)];
r = roots(derivative);
Vsmooth = r(real(r)>=0&imag(r)==0);
smootheta2 = Vcurve(Vsmooth,eta1,eta3);

% T plot
m=400;
Vlow=linspace(-2,0,m);
Vmid=linspace(0,Vsmooth,m);
Vup = linspace(Vsmooth,2,m);
Tlow=zeros(m,1);
Tmid = zeros(m,1);
Tup = zeros(m,1);

func=@(V)(eta1/(1+abs(V)));
for i=1:m
    Tlow(i) = func(Vlow(i));
    Tmid(i) = func(Vmid(i));
    Tup(i) = func(Vup(i));
end

close(figure(4))
figure(4)
plot(Vlow,Tlow,'r','linewidth',2)
hold on
plot(Vmid,Tmid, 'k-.')
plot(Vup,Tup,'r','linewidth',2)
xlabel('V')
ylabel('T')
print('-f4','T_equil','-djpeg');

function eta2 = Vcurve(V,eta1,eta3)
    T0 = eta1/(1+abs(V));
    eta2 = eta1+eta3*(T0-V)-T0-V*abs(V);
end