% Load function folder
addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\TwoD_Osc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Code Begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start on NS stable branch
%Works for n3<1
%Parameters
n = 600;
eta1 = 4; eta3 = 3/8;
nsbif = eta1*eta3;
A = 1; B = 1; Omega = 30;
h = min(.01, 2*pi/(10*Omega));

% Locate the smooth bifurcation
derivative = [-2 -(eta3+4) -2*(eta3+1) eta1-eta3*(eta1+1)];
r = roots(derivative);
Vsmooth = r(imag(r)==0);
smootheta2 = Vcurve(Vsmooth,eta1,eta3);

criteria = Vsmooth;

%Fixed parameter system
% Comparison plot
invOmegavec = linspace(.01,.2,20);
M = length(invOmegavec);
bifpredvec = zeros(1,M);
bifactualvec = zeros(1,M);
parfor i = 1:M 
    invOmega = invOmegavec(i);
    %Make sure to have a time step that will capture highly oscillatory
    %pieces
    h = min(.01,2*pi*invOmega/(5));
    start = nsbif; stop = nsbif+1;
    %Have a parameter grid fine enough to see the different bifurcation
    %locations for highly oscillatory cases
    eta2 = fliplr(linspace(start,stop,n));
    N = length(eta2);
    Vfinal = zeros(1,M);
    Tfinal = zeros(1,M);

    %Numerical Solution
    for j = 1:N
        vDEosc = @(t,V,T)((eta1-eta2(j))-V*abs(V)-T+eta3*(T-V)+A*sin(t/invOmega));
        tDEosc = @(t,V,T)(eta1-T*(1+abs(V))+B*sin(t/invOmega));
        
        vInit = -eta2(j)/eta1+eta3;
        tInit = eta1/(1-vInit);

        %Solve equations for long term equilibrium, step size h
        [~,Vosc,Tosc] = RK2sys(vDEosc,tDEosc,vInit,tInit,0,50,h);
        Vfinal(j) = Vosc(end);
        Tfinal(j) = Tosc(end);
    end

    mosc = 2*eta1*(1-eta3)/(pi*abs(A))-pi*abs(A)*eta3^2/(4*eta1*(1-eta3));
    bifpredvec(i) = mosc*invOmega;
    bifactualvec(i) = eta2(find(Vfinal>criteria,1))-nsbif;
end

close(figure(1))
figure(1)
plot(invOmegavec,bifpredvec,'k','linewidth',2);
hold on
plot(invOmegavec,bifactualvec,'r*');
set(gca,'fontsize',18)
xlabel('\Omega^{-1}','fontsize',32);
ylabel('\eta_2-\eta_{2ns}','fontsize',32)

print('-f1','osc_Omegacomp','-djpeg')

% Calculates the V curve
function eta2 = Vcurve(V,eta1,eta3)
    T0 = eta1/(1+abs(V));
    eta2 = eta1+eta3*(T0-V)-T0-V*abs(V);
end
