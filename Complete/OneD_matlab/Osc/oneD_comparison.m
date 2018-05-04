%Full analysis of OneD Osc.
%(1)A dynamics plot and a (2)comparison plot against frequency will be produced

criteria =.2;

% Load function folder
addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\OneD_w_Osc')
 

% Comparison plot
A=1;
M=20;
invOmegavec=linspace(.001,.2,M);
bifvec=zeros(1,length(invOmegavec));
bifactualvec=zeros(1,length(invOmegavec));
 
parfor i=1:M
    invOmega=invOmegavec(i);
    %Make sure to have a time step that will capture highly oscillatory
    %pieces
    h=min(.01,2*pi*invOmega/(10));
    start=.5;stop=0;
    %Have a parameter grid fine enough to see the different bifurcation
    %locations for highly oscillatory cases
    muper=fliplr(linspace(stop,start,500));
    N=length(muper);
    yfinal=zeros(1,N);
     
    %Numerical Solution
    for j=1:N
        yInit=1-sqrt(1+muper(j));
        yDE=@(t,y)(-muper(j)+2*abs(y)-y*abs(y)+A*sin(t/invOmega));
        [t,y]=RK2(yDE,yInit,0,15,h);
        yfinal(j)=y(end);
    end
 
    bifvec(i)=4*abs(A)*invOmega/(pi);
    bifactualvec(i)=muper(find(yfinal>criteria,1));
end

close(figure(1))
figure(1)
plot(invOmegavec,bifvec,'k','linewidth',2);
hold on
plot(invOmegavec,bifactualvec,'r*');
set(gca,'fontsize',14)
xlabel('\Omega^{-1}','fontsize',20);
ylabel('\mu','fontsize',20)


print('-f1','osc_Omegacomp','-djpeg')


