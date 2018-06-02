addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\Intro')

a = linspace(0,1,500);

lower = -sqrt(a);
upper = sqrt(a);

close(figure(1))
figure(1)
plot(a,lower,'k-.')
hold on
plot(a,upper,'r','linewidth',2)
set(gca,'fontsize',18)
xlabel('a','fontsize',32)
ylabel('x','fontsize',32)
xlim([-.5 1])

[a,x]=meshgrid(-1:0.2:1,-1:0.2:1);
dx=a-x.^2;
da=zeros(size(dx));
quiver(a,x,da,dx);
axis([-.5 1 -1 1])

print('-f1','saddlenode_bif_diagram','-djpeg')

%Initial values and set up
h=.01;
eps=.01;
start=1;stop=-.5;
time=abs(stop-start)/eps;
aInit=start;
xInit = sqrt(aInit);

xDE=@(t,x,a)(a-x^2);
aDE=@(t,x,a)(-eps);

%Numerical Solution
[~,xnum,anum]=RK2sys(xDE,aDE,xInit,aInit,0,time,h);

figure(1)
plot(anum,xnum,'k--','linewidth',2)
print('-f1','saddlenode_tipping','-djpeg')

% Zoom
xlim([-.1 .2])
ylim([-.2 .5])
print('-f1','saddlenode_tipping_zoom','-djpeg')

constvec = [1 0 -1];
for i = 1:3
    const= constvec(i);
    alimit = sqrt(abs(const));
    xlimit = abs(const);
    [a,x]=meshgrid(-alim-1:0.2:alim+1,-xlimit-1:0.2:xlimit+1);
    a=const*ones(size(a));
    dx=a-x.^2;
    da=zeros(size(dx));

    subplot(3,1,i)
    quiver(x,a,dx,da);
    axis([-xlimit-1 xlimit+1 alimit-.4 alimit+.1])
    axis off

    if const >0
        hold on
        plot(sqrt(const),const,'o',...
        'LineWidth',2,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[1 0 0],...
        'MarkerSize',10)
        plot(-sqrt(const),const,'ko','linewidth',2,'MarkerSize',10)
        title('a>0','fontsize',20)
        txt1 = '\surd{a}';
        text(sqrt(const)-.1,const-.2,txt1,'fontsize',18)
        txt2 = '-\surd{a}';
        text(-sqrt(const)-.1,const-.2,txt2,'fontsize',18)
    elseif const == 0
        hold on
        plot(0,0,'ko','linewidth',2,'MarkerSize',10)
        title('a=a_{bif}','fontsize',20)
        txt1 = '0';
        text(-.03,-.2,txt1,'fontsize',18)
    else
        title('a<0','fontsize',20)
        axis([-xlimit-1 xlimit+1 -alimit-.4 -alimit+.1])
    end
end
print('-f1','saddlenode','-djpeg')

