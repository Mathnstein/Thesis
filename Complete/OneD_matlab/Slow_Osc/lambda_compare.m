%Full analysis of one-dimensional slow+osc

%Load function folder
addpath('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\Thesis-Master\trunk\Complete\Graphs\OneD_slowosc')

% Set value to 1 for only lambda figure, 2 for epsilon, and 3 for both

value = 1;

%Hyperparameters
n=25;
criteria=.2;
A = 1;

if value == 1 || value == 3
    lower =.5;
    upper = 2;

    % This will compare lambda over a single choice in epsilon
    lambdavec = linspace(lower,upper,n);

    truetipvec = zeros(1,n);
    estimatedtipvec = zeros(1,n);
    slowtipvec = zeros(1,n);
    eps = .01;

    parfor i = 1:n
    %Initial values and set up
    lambda = lambdavec(i);
    Omega = eps^(-lambda);

    mudel = eps * (-2.33811);
    muosc = 4 * abs(A)/(pi * Omega);
    B = eps^(lambda-1) * abs(A);

    %Numerics
    h=min(.01, 2 * pi/(10 * Omega));
    start = .5; stop = -.5;
    time = abs(stop - start)/eps;
    muInit = start;
    yInit = 1-sqrt(1 + muInit);

    yDE = @(t, y, mu)(-mu + 2 * abs(y)-y * abs(y) + A * sin(Omega * t));
    muDE = @(t, y, mu)(-eps);

    %Numerical Solution
    [~, y, mu] = RK2sys(yDE, muDE, yInit, muInit, 0, time, h);

    %Comparison
    truetipvec(i) = mu(find(y > criteria, 1));
    estimatedtipvec(i) = muosc + mudel * eps^((lambda-1)/3)*(pi * abs(A)/2)^(1/3);
    slowtipvec(i) = eps*log(eps)/2;
    end

    %Check to see minimum lambda holds
    miny = min(truetipvec)-eps;
    maxy = max(truetipvec);
    %yvec = linspace(miny,maxy,50);
    %xvec = minlambda*ones(1,50);

    close(figure(1))
    figure(1)
    plot(lambdavec, truetipvec, 'r*')
    hold on
    plot(lambdavec, estimatedtipvec, 'k','linewidth',2)
    plot(lambdavec, slowtipvec, 'b--','linewidth',2)
    %plot(xvec,yvec,'r')
    set(gca,'fontsize',14)
    xlabel('\lambda','fontsize',20);
    ylabel('\mu','fontsize',20);
    xlim([lower upper])

    print('-f1','slowosc_lambdacomp','-djpeg')
end

if value == 2 || value == 3
    % This will compare epsilon over a single choice in lambda
    epsvec = linspace(.005,.02,n);
    truetipvec = zeros(1,n);
    estimatedtipvec = zeros(1,n);
    slowtipvec = zeros(1,n);

    %Set up for an example of each case
    lambdavec = [.6, 2];
    for k=2:3
        lambda = lambdavec(k-1);

        parfor i = 1:n
            %Initial values and set up
            eps = epsvec(i);
            Omega = eps^(-lambda);

            mudel = eps * (-2.33811);
            muosc = 4 * abs(A)/(pi * Omega);
            B = eps^(lambda-1) * abs(A);

            %Numerics
            h = min(.01, 2 * pi/(10 * Omega));
            start = .5; stop = -.5;
            time = abs(stop - start)/eps;
            muInit = start;
            yInit = 1 - sqrt(1 + muInit);
            yDE = @(t, y, mu)(-mu + 2 * abs(y)-y * abs(y) + A * sin(Omega * t));
            muDE = @(t, y, mu)(-eps);

            %Numerical Solution
            [~, y, mu] = RK2sys(yDE, muDE, yInit, muInit, 0, time, h);

            %Comparison
            truetipvec(i) = mu(find(y > criteria,1));
            estimatedtipvec(i) = muosc + mudel * (pi * B/2)^(1/3);
            slowtipvec(i) = eps * log(eps)/2;
        end

        close(figure(k))
        figure(k)
        plot(epsvec, truetipvec, 'r*')
        hold on
        plot(epsvec, estimatedtipvec, 'k','linewidth',2)
        plot(epsvec, slowtipvec, 'b--','linewidth',2)
        set(gca,'fontsize',14)
        xlabel('\epsilon','fontsize',20)
        ylabel('\mu','fontsize',20)
    end
        print('-f2', 'slowosc_epscomp_case2', '-djpeg')
        print('-f3', 'slowosc_epscomp_case3', '-djpeg')
end
