
% Set value to 1 for only lambda figure, 2 for epsilon, and 3 for both

value = 2;

%Hyperparameters
n=30;
criteria=.2;
A = 1;

if value == 1 || value == 3
    lower =.5;
    upper = 1.3;

    % This will compare lambda over a single choice in epsilon
    lambdavec = linspace(lower,upper,n);

    truetipvec = zeros(1,n);
    estimatedtipvec = zeros(1,n);
    slowtipvec = zeros(1,n);
    eps = .005;

    for i = 1:length(lambdavec)
    %Initial values and set up
    lambda = lambdavec(i);
    Omega = eps^(-lambda);

    mudel = eps * (-2.33811);
    muosc = 4 * abs(A)/(pi * Omega);
    B = eps^(lambda-1) * abs(A);

    %Numerics
    h=min(.01, 2 * pi/(10 * Omega));
    start = 1.5; stop = -.5;
    time = abs(stop - start)/eps;
    muInit = start;
    yInit = 1-sqrt(1 + muInit);

    yDE = @(t, y, mu)(-mu + 2 * abs(y)-y * abs(y) + A * sin(Omega * t));
    muDE = @(t, y, mu)(-eps);

    %Numerical Solution
    [~, y, mu] = RK2sys(yDE, muDE, yInit, muInit, 0, time, h);

    %Comparison
    truetipvec(i) = mu(find(y < criteria, 1, 'last'));
    estimatedtipvec(i) = muosc + mudel * (pi * B/2)^(1/3);
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
    plot(lambdavec, estimatedtipvec, 'k')
    plot(lambdavec, slowtipvec, 'b--')
    %plot(xvec,yvec,'r')
    xlabel('\lambda'); ylabel('\mu');
    xlim([lower upper])

    fileID = fopen('slowosc_lambdacomp_information.txt','w');
    fprintf(fileID,'epsilon=%f\n',lambda);
    fprintf(fileID,'A=%f\n',A);
    fprintf(fileID,'tipping criteria=%f\n',criteria);
    fclose(fileID);
    print('-f1','slowosc_lambdacomp','-djpeg')
end

if value == 2 || value == 3
    % This will compare epsilon over a single choice in lambda
    epsvec = linspace(.005,.01,n);
    truetipvec = zeros(1,n);
    estimatedtipvec = zeros(1,n);
    slowtipvec = zeros(1,n);

    %Set up for an example of each case
    lambdavec = [.5, 1.4];
    delete('slowosc_epscomp_information.txt')
    for k=2:3
        lambda = lambdavec(k-1);

        for i = 1:length(epsvec)
            %Initial values and set up
            eps = epsvec(i);
            Omega = eps^(-lambda);

            mudel = eps * (-2.33811);
            muosc = 4 * abs(A)/(pi * Omega);
            B = eps^(lambda-1) * abs(A);

            %Numerics
            h = min(.01, 2 * pi/(20 * Omega));
            start = 1.5; stop = -.5;
            time = abs(stop - start)/eps;
            muInit = start;
            yInit = 1 - sqrt(1 + muInit);
            yDE = @(t, y, mu)(-mu + 2 * abs(y)-y * abs(y) + A * sin(Omega * t));
            muDE = @(t, y, mu)(-eps);

            %Numerical Solution
            [~, y, mu] = RK2sys(yDE, muDE, yInit, muInit, 0, time, h);

            %Comparison
            truetipvec(i) = mu(find(y < criteria, 1, 'last'));
            estimatedtipvec(i) = muosc + mudel * (pi * B/2)^(1/3);
            slowtipvec(i) = eps * log(eps)/2;
        end

        fileID = fopen('slowosc_epscomp_information.txt','a');
        fprintf(fileID,'Case= %f\n',k);
        fprintf(fileID,'lambda=%f\n',lambda);
        fprintf(fileID,'A=%f\n',A);
        fprintf(fileID,'tipping criteria=%f\n \n',criteria);
        fclose(fileID);

        close(figure(k))
        figure(k)
        plot(epsvec, truetipvec, 'r*')
        hold on
        plot(epsvec, estimatedtipvec, 'k')
        plot(epsvec, slowtipvec, 'b--')
        xlabel('\epsilon'); ylabel('\mu')
    end
        print('-f2', 'slowosc_epscomp_case2', '-djpeg')
        print('-f3', 'slowosc_epscomp_case3', '-djpeg')
end
