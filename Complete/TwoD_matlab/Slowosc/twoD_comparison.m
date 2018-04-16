% Load function folder
addpath('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Functions')
cd('C:\Users\codyg\Desktop\MSc_Thesis\Cody\trunk\Complete\Graphs\TwoD_SlowOsc')

criteria = .2;

% Set value to 1 for only lambda figure, 2 for epsilon, and 3 for both
value = 3;

if value == 1 || value == 3
    %Here we will fix lambda and compare across epsilon
    k = 20;
    epsvec = linspace(.005,.01,k);
    tippredveclarge = zeros(1,k);
    tippredREDveclarge = zeros(1,k);
    tippredREDvecsmall = zeros(1,k);
    tipactualvec = zeros(1,k);
    
    lambdavec = [.6, 1.4];
    for l=2:3
        lambda = lambdavec(l-1);

        parfor i =1:k
            eps = epsvec(i);
            Omega = eps^(-lambda);
            A=2;
            B=2;
            h = min(.01,2*pi/(10*Omega));

            %-------------------------------------------------
            % System 
            start = eta1*eta3-1; stop=eta1*eta3+1;
            time = (stop-start)/(eps);

            % Non-dim Stommel Equations
            vDE = @(t,V,T,eta2)((eta1-eta2)-V*abs(V)-T+eta3*(T-V)+A*sin(Omega*t));
            tDE = @(t,V,T,eta2)(eta1-T*(1+abs(V))+B*sin(Omega*t));
            eta2DE = @(t,V,T,eta2)(-eps);

            eta2Init = stop;
            vInit = -stop/eta1+eta3;
            tInit = eta1/(1-vInit);

            % Solve equations for numerical soln, step size h
            [~,Vnum,Tnum,eta2num] = RK2sys3(vDE,tDE,eta2DE,vInit,tInit,eta2Init,0,time,h);

            % Ignore transient behavior
            Vnum = Vnum(100:end);
            Tnum = Tnum(100:end);
            eta2num = eta2num(100:end);

            Min = min(Vnum);
            Max = max(Vnum);

            % Determine the tipping points for each case
            Mat = [eta3 1-eta3;eta1 1];
            tippredlarge = eta1*eta3-eps*log(eps)/min(eig(Mat));
            tippredREDlarge = eta1*eta3+eps*log(eps)/(eta1-eta3-eta1*eta3);
            
            scaling = eps^(lambda-1);
            c0 = pi*abs(A)/(eta1*(1-eta3));
            c1 = eta1*(1-eta3)*abs(A)/pi;
            c2 = eta3*pi/(2*eta1*(1-eta3));
            smooth = scaling^(1/3)*c0^(1/3)*(-eps*2.33811);
            oscillatory = c1*(2-c2^2)/Omega;
            
            tippredREDsmall = eta1*eta3+smooth+oscillatory
            tipactual=eta2num(find(Vnum>criteria,1));

            tippredveclarge(i) = tippredlarge;
            tippredREDveclarge(i) = tippredREDlarge;
            tippredREDvecsmall(i) = tippredREDsmall;
            tipactualvec(i) = tipactual;
        end
    
    close(figure(l+1))
    figure(l+1)
    plot(epsvec,tipactualvec,'r*')
    hold on
    plot(epsvec,tippredveclarge,'g')
    plot(epsvec,tippredREDveclarge,'k')
    plot(epsvec,tippredREDvecsmall,'b')
    xlabel('\epsilon')
    ylabel('\eta_2')
    end
    print('-f3','slowosc_epscomp_mixed','-djpeg');
    print('-f4','slowosc_epscomp_slow','-djpeg');
end

if value == 2 || value == 3
    %Here we will fix epsilon and compare across lambda
    k = 20;
    lambdavec = linspace(.5,1.5,k);
    tippredveclarge = zeros(1,k);
    tippredREDveclarge = zeros(1,k);
    tippredREDvecsmall = zeros(1,k);
    tipactualvec = zeros(1,k);
    eps = .01;

    parfor i =1:k
        lambda = lambdavec(i);
        Omega = eps^(-lambda);
        A=2;
        B=2;
        h = min(.01,2*pi/(10*Omega));

        %-------------------------------------------------
        % System 
        start = eta1*eta3-1; stop=eta1*eta3+1;
        time = (stop-start)/(eps);

        % Non-dim Stommel Equations
        vDE = @(t,V,T,eta2)((eta1-eta2)-V*abs(V)-T+eta3*(T-V)+A*sin(Omega*t));
        tDE = @(t,V,T,eta2)(eta1-T*(1+abs(V))+B*sin(Omega*t));
        eta2DE = @(t,V,T,eta2)(-eps);

        eta2Init = stop;
        vInit = -stop/eta1+eta3;
        tInit = eta1/(1-vInit);

        % Solve equations for numerical soln, step size h
        [~,Vnum,Tnum,eta2num] = RK2sys3(vDE,tDE,eta2DE,vInit,tInit,eta2Init,0,time,h);

        % Ignore transient behavior
        Vnum = Vnum(100:end);
        Tnum = Tnum(100:end);
        eta2num = eta2num(100:end);

        Min = min(Vnum);
        Max = max(Vnum);

        % Determine the tipping points for each case
        % Determine the tipping points for each case
        Mat=[eta3 1-eta3;eta1 1];
        tippredlarge = eta1*eta3-eps*log(eps)/min(eig(Mat));
        tippredREDlarge = eta1*eta3+eps*log(eps)/(eta1-eta3-eta1*eta3);
        
        scaling = eps^(lambda-1);
        c0 = pi*abs(A)/(eta1*(1-eta3));
        c1 = eta1*(1-eta3)*abs(A)/pi;
        c2 = eta3*pi/(2*eta1*(1-eta3));
        smooth = scaling^(1/3)*c0^(1/3)*(-eps*2.33811);
        oscillatory = c1*(2-c2^2)/Omega;
        
        tippredREDsmall = eta1*eta3+smooth+oscillatory;
        tipactual=eta2num(find(Vnum>criteria,1));

        tippredveclarge(i) = tippredlarge;
        tippredREDveclarge(i) = tippredREDlarge;
        tippredREDvecsmall(i) = tippredREDsmall;
        tipactualvec(i) = tipactual;
    end
    
    close(figure(5))
    figure(5)
    plot(lambdavec,tipactualvec,'r*')
    hold on
    plot(lambdavec,tippredveclarge,'g')
    plot(lambdavec,tippredREDveclarge,'k')
    plot(lambdavec,tippredREDvecsmall,'b')
    xlabel('\lambda')
    ylabel('\eta_2')
    
    print('-f5','slowosc_lambdacomp','-djpeg');
end
