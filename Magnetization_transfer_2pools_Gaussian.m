%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2-pool Magnetization Transfer model, assuming Gaussing or super-Lorenzian
%lineshapes in the style of Henkelman 1991 and later papers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%data to fit
%experimental data at 3 different pulse amplitudes, expressed in flip angles,
%as constant wave-equivalent, as in the Ramani 2002 papers
amp_eql = [194 129 65]; 
%frequency offsets compared to Larmor frequency for each pulse amplitude set (in Hz)
fr1 = [10000 1710 210 1280 4280 430]';
fr2= [860 1280 4280 10000 210 1710 430]';
fr3= [210 10000 860 1280 4280 430 1710]';
%Magnetization transfer ratio from experimental data
mtr1=[0.97 0.89 0.23 0.85 0.95 0.49]';
mtr2=[0.88 0.92 0.98 0.99 0.33 0.95 0.67]';
mtr3=[0.63 1 0.97 0.98 1 0.90 0.99]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%parameter tidy up
amp1=amp_eql(1).*ones(length(mtr1),1);
amp2=amp_eql(2).*ones(length(mtr2),1);
amp3=amp_eql(3).*ones(length(mtr3),1);
fr=[fr1; fr2 ;fr3];
amp=[amp1; amp2; amp3];
mtr=[mtr1; mtr2; mtr3];
x = [fr, amp];  
ydata=mtr;
times = logspace(0, 10000,1e5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%fitting of the data

%Gaussian lineshape
Rb=1;
testfun =  @(z,x) z(5).*(Rb*z(1) + (x(:,2).^2.*sqrt(pi/2).*z(4).*exp( -(2.*pi.*x(:,1).*z(4)).^2/2))...
    + Rb + z(2))./ (z(1).*(Rb + (x(:,2).^2.*sqrt(pi/2).*z(4).*exp( -(2.*pi.*x(:,1).*z(4)).^2/2)))...
    + (1 + (x(:,2)./(2.*pi.*x(:,1))).^2.*z(3)).*(1 + (x(:,2).^2.*sqrt(pi/2).*z(4).*...
    exp( -(2.*pi.*x(:,1).*z(4)).^2/2)) + z(2)));

%initial start values based on the Henkelman papers
b1 = [2 176 20 12.9E-6 400];  

%fit function
opts = statset('nlinfit');
opts.MaxIter = 10e5;
opts.TolFun = 1e-30;
opts.TolX = 1e-30;
[par3,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x, ydata, testfun, b1, opts); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%plot results
%estimate z-spectrum for each amplitude
w=amp_eql(1);
Rrfb1 = w^2*sqrt(pi/2)*par3(4)*exp(-(2*pi*par3(4).*times).^2/2);
testfun1 =  par3(5).*(1*par3(1) + Rrfb1 + 1 + par3(2)) ./ (par3(1)*(1 + Rrfb1) + (1 + (w./(2*pi*times)).^2.*(par3(3))).*(1 + Rrfb1 + par3(2)));

w=amp_eql(2);
Rrfb2 = w^2*sqrt(pi/2)*par3(4)*exp(-(2*pi*par3(4).*times).^2/2);
testfun2 =  par3(5).*(1*par3(1) + Rrfb2 + 1 + par3(2)) ./ (par3(1)*(1 + Rrfb2) + (1 + (w./(2*pi*times)).^2.*(par3(3))).*(1 + Rrfb2 + par3(2)));

w=amp_eql(3);
Rrfb3 = w^2*sqrt(pi/2)*par3(4)*exp(-(2*pi*par3(4).*times).^2/2);
testfun3 =  par3(5).*(1*par3(1) + Rrfb3 + 1 + par3(2)) ./ (par3(1)*(1 + Rrfb3) + (1 + (w./(2*pi*times)).^2.*(par3(3))).*(1 + Rrfb3 + par3(2)));

%plot it in log space
figure()
semilogx(fr1,mtr1, 'ro', times, testfun1, 'r-',...
    fr2, mtr2, 'bo', times, testfun2, 'b-',...
    fr3,mtr3, 'go', times, testfun3, 'g-');
title('Measurements and fitted model')
legend('Measurements', 'Fitted model', 'Location', 'SouthEast')
xlabel('Frequency offset [Hz]')
ylabel('MT/M0')
axis([0 1e5 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
