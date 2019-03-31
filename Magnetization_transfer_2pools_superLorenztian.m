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
times = logspace(0, 10000,1e3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%fitting of the data

%super Lorentzian lineshape


%Fit function 2
opts = statset('nlinfit');
opts.MaxIter = 10e5;
opts.TolFun = 1e-30;
opts.TolX = 1e-30;


for i = 1:length(x)
    S(i).v1 = x(i,1);
    S(i).v2 = x(i,2);
end
  
%initial start values based on the Henkelman papers
b1 = [2 176 20 12.9E-6];  %based on henkelman

%fit function
testfun =  @(z,x)   100.*(1*z(1) + ...         
    (x(:,2).^2*pi.*(arrayfun( @(t) (integral(@(theta) (sin(theta) *sqrt(2/pi) *z(4)/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t.v1*z(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), S))')...      
    + 1 + z(2))./(z(1).*(1 + ...
    x(:,2).^2*pi.*(arrayfun( @(t) (integral(@(theta) (sin(theta) *sqrt(2/pi) *z(4)/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t.v1*z(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), S))') +...
    (1 + (x(:,2)./(2.*pi.*x(:,1))).^2.*z(3)).*(1 +  ...   
    (x(:,2).^2*pi.*(arrayfun( @(t) (integral(@(theta) (sin(theta) *sqrt(2/pi) *z(4)/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t.v1*z(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), S))')...   
    + z(2)))...
    *((z(1)>1)&(z(1)<4))...
    *((z(2)>5)&(z(2)<40))...
    *((z(3)>20)&(z(3)<60))...
    *((z(4)>5e-6)&(z(4)<1.1e-5));

[par3,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x, ydata, testfun, b1, opts); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%estimate z-spectrum for each amplitude

w=amp_eql(1);

testfun1 =  100*(1*par3(1) + ...         
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')...      
    + 1 + par3(2))./(par3(1)*(1 + ...
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')) +...   
(1 + (w./(2.*pi.*times)).^2.*(par3(3)))*(1 +  ...   
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')...   
    + par3(2)));


w=amp_eql(2);
testfun2 =  100*(1*par3(1) + ...         
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')...      
    + 1 + par3(2))./(par3(1)*(1 + ...
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')) +...   
(1 + (w./(2.*pi.*times)).^2.*(par3(3)))*(1 +  ...   
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')...   
    + par3(2)));

w=amp_eql(3);
testfun3 =  100*(1*par3(1) + ...         
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')...      
    + 1 + par3(2))./(par3(1)*(1 + ...
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')) +...   
(1 + (w./(2.*pi.*times)).^2.*(par3(3)))*(1 +  ...   
    (w^2.*pi.*(arrayfun(@(t) (integral(@(theta)((sin(theta) *sqrt(2/pi) *par3(4))/(abs(3*(cos(theta)).^2-1))*exp(-2* (2*pi*t*par3(4)./abs(3*(cos(theta)).^2-1)).^2)), 0, 0.5*pi)), times))')...   
    + par3(2)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
