rp=bpfiltered_signal; % construct carrier = rBPF
t=0:T_s:length(rp)*T_s-T_s;      % time vector
fl = 100; ff = [0 0.01 .02 1]; fa = [ 1 1 0 0];
h = firpm(fl,ff,fa);
mu1=.001; mu2=.005;                       % algorithm stepsizes
lent=length(r);
% th1=zeros(1,lent);  % initialize estimates
th1(1) = -0.1;
th2(1) = 0.3;
% th2=zeros(1,lent); 
z1 = zeros(1,fl+1); z2 = zeros(1,fl+1);
for k=1:length(rp)-1
  z1= [z1(2:fl+1) , rp(k)*sin(4*pi*f0*(k*Ts)+2*th1(k))];
  update1=fliplr(h)*z1'; % new output of LPF
  th1(k+1)=th1(k)- mu1 *update1 ; % algorithm update          % top PLL
  z2= [z2( 2:fl+1) , rp(k)*sin(4*pi*f0*(k*Ts)+2*th2(k))]  ;
  update2=fliplr(h)*z2'; % new output of LPF
  th2(k+1)=th2(k)- mu2 *update2 ; % algorithm update   % bottom PLL
  carest(k)=cos(2*pi*f0*t(k)+2*th1(k)+2*th2(k));                  % carrier estimate
end
figure('Name',"DPLL")
subplot(3,1,1), plot(t,th1)              % plot first theta
title('output of first PLL')
ylabel('\theta_1')
subplot(3,1,2), plot(t,th2)              % plot second theta
title('output of second PLL')
ylabel('\theta_2')
subplot(3,1,3), plot(bpfiltered_signal(1:length(carest))-carest') % plot difference between estimate
                                % and preprocesssed carrier rp
title('error between preprocessed carrier and estimated carrier')
xlabel('time')
ylabel('f_0 - f')