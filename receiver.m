clc 
clear 
close all
  % mysteryA:   SRRCLength:     4
  %             SRRCrolloff:    0.33
  %             T_t:            8.9e-6 s
  %             f_if:           1.6 MHz
  %             f_s:            700 kHz 
  % 
  %  mysteryB:  SRRCLength:     5
  %             SRRCrolloff:    0.4
  %             T_t:            7.5e-6 s
  %             f_if:           1.2A MHz
  %             f_s:            950 kHz
  % 
  %  mysteryC:  SRRCLength:     3
  %             SRRCrolloff:    0.14
  %             T_t:            8.14e-6 s
  %             f_if:           2.2 MHz
  %             f_s:            819 kHz
str = input('Input the name of mystery.(ex.A for mysteryA):','s');
strMystery = append('mystery',str,'.mat');
strTfParams = append('mystery',str,'_tf_params.mat');
if(isfile(strMystery)==0 || isfile(strTfParams)==0)
    f = msgbox('File does not exist.');
    return;
end
load(strMystery);
load(strTfParams);
%% Received signal plot
a = r;
T_s = 1/f_s;                                        %Sampling time of the recieved signal                 
figure("Name","Received Signal A"), plotspec(r(1:100),T_s)
%% aliasing condition
fc = zeros(1,10);
for k = 1:length(fc)
    f = abs(f_if - k * f_s);
    fc(k) = f;
end
fc = min(fc);
%% Upsampling
m = 2;                      % oversampling factor
N = length(r);              % total length of the received signal
rup = zeros(1,N*m);         % zero padding
rup(1:m:N*m) = r;           % oversampling the received signal by m times
fs = m * f_s;               % increasing sampling frequency by m times
Ts = 1/fs;                  % new sampling period 
% f0 = 2*fc/fs;             % center frequency of band passed filter
% fl = 500; ff = [0, f0-0.03, f0-0.02, f0+0.02, f0+0.03, 1]; fa = [0, 0, 1, 1, 0, 0];             % bandpass filter parameters
% h = firpm(fl,ff,fa);
% y = filter(h,1,rup);
rsc = rup.^2;          % squaring the oversampled signal(non linearity)

figure("Name","Square of Received Signal"), plotspec(rsc,Ts);
%return
%% PLL preprocessing
f0 = 4*fc/fs;
fl = 700; ff = [0, f0-0.02, f0-0.01, f0+0.01, f0+0.02, 1]; fa = [0, 0, 1, 1, 0, 0];             % bandpass filter parameters
b = firpm(fl,ff,fa);
rp = filter(b,1,rsc);
figure("Name", "Estimated Carrier Frequency"),plotspec(rp,Ts);
figure("Name","Frequency Response of Bandpass Filter"), freqz(b)
%% dual pll: Adapt the duallplls.m
f0=fc;                                   % center frequency of bandpass filter
t=0:Ts:length(rp)*Ts-Ts; 
mu1=.01; mu2=.0002;                       % algorithm stepsizes
lent=length(t); th1=zeros(1,lent);       % initialize estimates
th2=zeros(1,lent); carest=zeros(1,lent);
for k=1:lent-1
  th1(k+1)=th1(k)-mu1*rp(k)*sin(4*pi*f0*t(k)+2*th1(k));           % top PLL
  th2(k+1)=th2(k)-mu2*rp(k)*sin(4*pi*f0*t(k)+2*th1(k)+2*th2(k));  % bottom PLL
  carest(k)=cos(2*pi*f0*t(k)+th1(k)+th2(k));                  % carrier estimate
end
figure("Name","Dual PLL"), 
subplot(2,1,1), plot(t,th1)              % plot first theta
title('output of first PLL')
ylabel('\theta_1')
subplot(2,1,2), plot(t,th2)              % plot second theta
title('output of second PLL')
ylabel('\theta_2')

carest=downsample(carest,2);  %downsampling
figure("Name","Recovered Carrier"), plotspec(carest, T_s)
demod_sig = a .* carest';                  % Demodulate the signal using carrier
figure("Name","Demodulated Signal"), plotspec(demod_sig, T_s)
%% Match filtering &&  LPF
M = f_s*T_t;   
matchfilt = srrc(SrrcLen, beta, M, 0);
filter_signal= filter(matchfilt,1,demod_sig);
figure('Name',"Signal after Match filtering"), plotspec(filter_signal(1:100),T_s);
% %return
% %% interpolation/downsampler with timing recovery          
len = round(length(r)/M);                                        % Length of the received Signal                      
r_tim = zeros(1,len);                           % initialize down smapled signal 
l = SrrcLen;
% timing-recovery algorithm using the recursive output-power-maximization algorithm
tnow = l*M+1; tau=0;                                    % initialized Variables
tausave = zeros(1,len); tausave(1)=tau; i=0;                            
mu = 0.06;                                             % algorithm stepsize
delta = 0.3;                                           % time for derivate
r_mt = filter_signal;
while tnow < length(r)-l*M                                   % run iteration
  i=i+1;
  r_tim(i) = interpsinc(r_mt, tnow+tau, l);             % interpolated value at tnow+tau
  x_deltap = interpsinc(r_mt, tnow+tau+delta, l);       % value to right
  x_deltam = interpsinc(r_mt, tnow+tau-delta, l);       % value to left
  dx = x_deltap - x_deltam;                             % numerical derivative                      % alg update
  tau = tau + mu*dx*r_tim(i);                           % alg update (energy)
  tnow = tnow + M; 
  tausave(i)=tau;                                       % save for plotting
end
% plot results
figure('Name',"Constellation diagram"), subplot(2,1,1), plot(r_tim(1:i-2),'b.')        
title('Constellation diagram');
ylabel('Estimated Symbols');
grid on
subplot(2,1,2), plot(tausave(1:i-2))                      
ylabel('off_set estimates'), xlabel('iterations')

 %% DD Timing recovery
% l= SrrcLen;
% M = f_s * T_t;
% n = round(length(r)/M);
% tnow=l*M+1; tau=0; r_tim=zeros(1,n);           % initialize variables
% r_mt = filter_signal;
% tausave=zeros(1,n); tausave(1)=tau; i=0;
% mu=0.008;                                    % algorithm stepsize
% delta=0.3;                                  % time for derivative
% while tnow<length(r)-2*l*M                  % run iteration
%   i=i+1;
%   r_tim(i)=interpsinc(r_mt,tnow+tau,l);           % interpolated value at tnow+tau
%   x_deltap=interpsinc(r_mt,tnow+tau+delta,l);  % get value to the right
%   x_deltam=interpsinc(r_mt,tnow+tau-delta,l);  % get value to the left
%   dx=x_deltap-x_deltam;                     % calculate numerical derivative
%   qx=quantalph(r_tim(i),[-3,-1,1,3]);          % quantize xs to nearest 4-PAM symbol
%   tau=tau+mu*dx*(qx-r_tim(i));                 % alg update: DD
%   tnow=tnow+M; tausave(i)=tau;              % save for plotting
% end
% figure('Name',"Constellation diagram"), subplot(2,1,1), plot(r_tim(1:i-2),'b.')        
% title('Constellation diagram');
% ylabel('Estimated Symbols');
% grid on
% subplot(2,1,2), plot(tausave(1:i-2))                      
% ylabel('off_set estimates'), xlabel('iterations')
%  %% Training localizer & DDequalizer.m find a LMS equalizer f.
% preamble = '0x0 This is is the Frame Header 1y1';                     % preamble letters used as training data for equalizer
% preambleCode = letters2pam(preamble);
% preambleLen = length(preambleCode);
% userDataLength = 125;                                                 % 125 PAM symbol
% frameLength = userDataLength * 4;
% totalFrameLength = frameLength + preambleLen;
% 
% % LMS
% n = 11;
% mu=.02; delta=5;             % stepsize and delay delta
% f = zeros(n,1);
% for i=n+1:preambleLen                 % iterate
%   rr=r_tim(i:-1:i-n+1)' ;        % vector of received signal
%   e=preambleCode(i-delta)-rr'*f;        % calculate error
%   f = f+mu*e*rr; % update equalizer coefficients
% end
% 
% equalizedSignal = filter(f,1,r_tim);
% equalizedSignal = equalizedSignal(delta+1:end);
% % figure("Name","Equalized Signal"), plotspec(equalizedSignal,T_s)
% figure('Name',"Constellation diagram after LMS"), plot(equalizedSignal,'b.')
% grid on
% title('Constellation diagram after LMS');

%% DD
n=20; f=zeros(n,1);           % initialize equalizer
f(round(n/2)) = 1;
mu=.00001;                       % stepsize
for i=n+1:m                  % iterate
  rr=r_tim(i:-1:i-n+1)';         % vector of received signal
  e=quantalph(f'*rr,[-3, -1, 1, 3 ])-f'*rr;       % calculate error
  f=f+mu*e*rr;               % update equalizer coefficients
  % f(i) = f;
end
equalizedSignal = filter(f,1,r_tim);
equalizedSignal = equalizedSignal(delta+1:end);
% figure("Name","Equalized Signal"), plotspec(equalizedSignal,T_s)
figure('Name',"Constellation diagram after DD"), plot(equalizedSignal,'b.')
grid on
title('Constellation diagram after DD');
% Equalizer Testing
% % finaleq=f;                   % test final filter f        % output of channel
% m = 1000;
% yt=filter(f,1,r_tim);            % use final filter f to test
% dec=sign(real(yt));          % quantization
% for sh=0:n                   % if equalizer is working, one
%   err(sh+1)=0.5*sum(abs(dec(sh+1:m)-r_tim(1:m-sh)));
% end                          % of these delays has zero error
% err
% 
% [hb,w]=freqz(b,1);
% [hf,w]=freqz(f,1);
% [hc,w]=freqz(conv(b,f),1);
% figure("Name","Equalizer Testing")
% semilogy(w,abs(hb))
% hold on
% semilogy(w,abs(hf),'r')
% semilogy(w,abs(hc),'g')
% semilogy(w,abs(hb).*abs(hf),'k')
% grid on
% hold off
%% frame synchronization
preamble = '0x0 This is is the Frame Header 1y1';                     % preamble letters used as training data for equalizer
preambleCode = letters2pam(preamble);
preambleLen = length(preambleCode);
userDataLength = 125;                                                 % 125 PAM symbol
frameLength = userDataLength * 4;
totalFrameLength = frameLength + preambleLen;
quantizedSignal = quantalph(equalizedSignal, [-3, -1, 1, 3])';
y = conv(flip(preambleCode), quantizedSignal);
[m,ind] = max(y);
%headstart = length(quantizedSignal(i))-ind+1+preambleLen;
% figure('Name',"Correlation"), subplot(3,1,1), stem(preambleCode)             % plot header
% title('Header')
% subplot(3,1,2), stem(r_tim(1:500))             % plot data sequence
% title('Data with embedded header')
% subplot(3,1,3), stem(y(1:500))                % plot correlation
% title('Correlation of header with data')
figure('Name',"Correlation"), stem(y)
grid on
p = conv(flip(preambleCode), preambleCode);     % finding the peak of conolution of preamble
[value,firstIndex] = max(abs(y(1:200)));                % taking first peak in correlated signal 
% secondIndex = firstIndex + totalFrameLength;    % taking second index
ind = [firstIndex, ];
for i = 2:(length(y))/totalFrameLength
    firstIndex = firstIndex + totalFrameLength;
    ind(i) = firstIndex;    
end
%count = 0;
% for i = 1:length(y)
%     if (abs(y(i)) > (value-10))
%         count = count+1;
%         ind(count) = i;
%     end
% end
updatedIndex = y(ind);                       % to work with phase shift of 180
msgLength = ind(2)-ind(1)-preambleLen;
receivedMsg = [];

for i = 1:length(ind)-1
    updateMsg = equalizedSignal(ind(i)+1:ind(i)+msgLength);
    if(updatedIndex(i) > 0)
        receivedMsg = [receivedMsg, updateMsg];                % integrate the frame data                                                                            %     in the correlated signal (A & B) half of signal the peamble siignal has minus values  
      
    else 
        receivedMsg = [receivedMsg, -updateMsg];
       
    end
end
receivedMsg = quantalph(receivedMsg, [-3 -1 1 3])';
% %% DD
% n=8; f=zeros(n,1);           % initialize equalizer
% f(round(n/2)) = 1;
% mu=.006;                       % stepsize
% for i=n+1:length(receivedMsg)                  % iterate
%   rr=receivedMsg(i:-1:i-n+1)';         % vector of received signal
%   e=quantalph(f'*rr,[-3, -1, 1, 3 ])-f'*rr;       % calculate error
%   f=f+mu*e*rr;               % update equalizer coefficients
%   % f(i) = f;
% end
% equalizedSignal = filter(f,1,receivedMsg);
% equalizedSignal = equalizedSignal(delta+1:end);
% % figure("Name","Equalized Signal"), plotspec(equalizedSignal,T_s)
% figure('Name',"Constellation diagram after DD"), plot(equalizedSignal,'b.')
% % cvar = (receivedMsg -z) * (receivedMsg -z)'/length(receivedMsg); 
% % lmp = length(receivedMsg);
% % d = r(1:lmp);
% % pererr = 100 * sum(abs(sign(receivedMsg - d)))/lmp;
% receivedMsg = quantalph(receivedMsg, [-3 -1 1 3])';
receivedmsg = pam2letters(receivedMsg)