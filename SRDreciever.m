clc, clear, close all
   % mysteryA:  SRRCLength:     4
   %            SRRCrolloff:    0.33
   %            T_t:            8.9e-6 s
   %            f_if:           1.6 MHz
   %            f_s:            700 kHz 
   % 
   % mysteryB:  SRRCLength:     5
   %            SRRCrolloff:    0.4
   %            T_t:            7.5e-6 s
   %            f_if:           1.2 MHz
   %            f_s:            950 kHz
   % 
   % mysteryC:  SRRCLength:     3
   %            SRRCrolloff:    0.14
   %            T_t:            8.14e-6 s
   %            f_if:           2.2 MHz
   %            f_s:            819 kHz

str = input('Input the name of mystery (ex. A for mysteryA): ', 's');

% Construct file names based on user input
strMystery = strcat('mystery', str, '.mat');
strTfParams = strcat('mystery', str, '_tf_params.mat');

% Check if the files exist
if ~isfile(strMystery) || ~isfile(strTfParams)
    % Display a message box if either file does not exist
    f = msgbox('One or more files do not exist.');
    return;
end
 load(strMystery);
 load(strTfParams);
rsc = r;
Ts = 1/f_s;
figure('Name',"Received  signal"), plotspec(r, Ts);
%%Using a BPF to get rid of broadband noise
f0 = abs(f_if - round(f_if / f_s)*f_s);
%Since 2*f0 is less than Nyquist frequency , we have to upsample the signal
%so that we can fulfill the Nyquist criteria.

%% Preprocessing with upsampling
if 4*f0 > f_s
    N = 2;
    a = upsample(r,N);
    fs = 2*f_s;
    T_s = 1/fs;
    f1= (2*f0/fs)-0.16;
    f2= (2*f0/fs)-0.15;
    f3= (2*f0/fs)+0.15;
    f4= (2*f0/fs)+0.16;

    fl=500;         % length of the BP filter
    ff=[0 f1 f2 f3 f4 1];  % BPF center frequency at .4
    fa=[0 0 1 1 0 0];                  
    h=firpm(fl,ff,fa);                 % BPF design via firpm
    y=filter(h,1,a);                  % filter for preprocessed r
else 
    y = r;
end
upsampled_signal = y.^2;
figure('Name',"Squared Signal"), plotspec(upsampled_signal, T_s);
%return
% Center frquency used in BPF
fa = (4*f0) / fs 
b = firpm(200,[0 fa-0.02 fa-0.01 fa+0.01 fa+0.02 1], [0 0 1 1 0 0]);
bpfiltered_signal = filter(b,1,upsampled_signal);
figure('Name',"Frequency response of BPF"), freqz(b);
figure('Name',"Estimated Carrier"), plotspec(bpfiltered_signal,T_s);
%return
%% Dual Pll
rp=bpfiltered_signal; % construct carrier = rBPF
t=T_s:T_s:length(rp)*T_s;      % time vector
fl = 100; ff = [0 0.01 .02 1]; fa = [ 1 1 0 0];
h = firpm(fl,ff,fa);
mu1=.04; mu2=0.012;                       % algorithm stepsizes
lent=length(t);
th1=zeros(1,lent);  % initialize estimates
th2=zeros(1,lent); 
z1 = zeros(1,fl+1); z2 = zeros(1,fl+1);
carest = zeros(1,lent);
for k=1:lent-1
  z1= [z1(2:fl+1) , rp(k)*sin(4*pi*f0*t(k)+2*th1(k))];
  update1=fliplr(h)*z1'; % new output of LPF
 
  th1(k+1)=th1(k)- mu1 *update1 ; % algorithm update          % top PLL
  z2= [z2(2:fl+1) , rp(k)*sin(4*pi*f0*t(k)+2*th1(k)+2*th2(k))]; % bottom PLL
  update2=fliplr(h)*z2'; % new output of LPF

  th2(k+1)=th2(k)- mu2 *update2 ; % algorithm update  

  
  carest(k)=cos(2*pi*f0*t(k)+th1(k)+th2(k));                   %  +th2(k)); carrier estimate
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
carest = downsample(carest, 2); 


%% Demodulating the signal to find baseband signal after finding carrier
% if interp == true
%     carest = carest(1:2:length(carest)) ;
% end

figure('Name',"Recovered Carrier Signal"), plotspec(carest, Ts);
demod_sig = r.* carest';
figure('Name',"Demodulated Signal"), plotspec(demod_sig, Ts); grid on;
%Match filtering aswell as LPF
% b = firpm(100,[0 0.15 0.2 1], [1 1 0 0]);
% after_lowpass = filter(b,1,demod_sig);
% figure('Name',"After lowpass"), plotspec(after_lowpass,Ts);
M = T_t * f_s ;   
matchfilt = srrc(srrcLen, beta, M, 0);
filter_signal= filter(matchfilt,1,demod_sig);
figure('Name',"Signal after Match filtering"), plotspec(filter_signal,Ts);


%% interpolation/downsampler with timing recovery
          
len = round(length(r)/M);                                        % Length of the received Signal                                                
l = srrcLen;
M = T_t / Ts;
r_mt = filter_signal;
tnow=l*M+1; tau=0; r_tim=zeros(1,length(r));           % initialize variables
tausave=zeros(1,length(r)); tausave(1)=tau; i=0;
mu=0.005;                                    % algorithm stepsize
delta=0.3;                                  % time for derivative
while tnow<length(r)-l*M                    % run iteration
  i=i+1;
  r_tim(i)=interpsinc(r_mt,tnow+tau,l);           % interpolated value at tnow+tau
  x_deltap=interpsinc(r_mt,tnow+tau+delta,l);  % get value to the right
  x_deltam=interpsinc(r_mt,tnow+tau-delta,l);  % get value to the left
  dx=x_deltap-x_deltam;                     % calculate numerical derivative

  tau=tau+mu*dx*r_tim(i);                      % alg update (energy)
  tnow=tnow+ M; tausave(i)=tau;              % save for plotting
end

% plot results
figure('Name',"Constellation diagram")
subplot(2,1,1), plot(r_tim(1:i-2),'b.')    % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))               % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

%return
%% Training localizer & DDequalizer.m find a LMS equalizer f.
preamble = '0x0 This is is the Frame Header 1y1';       % preamble letters used as training data for equalizer
preambleCode = letters2pam(preamble);                   % conversion of preamble to 4 pam
preambleLen = length(preambleCode);                     % length of preamble
userDataLength = 125;                                   % data length
datalength = userDataLength * 4;                        % total message length
totalFrameLength = datalength + preambleLen;            % total length of one frame
% LMS
n=10;  f=zeros(n,1);         % initialize equalizer at 0
mu=.02; delta=8;             % stepsize and delay delta
for i=n+1:preambleLen        % iterate
  rr=r_tim(i:-1:i-n+1)' ;    % vector of received signal
  e=preambleCode(i-delta)-rr'*f;        % calculate error
  f=f+mu*e*rr;               % update equalizer coefficients
end
equalizedSignal = filter(f,1,r_tim);                    % applying filter coefficient to the downsampled signal
figure('Name',"Constellation diagram after equalization"), plot(equalizedSignal,'b.')
grid on
title('Constellation diagram');


%% frame synchronization
quantizedSignal = quantalph(equalizedSignal, [-3, -1, 1, 3])';
y = conv(flip(preambleCode), quantizedSignal);
[m,ind] = max(y);
headstart = length(quantizedSignal(i))-ind+1+preambleLen;
figure('Name',"Correlation"), subplot(3,1,1), stem(preambleCode)             % plot header
title('Header')
subplot(3,1,2), stem(r_tim(1:500))             % plot data sequence
title('Data with embedded header')
subplot(3,1,3), stem(y(1:500))                 % plot correlation
title('Correlation of header with data')
[value,firstIndex] = max(abs(y(1:200)));       % taking first peak in correlated signal 
ind = [firstIndex, ];       
for i = 2:(length(y))/totalFrameLength         % updating arrays for total number of frame and indices of peak
    firstIndex = firstIndex + totalFrameLength;
    ind(i) = firstIndex;    
end
updatedIndex = y(ind);                         % to work with phase shift of 180
msgLength = ind(2)-ind(1)-preambleLen;         % message length
receivedMsg = [];                              % initialization of array to receive message
for i = 1:length(ind)-1
    updateMsg = quantizedSignal(ind(i)+1:ind(i)+msgLength); % frame data
    if(updatedIndex(i) > 0)     
        receivedMsg = [receivedMsg, updateMsg];                % integrate the frame data                                                                            %     in the correlated signal (A & B) half of signal the peamble siignal has minus values        
    else 
        receivedMsg = [receivedMsg, -updateMsg];               % integrate -180 phase shifted data       
    end
end
reconstructed_msg = pam2letters(receivedMsg)      % decode message