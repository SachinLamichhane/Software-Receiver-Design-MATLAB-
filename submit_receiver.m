clc 
clear all
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
figure("Name","Received Signal"), plotspec(r,T_s)
%% aliasing condition
fc = zeros(1,10);                                   % initialize a zero padded array to store aliasing frequecies
for k = 1:length(fc)
    f = abs(f_if - k * f_s);                        % kth aliasing ffrequency
    fc(k) = f;                                      % store the aliasing frequency
end
fc = min(fc);                                       % minimum aliasing frequency
%% Upsampling
m = 2;                      % oversampling factor
N = length(r);              % total length of the receieved signal
rup = ones(1,N*m);         % zero padding
rup(1:m:N*m) = r;           % oversampling the received signal by m times
fs = m * f_s;               % increasing sampling frequency by m times
Ts = 1/fs;                  % new sampling period 
rsc = rup.^2;                 % squaring the oversampled signal(non linearity)
figure("Name","Square of Received Signal"), plotspec(rsc,Ts);       % plot spectrum of squared signal

%% PLL preprocessing
f0 = 4*fc/fs;               % center frequency of squared signal for bandpass filter(bpf)
fl = 700; ff = [0, f0-0.02, f0-0.01, f0+0.01, f0+0.02, 1]; fa = [0, 0, 1, 1, 0, 0];     % bpf parameters        % bandpass filter parameters
b = firpm(fl,ff,fa);        % bpf coefficient
rp = filter(b,1,rsc);       % apply bpf to squared signal
figure("Name", "Estimated Carrier Frequency"),plotspec(rp,Ts);      % plot spectrum
figure("Name","Frequency Response of Bandpass Filter"), freqz(b)    % plot frequency response


%% dual pll: Adapt the duallplls.m
f0=fc;                                   % estimated freqeuncy
t=0:Ts:length(rsc)*Ts-Ts;                 % time vector
mu1=.04; mu2=.00012;                     % algorithm stepsizes
lent=length(t); th1=zeros(1,lent);       % initialize phase estimate of firstt pll
th2=zeros(1,lent); carest=zeros(1,lent); % initialize phase estimate of 2nd pll
for k=1:lent-1
  th1(k+1)=th1(k)-mu1*rp(k)*sin(4*pi*f0*t(k)+2*th1(k));           % top PLL
  th2(k+1)=th2(k)-mu2*rp(k)*sin(4*pi*f0*t(k)+2*th1(k)+2*th2(k));  % bottom PLL
  carest(k)=cos(2*pi*f0*t(k)+th1(k)+th2(k));                  % carrier estimate
end
figure("Name","Dual PLL"), 
subplot(3,1,1), plot(t,th1)              % plot first theta
title('output of first PLL')
ylabel('\theta_1')
subplot(3,1,2), plot(t,th2)              % plot second theta
title('output of second PLL')
ylabel('\theta_2')
subplot(3,1,3), plot(rsc(1:length(carest))-carest) % plot difference between estimate
                                % and preprocesssed carrier rp
title('error between preprocessed carrier and estimated carrier')
xlabel('time')
ylabel('f_0 - f')
% downsampling of upsampled signal an demodulation
carest=carest(1:m:length(carest));  %downsampling
figure("Name","Recovered Carrier"), plotspec(carest, T_s)
carest = carest';
demod_sig = r.* carest;                  % Demodulate the signal using carrier
figure("Name","Demodulated Signal"), plotspec(demod_sig, T_s)
%% Match filtering &&  LPF
M = f_s*T_t;                                % oversampling factor
matchfilt = srrc(SrrcLen, beta, M, 0);       
filter_signal= filter(matchfilt,1,demod_sig);
figure('Name',"Signal after Match filtering"), plotspec(filter_signal,T_s);
%% interpolation/downsampler with timing recovery          
len = floor(length(r)/M);                                        % Length of the received Signal                      
r_tim = zeros(1,len);                           % initialize down smapled signal 
l = SrrcLen;
% timing-recovery algorithm using the recursive output-power-maximization algorithm
tnow = l*M+1; tau=0;                                    % initialized Variables
tausave = zeros(1,len); tausave(1)=tau; i=0;                            
mu = 0.009;                                             % algorithm stepsize
delta = 0.3;                                           % time for derivate
r_mt = filter_signal;
while tnow < length(r)-l*M                                   % run iteration
  i=i+1;
  r_tim(i) = interpsinc(r_mt, tnow+tau, l);             % interpolated value at tnow+tau
  x_deltap = interpsinc(r_mt, tnow+tau+delta, l);       % value to right
  x_deltam = interpsinc(r_mt, tnow+tau-delta, l);       % value to left
  dx = x_deltap - x_deltam;                             % numerical derivative                     
  tau = tau + mu*dx*r_tim(i);                           % alg update (energy)
  tnow = tnow + M; 
  tausave(i)=tau;                                       % save for plotting
end
% plot results

figure('Name',"Constellation diagram"), subplot(2,1,1), plot(r_tim,'b.')        
title('Constellation diagram');
ylabel('Estimated Symbols');
grid on
subplot(2,1,2), plot(tausave(1:i-2))                      
ylabel('off_set estimates'), xlabel('iterations')

%% Training localizer & DDequalizer.m find a LMS equalizer f.
preamble = '0x0 This is is the Frame Header 1y1';       % preamble letters used as training data for equalizer
preambleCode = letters2pam(preamble);                   % conversion of preamble to 4 pam
preambleLen = length(preambleCode);                     % length of preamble
userDataLength = 125;                                   %125 PAM symbol
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
equalizedSignal = filter(f,1,r_tim);                  % applying filter coefficient to the downsampled signal
% equalizedSignal = equalizedSignal(delta+1:end);
%figure("Name","Frequency Response of Equalizer"), freqz(f)
%figure("Name","Equalized Signal"), plotspec(equalizedSignal,T_s)
figure('Name',"Constellation diagram"), plot(equalizedSignal,'b.')
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
% Apply equalizer to each frame
numFrames = floor(length(r_tim) / totalFrameLength);
receivedMsg = zeros(1, datalength);  % Initialize array to store the received message

for frameIndex = 1:numFrames
    startIndex = (frameIndex - 1) * totalFrameLength + 1;
    endIndex = frameIndex * totalFrameLength;
    
    % Extract the frame to be equalized
    frameToEqualize = r_tim(startIndex:endIndex);
    
    % Apply the equalizer to the frame
    equalizedFrame = filter(f, 1, frameToEqualize);
    
    % Store the equalized frame in the received message array
    receivedMsg((frameIndex - 1) * userDataLength + 1 : (frameIndex - 1) * userDataLength + length(equalizedFrame(preambleLen + 1 : end))) = equalizedFrame(preambleLen + 1 : end);
end

% ... (Further processing of the receivedMsg if needed)

% Example: Convert the received message to letters
quantizedSignal = quantalph(receivedMsg, [-3, -1, 1, 3]);
receivedMessage = pam2letters(quantizedSignal)'

% for i = 1:length(ind)-1
%     updateMsg = quantizedSignal(ind(i)+1:ind(i)+msgLength); % frame data
%     if(updatedIndex(i) > 0)     
%         receivedMsg = [receivedMsg, updateMsg];                % integrate the frame data                                                                            %     in the correlated signal (A & B) half of signal the peamble siignal has minus values        
%     else 
%         receivedMsg = [receivedMsg, -updateMsg];               % integrate -180 phase shifted data       
%     end
% end
     % decode message
