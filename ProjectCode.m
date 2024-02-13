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
  %             f_if:           1.2 MHz
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
    f = msgbox('File is not exist.');
    return;
end

load(strMystery);
load(strTfParams);

%% Received signal plot

T_s = 1/f_s;                                        %Sampling time of the recieved signal                 
figure("Name","Received Signal"), plotspec(r,T_s)

%% aliasing condition
fc = zeros(1,10);
for k = 1:length(fc)
    f = abs(f_if - k * f_s);
    fc(k) = f;
end
fc = min(fc);

%% Upsampling
if f_if >  2*f_s
m = 3;                      % oversampling factor
N = length(r);              % total length of the receieved signal
rup = zeros(1,N*m);         % zero padding
rup(1:m:N*m) = r;           % oversampling the received signal by m times
fs = m * f_s;               % increasing sampling frequency by m times
Ts = 1/fs;                  % new sampling period 

rsc = rup.^2;               % squaring the oversampled signal(non linearity)
figure("Name","Square of Received Signal"), plotspec(rsc,Ts)

elseif f_if < 2*f_s
    m = 2;                      % oversampling factor
    N = length(r);              % total length of the receieved signal
    rup = zeros(1,N*m);         % zero padding
    rup(1:m:N*m) = r;           % oversampling the received signal by m times
    fs = m * f_s;               % increasing sampling frequency by m times
    Ts = 1/fs;                  % new sampling period 

    rsc = rup.^2;               % squaring the oversampled signal(non linearity)
    figure("Name","Square of Received Signal"), plotspec(rsc,Ts)    
end
%% Bandpass filtering of squared signal
f0 = m*fc*2/fs;             % center frequency of band passed filter
fl = 500; ff = [0, f0-0.018, f0-0.0099, f0+0.0099, f0+0.018, 1]; fa = [0, 0, 1, 1, 0, 0];             % bandpass filter parameters
b = firpm(fl,ff,fa);
figure("Name","Frequency Response of Bandpass Filter"), freqz(b)
rp = filter(b,1,rsc);
figure("Name", "Estimated Carrier Frequency"),plotspec(rp,Ts);