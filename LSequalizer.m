%b = [0.5 1 -0.6];
% b = [1 1 -0.9 -0.2 0.2 1];
% pcode = letters2pam('0x0 This is is the Frame Header 1y1');
% m = length(pcode);
% %h = filter(b,1,s);
% r = r_tim;                   % output of channel
% %n=100;                               % length of equalizer - 1
% %delta=9;                           % use delay <=n*length(b)
% p=length(pcode)-delta;
% R=toeplitz(r(n+1:p),r(n+1:-1:1));  % build matrix R
% S=pcode(n+1-delta:p-delta)';           % and vector S
% f=inv(R'*R)*R'*S;                 % calculate equalizer f
% Jmin=S'*S-S'*R*inv(R'*R)*R'*S;     % Jmin for this f and delta
% y=filter(f,1,r);                   % equalizer is a filter
% dec=quantalph(y, [-3 ,-1, 1, 3])';                       % quantize and find errors
% err=0.5*sum(abs(dec(delta+1:m)-pcode(1:m-delta)))
% %figure('Name',"Constellation diagram after Equalizer"), plot(y,'b.')

% LMS test
preamble = '0x0 This is is the Frame Header 1y1';                     % preamble letters used as training data for equalizer
preambleCode = letters2pam(preamble);
preambleLen = length(preambleCode);
userDataLength = 125;                                                 % 125 PAM symbol
frameLength = userDataLength * 4;
totalFrameLength = frameLength + preambleLen;
% LMS
%n=15; f=zeros(n,1);           % initialize equalizer at 0
mu=.01; %delta=10;             % stepsize and delay delta
m = length(r_tim);
for i=n+1:preambleLen                  % iterate
  rr=r_tim(i:-1:i-n+1)' ;        % vector of received signal
  e=preambleCode(i-delta)-rr*f';        % calculate error
  f=f+mu*e'*rr;               % update equalizer coefficients
  p=length(preambleCode)-delta;
  S=preambleCode(n+1-delta:p-delta)';
  R=toeplitz(r(n+1:p),r(n+1:-1:1));
  Jmin=S'*S-S'*R*inv(R'*R)*R'*S;
end
equalizedSignal = filter(f,1,r_tim);
dec = quantalph(equalizedSignal,[-3, -1,1,3]);
err=0.5*sum(abs(dec(delta+1:m)-r_tim(1:m-delta)'))