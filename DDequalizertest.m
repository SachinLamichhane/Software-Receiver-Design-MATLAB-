
m = 10e7;
s = pam(m,4,5);
b = [1, 1, -0.9, -0.2, 0.2, 1];
r = filter(b,1,s);
n = 37; 
f = [zeros(1,floor(n/2)) 1 zeros(1,floor(n/2))]';
% f = zeros(n,1);
% f(round(n/2)) = 1; % center spike initialization of equalizer
mu=1*10e-6; % stepsize
%m = length(r);
for i=n+1:m                  % iterate
  rr=r(i:-1:i-n+1)';         % vector of received signal
  e=quantalph(f'*rr,[-3,-1,1,3])-f'*rr;       % calculate error
  f=f+mu*e*rr;               % update equalizer coefficients
end
equalizedSignal = filter(f,1,r);
%equalizedSignal = equalizedSignal(1:end);
figure('Name',"Constellation Diagram after Equalization"), 
plot(equalizedSignal(1:length(equalizedSignal)),'b.')