% close all
close all, clear all

 m = 33931;
 a = pam(m,4,5);
 h = [0.5 0.8 -0.6 1 -0.8];
% %h = [0 1 0];
 r_tim = filter(h,1,a);
n=15; f=zeros(n,1);           % initialize equalizer
f(round(n/2)) = 1;
mu=.005;                       % stepsize
for i=n+1:length(r_tim)                  % iterate
  rr=r_tim(i:-1:i-n+1)';         % vector of received signal
  e=quantalph(f'*rr,[-3, -1, 1, 3 ])-f'*rr;       % calculate error
  f=f+mu*e*rr;               % update equalizer coefficients
  %f(i) = f;
end
equalizedSignal = filter(f,1,r_tim);
%figure("Name","Equalized Signal"), plotspec(equalizedSignal,T_s)
figure('Name',"Constellation diagram after Equalizer"), plot(equalizedSignal,'b.')
%%
finaleq = f;

 m = 1000;
 s = pam(m,4,5);
a = filter(h, 1, s); % Output of channel
c = filter(f,1,a);


% Quantization
dec = quantalph(real(c),[-3, -1, 1,3])';

% Test different delays and calculate the error
err = zeros(1, n+1);
for sh = 0:n
    err(sh+1) = 0.5*sum(abs(dec(sh+1:m) - s(1:m-sh)));
    
end
err
% Display the errors for different delays
% disp('Equalizer Error for Different Delays:');
% disp(err);
figure('Name',"Frequency Responses")
[hb,w] = freqz(h,1);
[hf,w] = freqz(f,1);
[hc,w] = freqz(conv(h,f),1);
semilogy(w,abs(hb))
hold on
semilogy(w,abs(hf),'r')
semilogy(w,abs(hc),'g')
%semilogy(w,abs(hb).*abs(hf),'k')
legend("channel","equalizer","combination")
hold off
figure('Name',"Convolution"), plot(conv(h,f));

