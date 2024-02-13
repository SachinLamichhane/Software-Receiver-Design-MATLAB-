%% to find the best mu
rp = bpfiltered_signal;  % construct carrier = rBPF
t = T_s:T_s:length(rp)*T_s;  % time vector
fl = 100; ff = [0 0.01 .02 1]; fa = [ 1 1 0 0];
h = firpm(fl, ff, fa);
lent = length(t);
th1 = zeros(1, lent);  % initialize estimates
th2 = zeros(1, lent);
carest = zeros(1, lent);

% Define the range and step size for mu1 and mu2
mu1_values = 0:0.0005:1;  % Modify the range and step size accordingly
mu2_values = 0:0.01:1;         % Modify the range and step size accordingly

best_mu1 = 0;
best_mu2 = 0;
best_error = inf;

for mu1 = mu1_values
    for mu2 = mu2_values
        z1 = zeros(1, fl+1);
        z2 = zeros(1, fl+1);

        for k = 1:lent-1
            z1 = [z1(2:fl+1) , rp(k)*sin(4*pi*f0*t(k)+2*th1(k))];
            update1 = fliplr(h) * z1';
            th1(k+1) = th1(k) - mu1 * update1;

            z2 = [z2(2:fl+1) , rp(k)*sin(4*pi*f0*t(k)+2*th2(k)+2*th1(k))];
            update2 = fliplr(h) * z2';
            th2(k+1) = th2(k) - mu2 * update2;

            carest(k) = cos(2*pi*f0*t(k)+th1(k)+th2(k));
        end

        % Calculate the error between preprocessed carrier and estimated carrier
        error = bpfiltered_signal(1:length(carest)) - carest';

        % Update best values if the current combination is better
        if norm(error) < best_error
            best_error = norm(error);
            best_mu1 = mu1;
            best_mu2 = mu2;
        end
    end
end

fprintf('Best mu1: %f\n', best_mu1);
fprintf('Best mu2: %f\n', best_mu2);
return