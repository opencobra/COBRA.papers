function [ess] = effectiveSampleSize(x)
%ess = effectiveSampleSize(x)
%compute the effective sample sizes of each parameter in x
%
%Input:
% x - a dim x N vector where N is the length of the chain
%
%Output:
% ess - either dim x 1 vector where ess(i) is the effective sample size of x(i,:)

N = size(x, 2); d = size(x, 1);
Neven = N - mod(N,2);
ess = zeros(d, 1);

% compute ess over each row
for i = 1:d
   % normalize i-th row
   x_ = x(i,:);
   x_ = x_ - mean(x_);
   v = mean(x_.^2);
   if (v < 1e-16 * (mean(x(i,:).^2)+1)) % when x is a constant, ess = N
      ess(i) = N;
      continue;
   end
   
   x_ = x_ / sqrt(v);
   
   % compute autocorrelation via Wiener-Khinchin theorem
   r = ifft(abs(fft(x_,2*N)).^2); % power spectral density
   ac = real(r(1:Neven))/N; % ess formula assume vector length is even
   
   % Geyer's monotone estimator
   minAC = cummin(ac(:, 1:2:end) + ac(:, 2:2:end));
   ess(i) = N/max(1,2*sum(minAC.*(minAC>0)) -1);
end