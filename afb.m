function [lo, hi] = afb(x, H, G, p, q, s)
% [lo, hi] = afb(x, H, G, p, q, s)
% Analsysis Filter Bank: implements the analysis filter bank for
% the overcomplete rational-dilation wavelet transform using the FFT.
% INPUT
%   x - input signal
%   H, G - low-pass and high-pass frequency response
%   q/p - dilation factor
%   s - high-pass channel rate
% OUTPUT
%   lo, hi - low-pass, high-pass output signals
% NEED
%   length(x) = multiple of lcm(q,s)
%   length(H) = length(x)
%   length(G) = length(x)
%   redundancy = 1/s * 1/(1-p/q) > 1
%   (No error checking done)
%
% Ilker Bayram and Ivan Selesnick
% Polytechnic Institute, New York
% November 2008
%
% % Example (check perfect reconstruction)
% p = 4;
% q = 5;
% s = 2;
% N = 7*q*s;
% [H,G] = MakeFreqResp(N, p, q, s);
% x = rand(1,N);
% [lo,hi] = afb(x, H, G, p, q, s);
% y = sfb(lo, hi, H, G, p, q, s);
% y = y(1:N);
% max(abs(x - y))    % reconstruction error ~ 0

% This function assumes the frequency responses, H and G,
% are the ones produced by 'MakeFreqResp'.

X = fft(x);

% low-pass subband
Xlo = H.*X;
N = length(X)*p/q;
if rem(N,2) == 0
    % N is even
    Xlo = [Xlo(1:N/2)  Xlo(end-N/2+1:end)];
else
    % N is odd
    Xlo = [Xlo(1:(N+1)/2)  Xlo(end-(N-1)/2+1:end)];
end
lo = ifft(Xlo/q);

% high-pass subband
Xhi = G.*X;
N = length(Xhi)/s;
% use the fact that G is zero except from (s-1)pi/s to (s+1)pi/s
if rem(s,2) == 0
    % s is even
    Xhi = Xhi((1:N)+(s/2-1)*N) + Xhi((1:N)+s/2*N);
else
    % s is odd
    Xhi = Xhi((1:N)+(s-1)/2*N);
end
hi = ifft(Xhi/s);

