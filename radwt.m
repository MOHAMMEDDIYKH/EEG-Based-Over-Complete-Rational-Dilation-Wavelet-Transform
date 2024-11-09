function w = radwt(x,J,p,q,s)
% w = radwt(x,J,p,q,s)
% RAtional-Dilation Wavelet Transform: an overcomplete rational-dilation
% wavelet transform implemented using the FFT.
% INPUT
%   x - input signal
%   J - number of levels
%   q/p - dilation factor
%   s - high-pass sampling factor
% OUTPUT
%   w - wavelet coefficients
% NOTES
%   w{j} is subband j for j = 1,...,J+1
%   w{1} is the highest-frequency subband signal
%   w{J+1} is the low-pass subband signal
% NEED
%   redundancy = 1/s * 1/(1-p/q) > 1
%   (No error checking done)
%  
% % Example (check perfect reconstruction)
% p = 4; q = 5; s = 2; J = 3;   % parameters
% N = 200;                      % signal length
% x = rand(1,N);                % test signal
% w = radwt(x, J, p, q, s);     % wavelet transform
% y = iradwt(w, p, q, s);       % inverse wavelet transform
% y = y(1:N);                   % remove trailing zeros
% max(abs(x - y))               % reconstruction error ~ 0
%
% Ilker Bayram and Ivan Selesnick
% Polytechnic Institute, New York
% November 2008
% 
% Reference: Frequency-Domain Design of Overcomplete
% Rational-Dilation Wavelet Transforms, submitted to 
% IEEE Transactions on Signal Processing, Nov 2008.
%
% http://taco.poly.edu/selesi/rational/

w = cell(1,J+1);
for j = 1:J,    
    % Zero-pad signal x to next multple of q*s
    N = length(x);    
    C = lcm(q,s);
    L = C * ceil(N/C);
    xp = [x zeros(1,L-N)];
    
    % Apply analysis filter bank
    [H, G] = MakeFreqResp(L, p, q, s);
    [x, w{j}] = afb(xp, H, G, p, q, s);   
end
w{J+1} = x;


