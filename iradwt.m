function y = iradwt(w,p,q,s)
% y = iradwt(w,p,q,s)
% Inverse RAtional-Dilation Wavelet Transform: an overcomplete
% rational-dilation wavelet transform implemented using the FFT.
%
% See also radwt

J = length(w)-1;
y = w{J+1};
for j = J:-1:1
    N = length(w{j})*s;
    [H, G] = MakeFreqResp(N, p, q, s);
    y = sfb(y(1:N*p/q), w{j}, H, G, p, q, s);
end
