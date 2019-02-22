function p = hygecdfTall(x,m,k,n,uflag)
%HYGECDF Hypergeometric cumulative distribution function.
%   P = HYGECDF(X,M,K,N) returns the hypergeometric cumulative
%   distribution function with parameters M, K, and N
%   at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   P = HYGECDF(X,M,K,N,'upper') returns the upper tail probability of 
%   the hypergeometric distribution function with parameters M, K, and N
%   at the values in X.
%
%   See also HYGEINV, HYGEPDF, HYGERND, HYGESTAT, CDF.

%   Copyright 1993-2014 The MathWorks, Inc.


if nargin > 4
    uflag = convertStringsToChars(uflag);
end

if nargin < 4,
    error(message('stats:hygecdf:TooFewInputs'));
end
if nargin > 4 
    if ~strcmpi(uflag,'upper')
        error(message('stats:cdf:UpperTailProblem'));
    end
    x = n - floor(x) - 1;
    k = m - k;
end
%[errorcode x m k n] = distchck(4,x,m,k,n);

% if errorcode > 0
%     error(message('stats:hygecdf:InputSizeMismatch'));
% end

%Initialize P to zero.
% if isa(x,'single') || isa(m,'single') || isa(k,'single') || isa(n,'single')
%    p = zeros(size(x),'single');
% else
%    p = zeros(size(x));
% end
p=min(x,0);

k4 = (isnan(x) | isnan(m) | isnan(k) | isnan(n));
p(k4) = NaN;

% Handle values of X for which P is zero by inspection.
k1 = (m - k - n + x + 1 <= 0 | x < 0);

% Handle values of X for which P is one by inspection.
k2 = (x >= n | x >= k);
p(k2) = 1;

% Return NaN for values of the parameters outside their respective limits.
k3 = (m < 0 | k < 0 | n < 0 | round(m) ~= m | round(k) ~= k | round(n) ~= n | ...
      n > m | k > m);
p(k3) = NaN;

kc = ~(k1|k2|k3|k4);

% Compute p when xx >= 0.
if gather(any(kc(:)))
    % For accuracy, for x values that are larger than the mean, compute the
    % upper tail 1-p instead of the lower tail p.  This will tend to
    % increase the number of significant bits in the result.
    lo = (x <= k.*n./m);
    t = kc & lo;
    if gather(any(t(:)))
        p(t) = F(floor(x(t)),m(t),k(t),n(t));
    end
    t = kc & ~lo;
    if gather(any(t(:)))
        p(t) = 1 - F(n(t)-floor(x(t))-1,m(t),m(t)-k(t),n(t));
    end
end

end

% -----------------------------------------
% Hypergeometric cdf without error checking
function p=F(x,m,k,n)
dens = hygepdfTall(x,m,k,n);

% compute hygecdf(x,m,k,n)/hygepdf(x,m,k,n) with a series whose terms can
% be computed recursively, backwards.
xmax = max(x(:));
ybig = repmat((0:xmax)', 1, length(x));
xbig = repmat(x(:)', xmax+1, 1);
mbig = repmat(m(:)', xmax+1, 1);
kbig = repmat(k(:)', xmax+1, 1);
nbig = repmat(n(:)', xmax+1, 1);

terms = ((ybig+1) .* (mbig-kbig-nbig+ybig+1)) ./ ((nbig-ybig) .* (kbig-ybig));
terms(ybig >= xbig) = 1;
terms = cumprod(terms, 'reverse');

terms(ybig > xbig) = 0;
ratio = sum(terms,1);
ratio = reshape(ratio,size(x));

p = ratio.*dens;
  
% Make sure that round-off errors never make P greater than 1.
p(p > 1) = 1;
end
