function y = hygepdfTall(x,m,k,n)
%HYGEPDF Hypergeometric probability density function.
%   Y = HYGEPDF(X,M,K,N) returns the hypergeometric probability
%   density function at X with integer parameters M, K, and N.
%   Note: The density function is zero unless X is an integer.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   See also HYGECDF, HYGEINV, HYGERND, HYGESTAT, PDF.

%   Reference:
%      [1]  Mood, Alexander M., Graybill, Franklin A. and Boes, Duane C.,
%      "Introduction to the Theory of Statistics, Third Edition", McGraw Hill
%      1974 p. 91.

%   Copyright 1993-2010 The MathWorks, Inc.


if nargin < 4,
    error(message('stats:hygepdf:TooFewInputs'));
end

%[errorcode x m k n] = distchck(4,x,m,k,n);
% 
% if errorcode > 0
%     error(message('stats:hygepdf:InputSizeMismatch'));
% end

% Initialize Y to zero.
% if isa(x,'single') || isa(m,'single') || isa(k,'single') || isa(n,'single')
%    y = zeros(size(x),'single');
% else
%    y = zeros(size(x));
% end
y=min(x,0);
y(isnan(x) | isnan(m) | isnan(k) | isnan(n)) = NaN;

% Return NaN for values of the parameters outside their respective limits.
k1 = (m < 0 | k < 0 | n < 0 | round(m) ~= m | round(k) ~= k ...
                    | round(n) ~= n | n > m | k > m);
y(k1) = NaN;

% Remove values of X for which Y is zero by inspection.
k2 = (m - k - n + x + 1 <= 0 | x < 0 | round(x) ~= x | x > n | x > k);

% At n=0 or n==m must have pdf=1
k3 = (n==0 | n==m) & ~k1 & ~k2;
y(k3) = 1;

% Find integer values of x that are within the correct range
kc = ~(k1 | k2 | k3);
if gather(any(kc(:)))
    x = x(kc);
    m = m(kc);
    k = k(kc);
    n = n(kc);
    if ~isfloat(x), x = double(x); end
    if ~isfloat(m), m = double(m); end
    if ~isfloat(k), k = double(k); end
    if ~isfloat(n), n = double(n); end
    
    lnsqr2pi = 0.9189385332046727;
    
    % Get the biggest term for m-choose-n.
    % Use binomial pdf with p=n/m. Since p=n/m, binomial deviance is zero.
    p = n./m;
    mn = -lnsqr2pi -0.5*log(n.*(1-n./m)) ...
        +stirlerr(m) -stirlerr(n) -stirlerr(m-n);
    
    ix0 = x==0;
    ixk = x==k;
    kx = zeros(size(x));
    kx(ix0) = k(ix0).*log(1-p(ix0));
    kx(ixk) = x(ixk).*log(p(ixk));
    ok = ~ix0 & ~ixk;
    kx(ok) = -lnsqr2pi -0.5*log(x(ok).*(1-x(ok)./k(ok))) ...
        +stirlerr(k(ok)) -stirlerr(x(ok)) -stirlerr(k(ok)-x(ok)) ...
        -binodeviance(x(ok),k(ok).*p(ok)) ...
        -binodeviance(k(ok)-x(ok),k(ok).*(1-p(ok)));
    
    inx0 = n==x;
    inxmk = (n-x)==(m-k);
    mknx = zeros(size(x));
    mknx(inx0) = (m(inx0)-k(inx0)).*log(1-p(inx0));
    mknx(inxmk) = (n(inxmk)-x(inxmk)).*log(p(inxmk));
    ok = ~inx0 & ~inxmk;
    mknx(ok) = -lnsqr2pi -0.5*log((n(ok)-x(ok)).*(1-(n(ok)-x(ok))./(m(ok)-k(ok)))) ...
        +stirlerr(m(ok)-k(ok)) -stirlerr(n(ok)-x(ok)) -stirlerr(m(ok)-k(ok)-n(ok)+x(ok)) ...
        -binodeviance(n(ok)-x(ok),(m(ok)-k(ok)).*p(ok)) ...
        -binodeviance(m(ok)-k(ok)-n(ok)+x(ok),(m(ok)-k(ok)).*(1-p(ok)));
    
    y(kc) = exp(kx + mknx - mn);
end
