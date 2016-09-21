function C = flt(FUN,n,DOMAIN)

% FLT.m Fast Legendre transform (part of QUINOPT)
%
% C = FLT(FUN,n,[a,b]) computes the first n Legendre coefficients of the
% expansion of the function FUN defined over the domain [a,b]. 
% Inputs:
%
%     FUN: a function handle to the function to be projected onto the
%          Legendre polynomials.
%       n: a positive integer.
%   [a,b]: a 2-by-1 real row vector, i.e. DOMAIN = [a,b] with a<b.
%
% This function implements Algorithm 2 from "A fast and simple algorithm 
% for the computation of Legendre coefficients", Iserles, Numer. Math.(2010)

% ----------------------------------------------------------------------- %
%        Author:    Giovanni Fantuzzi
%                   Department of Aeronautics
%                   Imperial College London
%       Created:    05/05/2016
% Last Modified:    05/05/2016
% ----------------------------------------------------------------------- %

if ~isa(FUN,'function_handle')
    error('Input FUN should be a valid function handle')
elseif ~isnumeric(n) || numel(n)~=1 || rem(n,1)
    error('Input n should be a positive integer')
elseif ~isnumeric(DOMAIN)||~isvector(DOMAIN)||numel(DOMAIN)~=2|| ...
        DOMAIN(1)>=DOMAIN(2)||any(DOMAIN==Inf)|| any(DOMAIN==-Inf)
    error('Domain must be a valid compact domain [a,b] with a<b.')
end

% Parameters from Iserles 2010
M = 100; 
N = 2^(ceil( log(n+3+2*M)/log(2) ));

% tolerance to clean coefficients from roundoff noise
tol = 1e-16;

% Find matrix g
g = zeros(N-3-2*M+1,M+1);
g(1,1)=1;
for m = 1:N-3-2*M
    % First col of g
    g(m+1,1) = m/(m-0.5)*g(m,1);
end
for j = 1:M
    % Rest of g
    m = (0:N-3-2*M)';
    g(m+1,j+1) = (m+j).*(j-0.5)./j./(m+j+0.5).*g(m+1,j);
end

% Form discrete cosine transform
% line with fft maybe faster, but less accurate from experiments
k = (0:N-1)';
zeta = cos(2*pi*k/N);                                         % points for the cos transform in [-1,1]
x =   0.5*( (DOMAIN(2)-DOMAIN(1))*zeta+DOMAIN(2)+DOMAIN(1) ); % rescale points to DOMAIN=[a,b]
fval = FUN(x);       
sigma = zeros(N,1);
for m = 0:N-1
    sigma(m+1) = 2/N*sum( fval.*cos(2*pi*k*m/N ));
end
%sigma = real(exp(-1i.* pi.*k./(2*N)).*fft(fval)).*(2/N);

% Compute coefficients
C = zeros(n,1);
for m = 0:n-1
    C(m+1)= g(m+1,1:M+1)*( sigma(m+2.*(0:M)+1)-sigma(m+2.*(0:M)+3) );
end
C = 0.5.*C(:);
C = C.*(C<=-tol|C>=tol);