function [f,A,phi,a,b,p1,p2] = qhm(s,f0,T,H,w)
% QHM    quasi-harmonic model
% Inputs:
% s      signal to be analysed.
% f0     K initial frequency estimates.
% T      1/2 window length - 1.
% H      Hop size.
% w      (optional) window of length 2*T + 1. If not provided, uses Hann window.
% Outputs:
% f      K x L frequency estimates where L = floor((N-(2*T+1))/H)+1 and N is the
%        number of samples.
% A      K x L amplitude estimates.
% phi    K x L phase estimates.
% a      K x L a coefficient estimates
% b      K x L b coefficient estimates
%
% if the first sample of s is time 0, the values in the columns of f, a and phi
% refer to times H, 2*H, ... L*H.
M=(2*T+1);
if (nargin < 5),
    if mod(M,2)==0,
        w=hamming(M+1);
        w=w(1:M);
    else
        w=hamming(M);
    end
else
    w=w(:);
    if (length(w) != M),
        error('Window "w" must have length 2*T + 1.');
    end
end
s=s(:).';
s=[zeros(1,T) s zeros(1,T)];
f0=f0(:);
K=length(f0);
N=length(s);
L=floor((N-M)/H)+1;
f=zeros(K,L);
A=zeros(K,L);
phi=zeros(K,L);
a=zeros(K,L);
b=zeros(K,L);
p1=zeros(K,L);
p2=zeros(K,L);
% QHM
for l=(1:L),
    t=(-T:T)+((l-1)*H+T+1);
    E=exp(2*pi*1i*(-T:T).'*f0.');
    Et=(-T:T).'.*E;
    E=w.*E;
    Et=w.*Et;
    E_Et=[E Et];
    x=(E_Et'*E_Et)^(-1)*(E_Et'*(s(t).'));
    a(:,l)=x(1:K);
    b(:,l)=x((K+1):(2*K));
    p1(:,l)=(real(a(:,l)).*real(b(:,l))+imag(a(:,l)).*imag(b(:,l)))./(abs(a(:,l)).^2);
    p2(:,l)=(real(a(:,l)).*imag(b(:,l))-imag(a(:,l)).*real(b(:,l)))./(abs(a(:,l)).^2);
    A(:,l)=abs(a(:,l));
    phi(:,l)=arg(a(:,l));
    f0=f0+p2(:,l)./(2*pi);
    f(:,l)=f0;
end
