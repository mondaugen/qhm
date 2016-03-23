function [f,A,phi,a,b,p1,p2] = aqhm(s,f_i,A_i,phi_meas,T,H,w)
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
if (nargin < 7),
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
K=size(f_i,1);
N=length(s);
n=(1:N);
L=floor(((N+2*T)-M)/H)+1;
m=((1:L)-1)*H+1;
%size(f)
%size(A)
phi=zeros(K,N);
phi_extr=zeros(K,N);
a=zeros(K,L);
b=zeros(K,L);
p1=zeros(K,L);
p2=zeros(K,L);
% aQHM
C=10;
c=1;
while (1),
    f=interp1(m,f_i,n,'cubic');
    A=interp1(m,A_i,n,'cubic');
    for l=(1:(L-1)),
        phi_i=phi_meas(:,l);
        for u=(1:H),
%            size(phi_extr)
%            size(phi_i)
            phi_extr(:,(l-1)*H+u)=phi_i;
            phi_i = phi_i + f(:,(l-1)*H+u)*2*pi;
        end
    end
    acc=phi_meas(:,1);
    for l=(1:(L-1)),
%        acc=phi_meas(:,l);
        m_=round(abs(phi_meas(:,l+1)-phi_extr(:,l*H+1))/(2*pi));
        a_=pi*(phi_meas(:,l+1)+2*pi*m_-phi_extr(:,l*H+1))/(2*H);
        for u=(1:H),
            phi(:,(l-1)*H+u)=acc;
            acc=acc+2*pi*f(:,(l-1)*H+u)+a_*sin(pi*(u-1)/H);
        end
%        mod(phi_meas(:,l+1),2*pi)
%        mod(acc,2*pi)
    end
    if c==(C+1),
        break;
    end
    % Here iterate only over hops where the window will have phase information
    % from which to derive the functions
    for l=((floor(T/H)+2):(floor((N-(T+1))/H)+1)),
        t=(-T:T)+((l-1)*H+1);
        E=exp(1i*(phi(:,t)-phi(:,(l-1)*H+1)).');
%        plot(t,E.');
        Et=(-T:T).'.*E;
        E=w.*E;
        Et=w.*Et;
        E_Et=[E Et];
        x=(E_Et'*E_Et)^(-1)*(E_Et'*(s(t).'));
        a(:,l)=x(1:K);
        b(:,l)=x((K+1):(2*K));
        p1(:,l)=(real(a(:,l)).*real(b(:,l))+imag(a(:,l)).*imag(b(:,l)))./(abs(a(:,l)).^2);
        p2(:,l)=(real(a(:,l)).*imag(b(:,l))-imag(a(:,l)).*real(b(:,l)))./(abs(a(:,l)).^2);
        A_i(:,l)=abs(a(:,l));
        phi_meas(:,l)=arg(a(:,l));
        p2(:,l)./(2*pi)
        f_i(:,l)=f_i(:,l)+p2(:,l)./(2*pi);
%        f_i(:,l)=f0;
    end
    c=c+1;
end
