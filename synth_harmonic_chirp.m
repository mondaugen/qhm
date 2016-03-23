function [x] = synth_harmonic_chirp(f_0,f_1,P,T,Fs),
% SYNTH_HARMONIC_CHIRP: Synthesizes a chirp with harmonically related partials
% Inputs:
% f_0   starting frequency.
% f_1   final frequency.
% P     number of partials
% T     length of output (in seconds)
% Fs    (optional) sample rate (in Hz)
if (nargin == 4),
    Fs=44100;
end
t=(0:1/Fs:T);
x=zeros(1,length(t));
for p=(1:P),
    phi_t=2*pi*p*(f_0*t+(f_1-f_0)/(2*T)*(t.^2));
    x=x+cos(phi_t)/P;
end;
