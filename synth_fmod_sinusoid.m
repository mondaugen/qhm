function [x]=synth_fmod_sinusoid(f0,phi0,fm,Am,phim_0,Fs,T)
t=(0:(1/Fs):T);
x=exp(1i*(phi0+Am./fm.*cos(phim_0)+2*pi*f0.*t-Am./fm.*cos(2*pi*fm*t+phim_0)));
