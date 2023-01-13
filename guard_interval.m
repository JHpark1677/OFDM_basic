function y=guard_interval(Ng,Nfft,NgType,ofdmSym)
if NgType==1
    y=[ofdmSym(Nfft-Ng+1:Nfft) ofdmSym(1:Nfft)]; % conducting cyclic prefix
elseif NgType==2
    y=[zeros(1,Ng) ofdmSym(1:Nfft)]; % conducting zero padding
end