function [ truncated_wavelets ] = zero_wavelets(vw,indices,waveletbasis,N,Jmin,Jmax )
% Takes a binary array of indices and an array of wavelets, same dimensions
% Zeros the wavelet in the wavelet array if the corresponding index is 0. 
% Basically, keeps the wavelets indicated within the indices array!

defval('N',7);
defval('Jmin',4);
defval('Jmax',4);
defval('waveletbasis','D4');
defval('precon',[1 1]);
numvals = 6*2^(2*N);
vw2 = vw;
for i = 1:numvals
    if indices(i) == 0
        vw2(i)=0;
    end
end

truncated_wavelets = vw2;

     
end