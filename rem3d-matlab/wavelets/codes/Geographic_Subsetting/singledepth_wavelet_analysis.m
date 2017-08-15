


%%%Localized Wavelet Analysis:
lfin = 7;
eo = 0;
write = 0;
[x,y,z,megalon,megalat] = Generate_lat_lon_CubedSphere(lfin,eo,write);
%[vw,depchk,vstats]=loris5AnantVer(wav,N,J,precon,depkm(index),xver,gnorjr,panel);

scale_of_interest = 4;
wavename = 'db4';
center_freq = centfrq(wavename);
pseudo_freq = scal2frq(scale_of_interest,wavename,1);


%%%Go through results 
dep= '0090';
N = '7';
MaxScale = '4';
Basis = 'D4';
homedir = '/home/anant/Software/rem3d/rem3d-matlab/wavelets/';
name = [homedir 'examples/' 'Exam-ple_Output/' 'loris5_GN_' dep '_' N '_' MaxScale '_' Basis '_1_1.mat'];
%name = ['/home/anant/Dropbox/Analyses/wavelets/loris5outputs/' 'loris5_GN_' dep '_' N '_' MaxScale '_' Basis '_1_1.mat'];
result = load(name);
dimensions = size(result.vw);
numvals = 6*(2^(2*lfin));

waveletcoeffs = [];
for i = 1:numvals
    waveletcoeffs = [waveletcoeffs result.vw(i)];
end
abswaveletcoeffs = abs(waveletcoeffs);
scatter(megalon,megalat,10,waveletcoeffs,'filled')
colorbar
caxis([-0.05 0.05])
title(['Absolute Value of Wavelet Coefficients for Depth =' dep 'km'])
xlabel('longitude')
ylabel('latitude')