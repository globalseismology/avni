%%%Localized Wavelet Analysis:
lfin = 7;
[x,y,z]=cube2sphere(lfin);
[vw,depchk,vstats]=loris5AnantVer(wav,N,J,precon,depkm(index),xver,gnorjr,panel);

scale_of_interest = 4;
wavename = 'db4'
center_freq = centfrq(wavename)
pseudo_freq = scal2frq(scale_of_interest,wavename,1)