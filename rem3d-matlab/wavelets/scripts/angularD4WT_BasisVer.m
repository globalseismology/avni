function f=angularD4WT_BasisVer(N,x,Jmax,tipe)

f = zeros(2^N,2^N);

if tipe == 'inverse'
    WaveletMaps_DV_Basis = load('/home/anant/Software/rem3d/rem3d/files/Precon_Mmaps.N8.J5.D4.mat');    
    for i = 1:2^(2*N)    
        CurrMap = x(i)*full(WaveletMaps_DV_Basis.Me(i).dv);    
        f = f+CurrMap;
    end    
end

if tipe == 'forward'
     = load('/home/anant/Software/rem3d/rem3d/files/Precon_WCoeffs.N8.J5.D4.mat');    
    for i = 1:2^(2*N)    
        CurrMap = x(i)*full(WaveletCoeffstoMakeDV.Me(i).dv);    
        f = f+CurrMap;
    end
    
end


end
