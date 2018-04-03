function f=angularD4WT_BasisVer(N,x,Jmax,tipe)

f = zeros(2^N,2^N);

VeloBasisName = ['/home/anant/Software/rem3d/rem3d/files/Precon_Mmaps.N' num2str(N) '.J' num2str(Jmax) '.D4.mat'];
wletBasisName = ['/home/anant/Software/rem3d/rem3d/files/Precon_WCoeffs.N' num2str(N) '.J' num2str(Jmax) '.D4.mat'];

if tipe == 'inverse'
    WaveletMaps_DV_Basis = load(VeloBasisName);    
    for i = 1:2^(2*N)    
        CurrMap = x(i)*full(WaveletMaps_DV_Basis.Me(i).dv);    
        f = f+CurrMap;
    end    
end

if tipe == 'forward'
     WaveletCoeffstoMakeDV = load(wletBasisName);    
    for i = 1:2^(2*N)    
        CurrMap = x(i)*full(WaveletCoeffstoMakeDV.Me(i).dv);    
        f = f+CurrMap;
    end
    
end


end
