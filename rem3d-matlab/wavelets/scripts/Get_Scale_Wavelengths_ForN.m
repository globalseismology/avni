%%Extract max wavelengths and degrees from d4 Wavelets. 
% Make Pretty? Plot!

N = 8;
depkm = 0;
Jmax = 5;

Res_Matrix = zeros(Jmax,Jmax);


for maxJ = 1:Jmax
    for curr_J = 1:maxJ
        [maxdeg,WavelengthKm] = Scale2Wavelength_Depth( N,curr_J,Jmax,depkm);
        Res_Matrix(maxJ,curr_J) = WavelengthKm;
    end
    maxJ
end


image(Res_Matrix)
title(['Wavelength for Scales with Differing Max Scale: N =' num2str(N) 'depth =' num2str(depkm)])
xlabel('Maximum J')
ylabel('Scale of Interest')
colorbar