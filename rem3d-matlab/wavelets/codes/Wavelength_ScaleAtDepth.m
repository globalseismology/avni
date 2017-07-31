function [ Wavelength ] = Wavelength_ScaleAtDepth( N,depth,Scale,wname )
% Calculate the wavelength of an equivalent wavelet.
% Do this given the resolution parameter,and a given depth(km) within the Earth.  
% Note that for Daubechies 4, wname = 'db4'
% Calculates wavelength in km.

% Anant Hariharan

NumPts = 2^(2*N) ;
re = 6371;
PerFaceSurface_Area = (1/6)*((4)*pi*(re*1000-depth*1000)^2);

BoxArea = PerFaceSurface_Area/NumPts;

Length_Dimension = sqrt(BoxArea); %This corresponds to the difference between two samples.

Linear_Wavenumber =scal2frq(Scale,wname,Length_Dimension);

Wavelength = 1/Linear_Wavenumber;
Wavelength = Wavelength/1000;
end

