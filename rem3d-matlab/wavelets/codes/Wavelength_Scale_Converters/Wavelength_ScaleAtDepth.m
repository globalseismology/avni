N = 7;
wname = 'db4';
depth = 0;
Scale = 4;
% Calculate the wavelength of an equivalent wavelet.
% Do this given the resolution parameter,and a given depth(km) within the Earth.  
% Note that for Daubechies 4, wname = 'db4'
% Calculates wavelength in km.
% Example: [ Wavelength,Linear_Wavenumber,GcDeg ] = Wavelength_ScaleAtDepth( 7,0,4,'db4')


% Anant Hariharan
FREQ = centfrq(wname);
NumBox = 2^(2*N) ;  %assuming odd gridding.  %number
re = 6371*1000; %m
PerFaceSurface_Area = (1/6)*((4)*pi*((re-depth*1000)^2)); %m^2

BoxArea = PerFaceSurface_Area/NumBox; %m^2

Length_Dimension = (BoxArea)^(0.5); %m %This corresponds to the difference between two samples.
Wavelength = [];
for i = 1:length(Scale)
    Wavelength(i) = (Scale*Length_Dimension)/FREQ
end

Wavelength = Wavelength/1000;
GcDeg = km2deg(Wavelength)


