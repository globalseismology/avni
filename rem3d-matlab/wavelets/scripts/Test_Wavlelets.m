Scale = [1 2 3 4 5 6 7]  ;
N = 7;
wname = 'db4'
Length_Dimension = 100000'

Linear_Wavenumber =scal2frq(Scale,wname,Length_Dimension);
Wavelength1 = 1./Linear_Wavenumber;
Wavelength = Wavelength1./1000

%%%%
depth= 100;
NumPts = 2^(2*N) ;
re = 6371;
PerFaceSurface_Area = (1/6)*((4)*pi*(re*1000-depth*1000)^2);

BoxArea = PerFaceSurface_Area/NumPts;

Length_Dimension = sqrt(BoxArea); %This corresponds to the difference between two samples.


for i = 1:length(Linear_Wavenumber)
   1/Linear_Wavenumber(i)
end