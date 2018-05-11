%%Prelude. 
%% Does Loading and Stuff. 

%% Load and Read Files... N  = 8, Csph. 
%Grids = load('/home/anant/Software/rem3d/rem3d/files/SC_Grid_N7_Jmax4_EulerConfig1.mat');
Grids = load('/home/anant/Software/rem3d/rem3d/files/Csph_Grid_N8_Jmax5_EulerConfig1.mat');

%WaveletMaps_DV_Basis = load('/home/anant/Software/rem3d/rem3d/files/SC_Mmaps.N7.J4.D4.mat');
WaveletMaps_DV_Basis = load('/home/anant/Software/rem3d/rem3d/files/Precon_Mmaps.N8.J5.D4.mat');

%WaveletCoeffstoMakeDV = load('/home/anant/Software/rem3d/rem3d/files/SC_WCoeffs.N7.J4.D4.mat');
WaveletCoeffstoMakeDV = load('/home/anant/Software/rem3d/rem3d/files/Precon_WCoeffs.N8.J5.D4.mat');

%ModelFile = load('SC_MIT_P08.N7.Jmax4.EulerConfig1.mat');
MITFile = load('/home/anant/Software/rem3d/rem3d/files/Csph_MIT_P08.N8.Jmax5.EulerConfig1.mat');
ME16File = load('/home/anant/Software/rem3d/rem3d/files/Csph_ME16_Vp.N8.Jmax5.EulerConfig1.mat');
USSLFile = load('/home/anant/Software/rem3d/rem3d/files/Csph_US-SL-2014.N8.Jmax5.EulerConfig1.mat');
YS_P_H15File = load('/home/anant/Software/rem3d/rem3d/files/Csph_YS-P-H15.N8.Jmax5.EulerConfig1.mat');
NWUSFile = load('/home/anant/Software/rem3d/rem3d/files/Csph_NWUS11-P.N8.Jmax5.EulerConfig1.mat');
DNA13File = load('/home/anant/Software/rem3d/rem3d/files/Csph_DNA13.N8.Jmax5.EulerConfig1.mat');
GYPSUMFile = load('/home/anant/Software/rem3d/rem3d/files/Csph_GYPSUM.N8.Jmax5.EulerConfig1.mat');
GAP_P4File = load('/home/anant/Software/rem3d/rem3d/files/Csph_GAP_P4.N8.Jmax5.EulerConfig1.mat');
THBUSFile = load('/home/anant/Software/rem3d/rem3d/files/Csph_THBUS.N8.Jmax5.EulerConfig1.mat');
