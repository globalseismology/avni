%% Get RMS of Different Scales for a Single Model- how much does incrementally 
%% adding shorter-wavelength structure change the model. 

% Anant Hariharan
clear all
clc
close all


%% Load and Read Files... N  = 8, Csph. 
%Grids = load('/home/anant/Software/rem3d/rem3d/files/SC_Grid_N7_Jmax4_EulerConfig1.mat');
Grids = load('/Users/Anant/Desktop/Education/Research/IRIS_Research/UMD_GlobalTomo/Thresholding_Implementation/Csph_Grid_N8_Jmax5_EulerConfig1.mat');

%WaveletMaps_DV_Basis = load('/home/anant/Software/rem3d/rem3d/files/SC_Mmaps.N7.J4.D4.mat');
WaveletMaps_DV_Basis = load('/Users/Anant/Desktop/Education/Research/IRIS_Research/UMD_GlobalTomo/Thresholding_Implementation/Precon_Mmaps.N8.J5.D4.mat');

%WaveletCoeffstoMakeDV = load('/home/anant/Software/rem3d/rem3d/files/SC_WCoeffs.N7.J4.D4.mat');
WaveletCoeffstoMakeDV = load('/Users/Anant/Desktop/Education/Research/IRIS_Research/UMD_GlobalTomo/Thresholding_Implementation/Precon_WCoeffs.N8.J5.D4.mat');

%ModelFile = load('SC_MIT_P08.N7.Jmax4.EulerConfig1.mat');
MITFile = load('/Users/Anant/Desktop/Education/Research/IRIS_Research/UMD_GlobalTomo/Thresholding_Implementation/Csph_MIT_P08.N8.Jmax5.EulerConfig1.mat');
ME16File = load('/Users/Anant/Desktop/Education/Research/IRIS_Research/UMD_GlobalTomo/Thresholding_Implementation/Csph_ME16_Vp.N8.Jmax5.EulerConfig1.mat');
USSLFile = load('/Users/Anant/Desktop/Education/Research/IRIS_Research/UMD_GlobalTomo/Thresholding_Implementation/Csph_US-SL-2014.N8.Jmax5.EulerConfig1.mat');

%% Set Params
% Note, Africa is on face 3. 
% Note, US is on face 1??
N=8;
LonMin = -125;
LonMax = -65;
LatMin = 31.7;
LatMax = 50;
Scale_List = [6 5 4 3 2 1]; 

redo = 1
if redo

for iiii = 1:length(Scale_List)
Desired_Scales = Scale_List(1:iiii)



%% Do Thresholding using the "Grids" file as a geographic reference. This example only does rectangular.
Geog_Indices = find(Grids.lon > LonMin & Grids.lon < LonMax & Grids.lat > LatMin & Grids.lat<LatMax);

re = 6371;


MITModelDepths = unique(MITFile.depth);
USSLModelDepths = unique(USSLFile.depth);
ModelDepthsFlip = flipud(USSLModelDepths);

MIT_RMS = [];
ME16_RMS = [];
USSL_RMS = [];



%% Now, find the model coefficients that correspond to the points you desire. Put in a structure. 
Data.Indices = Geog_Indices;
% First, extract the coefficients for the right depth. 
% MIT

for iiiii = 1:length(ModelDepthsFlip)
desired_radius = ModelDepthsFlip(iiiii)
MITwvAtDepth = MITFile.wvcoeffs(find(MITFile.depth == desired_radius));
MITveloAtDepth = MITFile.model(find(MITFile.depth == desired_radius));
%ME16
ME16wvAtDepth = ME16File.wvcoeffs(find(ME16File.depth == desired_radius));
ME16veloAtDepth = ME16File.model(find(ME16File.depth == desired_radius));
%Schmandt
USSLwvAtDepth = USSLFile.wvcoeffs(find(USSLFile.depth == desired_radius));
USSLveloAtDepth = USSLFile.model(find(USSLFile.depth == desired_radius));

%Data.velocoeffs

for i = 1:length(Data.Indices)
    Data.MITWvcoeffs(i) = MITwvAtDepth(Data.Indices(i));
    Data.ME16Wvcoeffs(i) = ME16wvAtDepth(Data.Indices(i));
    Data.USSLWvcoeffs(i) = USSLwvAtDepth(Data.Indices(i));
    
    Data.MITVelocoeffs(i) = MITveloAtDepth(Data.Indices(i));
    Data.ME16Velocoeffs(i) = ME16veloAtDepth(Data.Indices(i));
    Data.USSLVelocoeffs(i) = USSLveloAtDepth(Data.Indices(i));
end

%% Now, make a map of the wavelet coefficients that correspond to the thresholded region. 
%% Note: This is essentially a 'selective' forward transform into the wavelet domain. 

Reconstructed_MITWavelet_Map = zeros(2^N,2^N);
Reconstructed_ME16Wavelet_Map = zeros(2^N,2^N);
Reconstructed_USSLWavelet_Map = zeros(2^N,2^N);

for i = 1:length(Data.Indices)
    CurrMITMap = Data.MITVelocoeffs(i)*full(WaveletCoeffstoMakeDV.Me(Data.Indices(i)).dv);
    CurrME16Map = Data.ME16Velocoeffs(i)*full(WaveletCoeffstoMakeDV.Me(Data.Indices(i)).dv);
    CurrUSSLMap = Data.USSLVelocoeffs(i)*full(WaveletCoeffstoMakeDV.Me(Data.Indices(i)).dv);    
    Reconstructed_MITWavelet_Map = Reconstructed_MITWavelet_Map+CurrMITMap;
    Reconstructed_ME16Wavelet_Map = Reconstructed_ME16Wavelet_Map+CurrME16Map;
    Reconstructed_USSLWavelet_Map = Reconstructed_USSLWavelet_Map+CurrUSSLMap;
end

%% Now separate out by scales

MetaJmax =  ME16File.MetaJmax;
vwlev1=cube2scale(N,[MetaJmax MetaJmax]+1,1,1);

Reconstructed_MITWavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_ME16Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_USSLWavelet_Map_ScaleInspected = zeros(2^N,2^N);

for ii = 1:length(Desired_Scales)
    curr_scale = Desired_Scales(ii);
    Scale_check = find(vwlev1 == curr_scale);
    Reconstructed_MITWavelet_Map_ScaleInspected(Scale_check) = Reconstructed_MITWavelet_Map(Scale_check);
    Reconstructed_ME16Wavelet_Map_ScaleInspected(Scale_check) = Reconstructed_ME16Wavelet_Map(Scale_check);
    Reconstructed_USSLWavelet_Map_ScaleInspected(Scale_check) = Reconstructed_USSLWavelet_Map(Scale_check);
end

%% Now do the inverse transform! 
%%
Reconstructed_MITVeloMap = zeros(2^N,2^N);
Reconstructed_ME16VeloMap = zeros(2^N,2^N);
Reconstructed_USSLVeloMap = zeros(2^N,2^N);

for i = 1:2^(2*N)
    CurrMITMap = Reconstructed_MITWavelet_Map_ScaleInspected(i)*full(WaveletMaps_DV_Basis.Me(i).dv);
    CurrME16Map =  Reconstructed_ME16Wavelet_Map_ScaleInspected(i)*full(WaveletMaps_DV_Basis.Me(i).dv);
    CurrUSSLMap =  Reconstructed_USSLWavelet_Map_ScaleInspected(i)*full(WaveletMaps_DV_Basis.Me(i).dv);

    Reconstructed_MITVeloMap = Reconstructed_MITVeloMap+CurrMITMap;
    Reconstructed_ME16VeloMap = Reconstructed_ME16VeloMap+CurrME16Map;
    Reconstructed_USSLVeloMap = Reconstructed_USSLVeloMap+CurrUSSLMap;
end


Curr_MITRMS = sqrt(mean(mean(Reconstructed_MITVeloMap.^2)));
Curr_ME16RMS = sqrt(mean(mean(Reconstructed_ME16VeloMap.^2)));
Curr_USSLRMS = sqrt(mean(mean(Reconstructed_USSLVeloMap.^2)));
MIT_RMS = [MIT_RMS Curr_MITRMS];
ME16_RMS = [ME16_RMS Curr_ME16RMS];
USSL_RMS = [USSL_RMS Curr_USSLRMS];


end

RMS_Struct(iiii).MIT = MIT_RMS;
RMS_Struct(iiii).ME16 = ME16_RMS;
RMS_Struct(iiii).USSL = USSL_RMS;

end
else
    
%LOAD

end

figure()    

Scale6=plot(RMS_Struct(1).MIT,6371-ModelDepthsFlip,'-ro');
hold on
Scale5=plot(RMS_Struct(2).MIT,6371-ModelDepthsFlip,'-bo');
Scale4=plot(RMS_Struct(3).MIT,6371-ModelDepthsFlip,'-go');
Scale3=plot(RMS_Struct(4).MIT,6371-ModelDepthsFlip,'-co');
Scale2=plot(RMS_Struct(5).MIT,6371-ModelDepthsFlip,'-mo');
Scale1=plot(RMS_Struct(6).MIT,6371-ModelDepthsFlip,'-ko');

grid on
box on
hXLabel=xlabel('Tomographic Model RMS')
hYLabel=ylabel('Depth (km)')
set(gca,'ydir','reverse')
hTitle=title('Model RMS With Depth Over the US for MITP08')
p=patch([0 1 1 0],[410 410 660 660],'r','FaceAlpha',0.2);
hLegend = legend([Scale6,Scale5,Scale4,Scale3,Scale2,Scale1],{'Scaling Function','Up to Scale 5','Up to Scale 4','Up to Scale 3','Up to Scale 2','Up to Scale 1'})
xlim([0 0.2])
ylim([50 1200])
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:200:2500, ...
  'LineWidth'   , 1         );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

save('ScaleBasedRMS_Struct.mat','RMS_Struct')


