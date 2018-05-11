% clear; close all;
MIT_CumRMS = zeros(1,6);
ME16_CumRMS = zeros(1,6);
USSL_CumRMS = zeros(1,6);
DNA13_CumRMS = zeros(1,6);
NWUS_CumRMS = zeros(1,6);
%YS_P_H15_RMS = zeros(1,6);
Ved_Scott_CumRMS = zeros(1,6);
GYPSUM_CumRMS = zeros(1,6);
GAP_P4_CumRMS = zeros(1,6);

% Select the depth at which to do comparisons
% target_dep = 500; 
target_rad = 6371-target_dep;

% To what wavelet scale should I analyze?
Jmax = 5; %This isn't really an option!

% Load a colomap
load('ved_tomo.mat');
N=8;
% Load in a whole bunch of stuff
Prelude

Scale_List = [6 5 4 3 2 1]; 
ModelDepths = unique(USSLFile.depth);
ModelDepthsFlip = flipud(ModelDepths);

% load polygons...
Load_Continental_US_Polygon_Csph_Try2;
YellowStoneModelPolygonLoader;
NWUSModelPolygonLoader;

ME16veloAtDepth = ME16File.model(find(ME16File.depth == target_rad));
MITveloAtDepth = MITFile.model(find(MITFile.depth == target_rad));
USSLveloAtDepth = USSLFile.model(find(USSLFile.depth == target_rad));
YS_P_H15veloAtDepth = YS_P_H15File.model(find(YS_P_H15File.depth == target_rad));
NWUSveloAtDepth = NWUSFile.model(find(NWUSFile.depth == target_rad));
DNA13veloAtDepth = DNA13File.model(find(DNA13File.depth == target_rad));
GYPSUMveloAtDepth = GYPSUMFile.model(find(GYPSUMFile.depth == target_rad));
GAP_P4veloAtDepth = GAP_P4File.model(find(GAP_P4File.depth == target_rad));


%% Zero out the annoying stuff beyond model bounds for the regional models.

[in,on] = inpolygon( Grids.lon,Grids.lat,Continental_US_Lon,Continental_US_Lat);   % Logical Matrix
inon = in | on;     % Combine ?in? And ?on?
idx = find(inon(:));  
Geog_Indices = idx;

%%%%%%%%%%%%%%%%%%%%%% 
%Ensure model is zero outside bounds.
Empty_Depth_Vector = zeros(1,6*2^(2*8));
for i = 1:length(Geog_Indices)
    Empty_Depth_Vector(Geog_Indices(i)) = USSLveloAtDepth(Geog_Indices(i));
end
USSLveloAtDepth=Empty_Depth_Vector;
Empty_Depth_Vector = zeros(1,6*2^(2*8));
for i = 1:length(Geog_Indices)
    Empty_Depth_Vector(Geog_Indices(i)) = DNA13veloAtDepth(Geog_Indices(i));
end
DNA13veloAtDepth=Empty_Depth_Vector;

%%%%% Treat Yellowstone
[in,on] = inpolygon( Grids.lon,Grids.lat,YellowStoneModelLon,YellowStoneModelLat);   % Logical Matrix
inon = in | on;     % Combine ?in? And ?on?
idx = find(inon(:));  
Geog_Indices = idx;

Empty_Depth_Vector = zeros(1,6*2^(2*8));
for i = 1:length(Geog_Indices)
    Empty_Depth_Vector(Geog_Indices(i)) = YS_P_H15veloAtDepth(Geog_Indices(i));
end
YS_P_H15veloAtDepth=Empty_Depth_Vector;


%%%%% Treat NWUS
[in,on] = inpolygon( Grids.lon,Grids.lat,NWUSModelLon,NWUSModelLat);   % Logical Matrix
inon = in | on;     % Combine ?in? And ?on?
idx = find(inon(:));  
Geog_Indices = idx;

Empty_Depth_Vector = zeros(1,6*2^(2*8));
for i = 1:length(Geog_Indices)
    Empty_Depth_Vector(Geog_Indices(i)) = NWUSveloAtDepth(Geog_Indices(i));
end
NWUSveloAtDepth=Empty_Depth_Vector;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MIT_Face1 = MITveloAtDepth(find(Grids.face==1));
USSL_Face1 = USSLveloAtDepth(find(Grids.face==1));
ME16_Face1 = ME16veloAtDepth(find(Grids.face==1));
YS_P_H15_Face1 = YS_P_H15veloAtDepth(find(Grids.face==1));
DNA13_Face1 = DNA13veloAtDepth(find(Grids.face==1));
NWUS_Face1 = NWUSveloAtDepth(find(Grids.face==1));
GYPSUM_Face1 = GYPSUMveloAtDepth(find(Grids.face==1));
GAP_P4_Face1 = GAP_P4veloAtDepth(find(Grids.face==1));


MIT_Face1_reformat = reshape(MIT_Face1,[2^N 2^N]);
USSL_Face1_reformat = reshape(USSL_Face1,[2^N 2^N]);
ME16_Face1_reformat = reshape(ME16_Face1,[2^N 2^N]);
YS_P_H15_Face1_reformat = reshape(YS_P_H15_Face1,[2^N 2^N]);
DNA13_Face1_reformat = reshape(DNA13_Face1,[2^N 2^N]);
NWUS_Face1_reformat = reshape(NWUS_Face1,[2^N 2^N]);
GYPSUM_Face1_reformat = reshape(GYPSUM_Face1,[2^N 2^N]);
GAP_P4_Face1_reformat = reshape(GAP_P4_Face1,[2^N 2^N]);


[MIT_WaveletOut,MIT_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,MIT_Face1_reformat,Jmax);
[USSL_WaveletOut,USSL_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,USSL_Face1_reformat,Jmax);
[ME16_WaveletOut,ME16_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,ME16_Face1_reformat,Jmax);
[YS_P_H15_WaveletOut,YS_P_H15_PixelOut] = Threshold_Model_Geographic(N,Grids,YellowStoneModelLon,YellowStoneModelLat,YS_P_H15_Face1_reformat,Jmax);
[NWUS_WaveletOut,NWUS_PixelOut] = Threshold_Model_Geographic(N,Grids,NWUSModelLon,NWUSModelLat,NWUS_Face1_reformat,Jmax);
[DNA13_WaveletOut,DNA13_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,DNA13_Face1_reformat,Jmax);
[GYPSUM_WaveletOut,GYPSUM_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,GYPSUM_Face1_reformat,Jmax);
[GAP_P4_WaveletOut,GAP_P4_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,GAP_P4_Face1_reformat,Jmax);

vwlev1=cube2scale(N,[Jmax Jmax]+1,1,1);

for llll = 1:length(Scale_List)
llll
Desired_Scales = Scale_List(1:llll);

Reconstructed_MITWavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_ME16Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_USSLWavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_YS_P_H15Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_DNA13Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_NWUSWavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_GYPSUMWavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_GAP_P4Wavelet_Map_ScaleInspected = zeros(2^N,2^N);

for ii = 1:length(Desired_Scales)
    curr_scale = Desired_Scales(ii);
Scale_check = find(vwlev1 == curr_scale);
Reconstructed_MITWavelet_Map_ScaleInspected(Scale_check) = MIT_WaveletOut(Scale_check);
Reconstructed_ME16Wavelet_Map_ScaleInspected(Scale_check) = ME16_WaveletOut(Scale_check);
Reconstructed_USSLWavelet_Map_ScaleInspected(Scale_check) = USSL_WaveletOut(Scale_check);
Reconstructed_YS_P_H15Wavelet_Map_ScaleInspected(Scale_check) = YS_P_H15_WaveletOut(Scale_check);
Reconstructed_DNA13Wavelet_Map_ScaleInspected(Scale_check) = DNA13_WaveletOut(Scale_check);
Reconstructed_NWUSWavelet_Map_ScaleInspected(Scale_check) = NWUS_WaveletOut(Scale_check);
Reconstructed_GYPSUMWavelet_Map_ScaleInspected(Scale_check) = GYPSUM_WaveletOut(Scale_check);
Reconstructed_GAP_P4Wavelet_Map_ScaleInspected(Scale_check) = GAP_P4_WaveletOut(Scale_check);
end

%% Now do the inverse transform! 
%%
Reconstructed_MITVeloMap = angularD4WT(Reconstructed_MITWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_ME16VeloMap = angularD4WT(Reconstructed_ME16Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_USSLVeloMap = angularD4WT(Reconstructed_USSLWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_DNA13VeloMap = angularD4WT(Reconstructed_DNA13Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_NWUSVeloMap = angularD4WT(Reconstructed_NWUSWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_YS_P_H15VeloMap = angularD4WT(Reconstructed_YS_P_H15Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_GYPSUMVeloMap = angularD4WT(Reconstructed_GYPSUMWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_GAP_P4VeloMap = angularD4WT(Reconstructed_GAP_P4Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);


Curr_MITRMS = sqrt(mean(mean(Reconstructed_MITVeloMap.^2)));
Curr_ME16RMS = sqrt(mean(mean(Reconstructed_ME16VeloMap.^2)));
Curr_USSLRMS = sqrt(mean(mean(Reconstructed_USSLVeloMap.^2)));
Curr_DNA13RMS = sqrt(mean(mean(Reconstructed_DNA13VeloMap.^2)));
Curr_NWUSRMS = sqrt(mean(mean(Reconstructed_NWUSVeloMap.^2)));
Curr_YS_P_H15RMS = sqrt(mean(mean(Reconstructed_YS_P_H15VeloMap.^2)));
Curr_GYPSUMRMS = sqrt(mean(mean(Reconstructed_GYPSUMVeloMap.^2)));
Curr_GAP_P4RMS = sqrt(mean(mean(Reconstructed_GAP_P4VeloMap.^2)));



% RMS_Struct(llll).MIT(j) = [Curr_MITRMS];
% RMS_Struct(llll).ME16(j) = [Curr_ME16RMS];
% RMS_Struct(llll).USSL(j) = [Curr_USSLRMS];
MIT_CumRMS(llll)= Curr_MITRMS;
ME16_CumRMS(llll) = Curr_ME16RMS;
USSL_CumRMS(llll) = Curr_USSLRMS;
DNA13_CumRMS(llll) = Curr_DNA13RMS;
NWUS_CumRMS(llll) = Curr_NWUSRMS;
YS_P_H15_CumRMS(llll) = Curr_YS_P_H15RMS;
GYPSUM_CumRMS(llll) = Curr_GYPSUMRMS;
GAP_P4_CumRMS(llll) = Curr_GAP_P4RMS;

end
% figure()
% MITplot = plot(Scale_List,MIT_RMS,'linewidth',2);
% hold on
% ME16plot = plot(Scale_List,ME16_RMS,'linewidth',2);
% USSLplot = plot(Scale_List,USSL_RMS,'linewidth',2);
% DNA13plot = plot(Scale_List,DNA13_RMS,'linewidth',2);
% NWUSplot = plot(Scale_List,NWUS_RMS,'linewidth',2);
% %YS_P_H15plot = plot(Scale_List,YS_P_H15_RMS,'linewidth',2);
% GYPSUMplot = plot(Scale_List,GYPSUM_RMS,'linewidth',2);
% GAP_P4plot = plot(Scale_List,GAP_P4_RMS,'linewidth',2,'Color',[1,0.7,0.9]);
% 
% %title(['Scale-Based Model RMS at ' num2str(target_dep) ' km'],'fontsize',18)
% xlabel('Wavelet Scale','fontsize',12);  
% xticks([1 2 3 4 5 6])
% xticklabels({'1' '2' '3' '4' '5' '6'})
% 
% ylabel('RMS Amplitude','fontsize',12); 
% set(gca,'box','on','yscale','log'); 
% hLegend = legend( ...
%   [MITplot ME16plot USSLplot DNA13plot NWUSplot GYPSUMplot GAP_P4plot], ...
%   'MITP08' , 'ME16', 'USSL-2014','DNA13','NWUS','GYPSUM', 'GAP\_P4');
% grid on
% set(gca,'xdir','reverse') 
