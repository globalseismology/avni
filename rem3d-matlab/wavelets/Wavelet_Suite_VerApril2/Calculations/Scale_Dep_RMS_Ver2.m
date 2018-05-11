%% Make Scale-Dependent RMS Plot with Depth for a desired Geographic Region. 

% Anant Hariharan
clear all
clc
close all

redo=1;
if redo
    
Prelude

%% Set Params
% Note, Africa is on face 3. 
% Note, US is on face 1??
N=8;
Jmax=5;
Load_Continental_US_Polygon_Csph_Try2
Scale_List = [6 5 4 3 2 1]; 
LonMin = -125.6;
LonMax = -67;
LatMin = 24.4;
LatMax = 50.2;
Geog_Indices = find(Grids.lon > LonMin & Grids.lon < LonMax & Grids.lat > LatMin & Grids.lat<LatMax);
Data.Indices = Geog_Indices;

re = 6371;
% depth=100;
% desired_radius = re - depth;
MITModelDepths = unique(MITFile.depth);
USSLModelDepths = unique(USSLFile.depth);
ModelDepthsFlip = flipud(USSLModelDepths);

MIT_RMS = [];
ME16_RMS = [];
USSL_RMS = [];


for j = 1:length(ModelDepthsFlip)
Curr_Depth = ModelDepthsFlip(j)
%% Now, find the model coefficients that correspond to the points you desire. Put in a structure. 
% First, extract the coefficients for the right depth. 
%MIT
MITwvAtDepth = MITFile.wvcoeffs(find(MITFile.depth == Curr_Depth));
MITveloAtDepth = MITFile.model(find(MITFile.depth == Curr_Depth));
%ME16
ME16wvAtDepth = ME16File.wvcoeffs(find(ME16File.depth == Curr_Depth));
ME16veloAtDepth = ME16File.model(find(ME16File.depth == Curr_Depth));
%Schmandt
USSLwvAtDepth = USSLFile.wvcoeffs(find(USSLFile.depth == Curr_Depth));
USSLveloAtDepth = USSLFile.model(find(USSLFile.depth == Curr_Depth));

%Ensure model is zero outside bounds.
Empty_Depth_Vector = zeros(1,6*2^(2*8));
for i = 1:length(Data.Indices)
    Empty_Depth_Vector(Data.Indices(i)) = USSLveloAtDepth(Data.Indices(i));
end
USSLveloAtDepth=Empty_Depth_Vector;

MIT_Face1 = MITveloAtDepth(find(Grids.face==1));
USSL_Face1 = USSLveloAtDepth(find(Grids.face==1));
ME16_Face1 = ME16veloAtDepth(find(Grids.face==1));

MIT_Face1_reformat = reshape(MIT_Face1,[2^N 2^N]);
USSL_Face1_reformat = reshape(USSL_Face1,[2^N 2^N]);
ME16_Face1_reformat = reshape(ME16_Face1,[2^N 2^N]);

[MIT_WaveletOut,MIT_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,MIT_Face1_reformat,Jmax);
[USSL_WaveletOut,USSL_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,USSL_Face1_reformat,Jmax);
[ME16_WaveletOut,ME16_PixelOut] = Threshold_Model_Geographic(N,Grids,Continental_US_Lon,Continental_US_Lat,ME16_Face1_reformat,Jmax);




vwlev1=cube2scale(N,[Jmax Jmax]+1,1,1);
for llll = 1:length(Scale_List)

Desired_Scales = Scale_List(1:llll)
Reconstructed_MITWavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_ME16Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
Reconstructed_USSLWavelet_Map_ScaleInspected = zeros(2^N,2^N);

for ii = 1:length(Desired_Scales)
    curr_scale = Desired_Scales(ii);
    Scale_check = find(vwlev1 == curr_scale);
    Reconstructed_MITWavelet_Map_ScaleInspected(Scale_check) = MIT_WaveletOut(Scale_check);
    Reconstructed_ME16Wavelet_Map_ScaleInspected(Scale_check) = ME16_WaveletOut(Scale_check);
    Reconstructed_USSLWavelet_Map_ScaleInspected(Scale_check) = USSL_WaveletOut(Scale_check);
end


%% Now do the inverse transform! 
%%
Reconstructed_MITVeloMap = angularD4WT(Reconstructed_MITWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_ME16VeloMap = angularD4WT(Reconstructed_ME16Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_USSLVeloMap = angularD4WT(Reconstructed_USSLWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);

Curr_MITRMS = sqrt(mean(mean(Reconstructed_MITVeloMap.^2)));
Curr_ME16RMS = sqrt(mean(mean(Reconstructed_ME16VeloMap.^2)));
Curr_USSLRMS = sqrt(mean(mean(Reconstructed_USSLVeloMap.^2)));

RMS_Struct(llll).MIT(j) = [Curr_MITRMS];
RMS_Struct(llll).ME16(j) = [Curr_ME16RMS];
RMS_Struct(llll).USSL(j) = [Curr_USSLRMS];

end



end

else
    
    
load('ScaleBasedRMS_Struct.mat')



end   
    
 
figure()    

Scale6=plot(RMS_Struct(1).MIT,6371-ModelDepthsFlip,'-ro','LineWidth',4);
hold on
Scale5=plot(RMS_Struct(2).MIT,6371-ModelDepthsFlip,'-bo','LineWidth',4);
Scale4=plot(RMS_Struct(3).MIT,6371-ModelDepthsFlip,'-go','LineWidth',4);
Scale3=plot(RMS_Struct(4).MIT,6371-ModelDepthsFlip,'-co','LineWidth',4);
Scale2=plot(RMS_Struct(5).MIT,6371-ModelDepthsFlip,'-mo','LineWidth',4);
Scale1=plot(RMS_Struct(6).MIT,6371-ModelDepthsFlip,'-ko','LineWidth',4);

grid on
box on
hXLabel=xlabel('Tomographic Model RMS')
hYLabel=ylabel('Depth (km)')
set(gca,'ydir','reverse')
hTitle=title('Model RMS With Depth Over the US for MITP08')
p=patch([0 1 1 0],[410 410 660 660],'r','FaceAlpha',0.2);
hLegend = legend([Scale6,Scale5,Scale4,Scale3,Scale2,Scale1],{'Scaling Function','Up to Scale 5','Up to Scale 4','Up to Scale 3','Up to Scale 2','Up to Scale 1'})
xlim([0.08 0.25])
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
  'LineWidth'   , 2 , ...
  'FontSize', 23 );

set( hTitle                    , ...
    'FontSize'   , 27          , ...
    'FontWeight' , 'bold'      );


save('ScaleBasedRMS_Struct.mat','RMS_Struct')




figure()    

Scale6=plot(RMS_Struct(1).ME16,6371-ModelDepthsFlip,'-ro','LineWidth',4);
hold on
Scale5=plot(RMS_Struct(2).ME16,6371-ModelDepthsFlip,'-bo','LineWidth',4);
Scale4=plot(RMS_Struct(3).ME16,6371-ModelDepthsFlip,'-go','LineWidth',4);
Scale3=plot(RMS_Struct(4).ME16,6371-ModelDepthsFlip,'-co','LineWidth',4);
Scale2=plot(RMS_Struct(5).ME16,6371-ModelDepthsFlip,'-mo','LineWidth',4);
Scale1=plot(RMS_Struct(6).ME16,6371-ModelDepthsFlip,'-ko','LineWidth',4);

grid on
box on
hXLabel=xlabel('Tomographic Model RMS')
hYLabel=ylabel('Depth (km)')
set(gca,'ydir','reverse')
hTitle=title('Model RMS With Depth Over the US for ME16')
p=patch([0 1 1 0],[410 410 660 660],'r','FaceAlpha',0.2);
hLegend = legend([Scale6,Scale5,Scale4,Scale3,Scale2,Scale1],{'Scaling Function','Up to Scale 5','Up to Scale 4','Up to Scale 3','Up to Scale 2','Up to Scale 1'})
xlim([0.15 0.3])
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
  'LineWidth'   , 2 , ...
  'FontSize', 23 );

set( hTitle                    , ...
    'FontSize'   , 27          , ...
    'FontWeight' , 'bold'      );



figure()    

Scale6=plot(RMS_Struct(1).USSL,6371-ModelDepthsFlip,'-ro','LineWidth',4);
hold on
Scale5=plot(RMS_Struct(2).USSL,6371-ModelDepthsFlip,'-bo','LineWidth',4);
Scale4=plot(RMS_Struct(3).USSL,6371-ModelDepthsFlip,'-go','LineWidth',4);
Scale3=plot(RMS_Struct(4).USSL,6371-ModelDepthsFlip,'-co','LineWidth',4);
Scale2=plot(RMS_Struct(5).USSL,6371-ModelDepthsFlip,'-mo','LineWidth',4);
Scale1=plot(RMS_Struct(6).USSL,6371-ModelDepthsFlip,'-ko','LineWidth',4);

grid on
box on
hXLabel=xlabel('Tomographic Model RMS')
hYLabel=ylabel('Depth (km)')
set(gca,'ydir','reverse')
hTitle=title('Model RMS With Depth Over the US for USSL-2014')
p=patch([0 1 1 0],[410 410 660 660],'r','FaceAlpha',0.2);
hLegend = legend([Scale6,Scale5,Scale4,Scale3,Scale2,Scale1],{'Scaling Function','Up to Scale 5','Up to Scale 4','Up to Scale 3','Up to Scale 2','Up to Scale 1'})
xlim([0 0.5])
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
  'LineWidth'   , 2 , ...
  'FontSize', 23 );

set( hTitle                    , ...
    'FontSize'   , 27          , ...
    'FontWeight' , 'bold'      );
