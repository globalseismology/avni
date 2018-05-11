%% Make RMS Plot with Depth for a desired Geographic Region. 

% Anant Hariharan
clear all
clc
close all

redo=0;
if redo
    
Prelude

%% Set Params
% Note, Africa is on face 3. 
% Note, US is on face 1??
N=8;
Jmax=5;
Load_Continental_US_Polygon_Csph_Try2

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


%% Now, get the RMS!
Curr_MITRMS = sqrt(mean(mean(MIT_PixelOut.^2)));
Curr_ME16RMS = sqrt(mean(mean(ME16_PixelOut.^2)));
Curr_USSLRMS = sqrt(mean(mean(USSL_PixelOut.^2)));
MIT_RMS = [MIT_RMS Curr_MITRMS];
ME16_RMS = [ME16_RMS Curr_ME16RMS];
USSL_RMS = [USSL_RMS Curr_USSLRMS];
end
else
load('MIT_RMS.mat')
load('ME16_RMS.mat')
load('USSL_RMS.mat')
load('ModelDepthsFlip.mat')
end   
    
    
    
figure()
MIT=plot(MIT_RMS,6371-ModelDepthsFlip,'-ro','LineWidth',6);
hold on
ME16=plot(ME16_RMS,6371-ModelDepthsFlip,'-bo','LineWidth',6);
USSL=plot(USSL_RMS,6371-ModelDepthsFlip,'-ko','LineWidth',6);
grid on
box on
hXLabel=xlabel('Tomographic Model RMS')
hYLabel=ylabel('Depth (km)')
set(gca,'ydir','reverse')
hTitle=title('Model RMS With Depth Over the US')
p=patch([0 1 1 0],[410 410 660 660],'r','FaceAlpha',0.2);
hLegend = legend([MIT,ME16,USSL],{'MIT Model','ME16 Model','US-SL 2014 Model'})
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

save('MIT_RMS.mat','MIT_RMS')
save('ME16_RMS.mat','ME16_RMS')
save('USSL_RMS.mat','USSL_RMS')
save('ModelDepthsFlip.mat','ModelDepthsFlip')



