%% Calculate inter-model correlation over a specific region. 

Prelude
Load_Continental_US_Polygon_Csph_Try2
re = 6371;
depth=600;
desired_radius = re - depth;
N=8;
Jmax=5;
Scale_List = [6 5 4 3 2 1]; 
LonMin = -125.6;
LonMax = -67;
LatMin = 24.4;
LatMax = 50.2;
Geog_Indices = find(Grids.lon > LonMin & Grids.lon < LonMax & Grids.lat > LatMin & Grids.lat<LatMax);
Data.Indices = Geog_Indices;

redo = 1;
if redo
Corr_MITME16_Cum_Vec = [];
Corr_MITUSSSL_Cum_Vec = [];
Corr_ME16USSSL_Cum_Vec = [];

Corr_MITME16_Single_Vec = [];
Corr_MITUSSSL_Single_Vec = [];
Corr_ME16USSSL_Single_Vec = [];


%MIT    
MITwvAtDepth = MITFile.wvcoeffs(find(MITFile.depth == desired_radius));
MITveloAtDepth = MITFile.model(find(MITFile.depth == desired_radius));
%ME16
ME16wvAtDepth = ME16File.wvcoeffs(find(ME16File.depth == desired_radius));
ME16veloAtDepth = ME16File.model(find(ME16File.depth == desired_radius));
%Schmandt
USSLwvAtDepth = USSLFile.wvcoeffs(find(USSLFile.depth == desired_radius));
USSLveloAtDepth = USSLFile.model(find(USSLFile.depth == desired_radius));

%Ensure model is zero outside bounds.
Empty_Depth_Vector = zeros(1,6*2^(2*8));
for i = 1:length(Data.Indices)
    Empty_Depth_Vector(Data.Indices(i)) = USSLveloAtDepth(Data.Indices(i));
end
USSLveloAtDepth=Empty_Depth_Vector;


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

Reconstructed_MITWavelet_Map_SingleScaleInspected = zeros(2^N,2^N);
Reconstructed_ME16Wavelet_Map_SingleScaleInspected = zeros(2^N,2^N);
Reconstructed_USSLWavelet_Map_SingleScaleInspected = zeros(2^N,2^N);



for ii = 1:length(Desired_Scales)
    curr_scale = Desired_Scales(ii);
    Scale_check = find(vwlev1 == curr_scale);
    Reconstructed_MITWavelet_Map_ScaleInspected(Scale_check) = MIT_WaveletOut(Scale_check);
    Reconstructed_ME16Wavelet_Map_ScaleInspected(Scale_check) = ME16_WaveletOut(Scale_check);
    Reconstructed_USSLWavelet_Map_ScaleInspected(Scale_check) = USSL_WaveletOut(Scale_check);
end

Scale_check = find(vwlev1 == Scale_List(llll));

Reconstructed_MITWavelet_Map_SingleScaleInspected(Scale_check) = MIT_WaveletOut(Scale_check);
Reconstructed_ME16Wavelet_Map_SingleScaleInspected(Scale_check) = ME16_WaveletOut(Scale_check);
Reconstructed_USSLWavelet_Map_SingleScaleInspected(Scale_check) = USSL_WaveletOut(Scale_check);


%% Now do the inverse transform! 
%%
Reconstructed_MITVeloMap = angularD4WT(Reconstructed_MITWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_ME16VeloMap = angularD4WT(Reconstructed_ME16Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_USSLVeloMap = angularD4WT(Reconstructed_USSLWavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);

Reconstructed_MITVeloMap_SingleScale = angularD4WT(Reconstructed_MITWavelet_Map_SingleScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_ME16VeloMap_SingleScale = angularD4WT(Reconstructed_ME16Wavelet_Map_SingleScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
Reconstructed_USSLVeloMap_SingleScale = angularD4WT(Reconstructed_USSLWavelet_Map_SingleScaleInspected,[Jmax Jmax],[1 1],'inverse',1);




% Now get model correlation in model domain
CumCorrMatMITME16 = corrcoef(Reconstructed_MITVeloMap,Reconstructed_ME16VeloMap);
CumCorrMITME16 = CumCorrMatMITME16(1,2);
CumCorrMatMITUSSL = corrcoef(Reconstructed_MITVeloMap,Reconstructed_USSLVeloMap);
CumCorrMITUSSL = CumCorrMatMITUSSL(1,2);
CumCorrMatME16USSL = corrcoef(Reconstructed_ME16VeloMap,Reconstructed_USSLVeloMap);
CumCorrME16USSL =CumCorrMatME16USSL(1,2);

SingleScaleCorrMatMITME16 = corrcoef(Reconstructed_MITVeloMap_SingleScale,Reconstructed_ME16VeloMap_SingleScale);
SingleScaleCorrMITME16=SingleScaleCorrMatMITME16(1,2);
SingleScaleCorrMatMITUSSL = corrcoef(Reconstructed_MITVeloMap_SingleScale,Reconstructed_USSLVeloMap_SingleScale);
SingleScaleCorrMITUSSL=SingleScaleCorrMatMITUSSL(1,2);
SingleScaleCorrMatME16USSL = corrcoef(Reconstructed_ME16VeloMap_SingleScale,Reconstructed_USSLVeloMap_SingleScale);
SingleScaleCorrME16USSL=SingleScaleCorrMatME16USSL(1,2);

Corr_MITME16_Cum_Vec = [Corr_MITME16_Cum_Vec CumCorrMITME16];
Corr_MITUSSSL_Cum_Vec = [Corr_MITUSSSL_Cum_Vec CumCorrMITUSSL];
Corr_ME16USSSL_Cum_Vec = [Corr_ME16USSSL_Cum_Vec CumCorrME16USSL];

Corr_MITME16_Single_Vec = [Corr_MITME16_Single_Vec SingleScaleCorrMITME16];
Corr_MITUSSSL_Single_Vec = [Corr_MITUSSSL_Single_Vec SingleScaleCorrMITUSSL];
Corr_ME16USSSL_Single_Vec = [Corr_ME16USSSL_Single_Vec SingleScaleCorrME16USSL];
end
 
else
%load
end

%MIT-ME16
figure()
tot=plot(Scale_List,Corr_MITME16_Cum_Vec,'-ro','LineWidth',4)
hold on
single=plot(Scale_List,Corr_MITME16_Single_Vec,'-bo','LineWidth',4)
grid on
box on
hXLabel=xlabel('Scale of Model')
hYLabel=ylabel('Model Correlation')
hTitle=title(['Model Domain Cumulative and Scale-Based Model Correlation between MITP08 and ME16\_Vp Extracted over the US at ' num2str(depth) 'km'])
hLegend = legend([tot,single],{'Cumulative Model Correlation','Correlation at Each Scale'})
xlim([1 6])
ylim([0 0.7])
set(gca,'xdir','reverse')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:0.1:1, ...
  'XTick'       , 1:1:5, ...
  'LineWidth'   , 2 , ...
  'FontSize', 23 );

set( hTitle                    , ...
    'FontSize'   , 27          , ...
    'FontWeight' , 'bold'      );

text(6.1,-0.02,'Scaling Fn','FontSize', 20)

%MIT-USSL
figure()
tot=plot(Scale_List,Corr_MITUSSSL_Cum_Vec,'-ro','LineWidth',4)
hold on
single=plot(Scale_List,Corr_MITUSSSL_Single_Vec,'-bo','LineWidth',4)
grid on
box on
hXLabel=xlabel('Scale of Model')
hYLabel=ylabel('Model Correlation')
hTitle=title([{'Model Domain Cumulative and Scale-Based Model Correlation' ' between  US-Sl-2014 and MITP08 over the US'}])
hLegend = legend([tot,single],{'Cumulative Model Correlation','Correlation at Each Scale'})
xlim([1 6])
ylim([0 0.7])
set(gca,'xdir','reverse')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:0.1:1, ...
  'XTick'       , 1:1:5, ...
  'LineWidth'   , 2 , ...
  'FontSize', 23 );

set( hTitle                    , ...
    'FontSize'   , 23          , ...
    'FontWeight' , 'bold'      );

text(6.1,-0.23,'Scaling Fn','FontSize', 20)
text(5,0.1,'Depth = 600 km','FontSize', 28)

%ME16-USSL
figure()
tot=plot(Scale_List,Corr_ME16USSSL_Cum_Vec,'-ro','LineWidth',4)
hold on
single=plot(Scale_List,Corr_ME16USSSL_Single_Vec,'-bo','LineWidth',4)
grid on
box on
hXLabel=xlabel('Scale of Model')
hYLabel=ylabel('Model Correlation')
hTitle=title([{'Model Domain Cumulative and Scale-Based Model Correlation' ' between US-Sl-2014 and ME16\_Vp'}])
hLegend = legend([tot,single],{'Cumulative Model Correlation','Correlation at Each Scale'})
xlim([1 6])
ylim([-0.2 0.2])
set(gca,'xdir','reverse')
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , -0.5:0.1:1, ...
  'XTick'       , 1:1:5, ...
  'LineWidth'   , 2 , ...
  'FontSize', 23 );

set( hTitle                    , ...
    'FontSize'   , 23          , ...
    'FontWeight' , 'bold'      );

text(6.1,-0.23,'Scaling Fn','FontSize', 20)
text(5,0.1,'Depth = 600 km','FontSize', 28)
