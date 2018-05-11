%% Calculate model correlation across a depth. 

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

re = 6371;
depth=2800;
desired_radius = re - depth;

Corr_MITME16_Cum_Vec = [];
Corr_MITUSSSL_Cum_Vec = [];
Corr_ME16USSSL_Cum_Vec = [];

Corr_MITME16_Single_Vec = [];
Corr_MITUSSSL_Single_Vec = [];
Corr_ME16USSSL_Single_Vec = [];

N=8;
Jmax=5;
%MIT    
MITwvAtDepth = MITFile.wvcoeffs(find(MITFile.depth == desired_radius));
MITveloAtDepth = MITFile.model(find(MITFile.depth == desired_radius));
%ME16
ME16wvAtDepth = ME16File.wvcoeffs(find(ME16File.depth == desired_radius));
ME16veloAtDepth = ME16File.model(find(ME16File.depth == desired_radius));


MITModelCsph= reshape(MITveloAtDepth,[2^8 2^8 6]);
MITWaveletsCsph= reshape(MITwvAtDepth,[2^8 2^8 6]);

ME16ModelCsph= reshape(ME16veloAtDepth,[2^8 2^8 6]);
ME16WaveletsCsph= reshape(ME16wvAtDepth,[2^8 2^8 6]);
Scale_List = [6 5 4 3 2 1]; 
[vwlev1,vwlevs1]=cube2scale(N,[Jmax Jmax]+1,1,1);  



for iiii = 1:length(Scale_List)
iiii

Desired_Scales = Scale_List(1:iiii);
SingleScale = Scale_List(iiii);

Reconstructed_MITWavelet_Map_ScaleInspected = zeros(2^N,2^N,6);
Reconstructed_ME16Wavelet_Map_ScaleInspected = zeros(2^N,2^N,6);

for ii = 1:length(Desired_Scales)
    curr_scale = Desired_Scales(ii);
    Scale_check = find(vwlevs1 == curr_scale);
    Reconstructed_MITWavelet_Map_ScaleInspected(Scale_check) = MITWaveletsCsph(Scale_check);
    Reconstructed_ME16Wavelet_Map_ScaleInspected(Scale_check) = ME16WaveletsCsph(Scale_check);
end

Reconstructed_MITWavelet_Map_SingleScaleInspected = zeros(2^N,2^N,6);
Reconstructed_ME16Wavelet_Map_SingleScaleInspected = zeros(2^N,2^N,6);

curr_scale = SingleScale;
Scale_check = find(vwlevs1 == curr_scale);
Reconstructed_MITWavelet_Map_SingleScaleInspected(Scale_check) = MITWaveletsCsph(Scale_check);
Reconstructed_ME16Wavelet_Map_SingleScaleInspected(Scale_check) = ME16WaveletsCsph(Scale_check);


%% Now do the inverse transform for cumulative scales! 
%%
Reconstructed_MITVeloMap = zeros(2^N,2^N);
Reconstructed_ME16VeloMap = zeros(2^N,2^N);

for i = 1:2^(2*N)
    CurrMITMap = Reconstructed_MITWavelet_Map_ScaleInspected(i)*full(WaveletMaps_DV_Basis.Me(i).dv);
    CurrME16Map =  Reconstructed_ME16Wavelet_Map_ScaleInspected(i)*full(WaveletMaps_DV_Basis.Me(i).dv);

    Reconstructed_MITVeloMap = Reconstructed_MITVeloMap+CurrMITMap;
    Reconstructed_ME16VeloMap = Reconstructed_ME16VeloMap+CurrME16Map;
end



%% Now do the inverse transform for single scales! 
%%
Reconstructed_MITSingleScaleVeloMap = zeros(2^N,2^N);
Reconstructed_ME16SingleScaleVeloMap = zeros(2^N,2^N);

for i = 1:2^(2*N)
    CurrMITMap = Reconstructed_MITWavelet_Map_SingleScaleInspected(i)*full(WaveletMaps_DV_Basis.Me(i).dv);
    CurrME16Map = Reconstructed_ME16Wavelet_Map_SingleScaleInspected(i)*full(WaveletMaps_DV_Basis.Me(i).dv);

    Reconstructed_MITSingleScaleVeloMap = Reconstructed_MITSingleScaleVeloMap+CurrMITMap;
    Reconstructed_ME16SingleScaleVeloMap = Reconstructed_ME16SingleScaleVeloMap+CurrME16Map;
end

% Now get model correlation in model domain
CumCorrMatMITME16 = corrcoef(Reconstructed_MITVeloMap,Reconstructed_ME16VeloMap);
CumCorrMITME16 = CumCorrMatMITME16(1,2);


SingleScaleCorrMatMITME16 = corrcoef(Reconstructed_MITSingleScaleVeloMap,Reconstructed_ME16SingleScaleVeloMap);
SingleScaleCorrMITME16=SingleScaleCorrMatMITME16(1,2);


Corr_MITME16_Cum_Vec = [Corr_MITME16_Cum_Vec CumCorrMITME16];
Corr_MITME16_Single_Vec = [Corr_MITME16_Single_Vec SingleScaleCorrMITME16];



end





%MIT-ME16
figure()
tot=plot(Scale_List,Corr_MITME16_Cum_Vec,'-ro','LineWidth',6)
hold on
single=plot(Scale_List,Corr_MITME16_Single_Vec,'-bo','LineWidth',6)
grid on
box on
hXLabel=xlabel('Scale of Model')
hYLabel=ylabel('Model Correlation')
hTitle=title(['Global Cumulative and Scale-Based Model Correlation between MITP08 and ME16\_Vp at ' num2str(depth) 'km'])
hLegend = legend([tot,single],{'Cumulative Model Correlation','Correlation at Each Scale'})
xlim([1 6])
ylim([-0.05 0.7])
Scaltext = text(6.1,-0.09,'Scaling Fn');
Scaltext.FontSize = 23;
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
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
hLegend.FontSize = 23;
