%Pretty much the culmination of everything we've done so far involving
%wavelets

addpath(genpath('/home/moulik/Software/fjsimons-MMASS-v1.0.1'))
inferno; %load inferno colormap

clear; close all;


% Select the depth at which to do comparisons
target_dep = 100; 
target_rad = 6371-target_dep;

% To what wavelet scale should I analyze?
Jmax = 5; %This isn't really an option!

N=8;
% Load in a whole bunch of stuff
Prelude

Scale_List = [6 5 4 3 2 1]; 
ModelDepths = unique(USSLFile.depth);
ModelDepthsFlip = flipud(ModelDepths);

% load polygons...
Load_Continental_US_Polygon_Csph_Try2
%YellowStoneModelPolygonLoader
NWUSModelPolygonLoader

ME16veloAtDepth = ME16File.model(find(ME16File.depth == target_rad));
MITveloAtDepth = MITFile.model(find(MITFile.depth == target_rad));
USSLveloAtDepth = USSLFile.model(find(USSLFile.depth == target_rad));
%YS_P_H15veloAtDepth = YS_P_H15File.model(find(YS_P_H15File.depth == target_rad));
NWUSveloAtDepth = NWUSFile.model(find(NWUSFile.depth == target_rad));
DNA13veloAtDepth = DNA13File.model(find(DNA13File.depth == target_rad));
GYPSUMveloAtDepth = GYPSUMFile.model(find(GYPSUMFile.depth == target_rad));
GAP_P4veloAtDepth = GAP_P4File.model(find(GAP_P4File.depth == target_rad));

THBUSveloAtDepth = THBUSFile.model(find(THBUSFile.depth == target_rad));


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


Empty_Depth_Vector = zeros(1,6*2^(2*8));
for i = 1:length(Geog_Indices)
    Empty_Depth_Vector(Geog_Indices(i)) = THBUSveloAtDepth(Geog_Indices(i));
end
THBUSveloAtDepth=Empty_Depth_Vector;


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

%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MIT_Face1 = MITveloAtDepth(find(Grids.face==1));
USSL_Face1 = USSLveloAtDepth(find(Grids.face==1));
ME16_Face1 = ME16veloAtDepth(find(Grids.face==1));
THBUS_Face1 = THBUSveloAtDepth(find(Grids.face==1));
DNA13_Face1 = DNA13veloAtDepth(find(Grids.face==1));
NWUS_Face1 = NWUSveloAtDepth(find(Grids.face==1));
GYPSUM_Face1 = GYPSUMveloAtDepth(find(Grids.face==1));
GAP_P4_Face1 = GAP_P4veloAtDepth(find(Grids.face==1));


MIT_Face1_reformat = reshape(MIT_Face1,[2^N 2^N]);
USSL_Face1_reformat = reshape(USSL_Face1,[2^N 2^N]);
ME16_Face1_reformat = reshape(ME16_Face1,[2^N 2^N]);
THBUS_Face1_reformat = reshape(THBUS_Face1,[2^N 2^N]);
DNA13_Face1_reformat = reshape(DNA13_Face1,[2^N 2^N]);
NWUS_Face1_reformat = reshape(NWUS_Face1,[2^N 2^N]);
GYPSUM_Face1_reformat = reshape(GYPSUM_Face1,[2^N 2^N]);
GAP_P4_Face1_reformat = reshape(GAP_P4_Face1,[2^N 2^N]);


%%%%%% Start doing the comparisons over overlapping regions.
%First, ME16-ME16- do comparison over US.

%%Choose the main model!
MainModel = ME16_Face1_reformat;

[ME16_MIT_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MainModel,MIT_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_GYPSUM_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MainModel,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[ME16_USSL_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MainModel,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_DNA13_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MainModel,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_NWUS_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MainModel,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[ME16_GAP_P4_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MainModel,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_GYPSUM_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MIT_Face1_reformat,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[MITP08_USSL_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MIT_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_DNA13_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MIT_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_NWUS_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MIT_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[MITP08_GAP_P4_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MIT_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_USSL_Corr] = D4wt_Localized_Cumulative_Model_Correlation(GYPSUM_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_DNA13_Corr] = D4wt_Localized_Cumulative_Model_Correlation(GYPSUM_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_NWUS_Corr] = D4wt_Localized_Cumulative_Model_Correlation(GYPSUM_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[GYPSUM_GAP_P4_Corr] = D4wt_Localized_Cumulative_Model_Correlation(GYPSUM_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[USSL_DNA13_Corr] = D4wt_Localized_Cumulative_Model_Correlation(USSL_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[USSL_NWUS_Corr] = D4wt_Localized_Cumulative_Model_Correlation(USSL_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[USSL_GAP_P4_Corr] = D4wt_Localized_Cumulative_Model_Correlation(USSL_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[DNA13_NWUS_Corr] = D4wt_Localized_Cumulative_Model_Correlation(DNA13_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[DNA13_GAP_P4_Corr] = D4wt_Localized_Cumulative_Model_Correlation(DNA13_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[NWUS_GAP_P4_Corr] = D4wt_Localized_Cumulative_Model_Correlation(NWUS_Face1_reformat,GAP_P4_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
%
[THBUS_ME16_Corr] = D4wt_Localized_Cumulative_Model_Correlation(MainModel,THBUS_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_MIT_Corr] = D4wt_Localized_Cumulative_Model_Correlation(THBUS_Face1_reformat,MIT_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_GYPSUM_Corr] = D4wt_Localized_Cumulative_Model_Correlation(THBUS_Face1_reformat,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[THBUS_USSL_Corr] = D4wt_Localized_Cumulative_Model_Correlation(THBUS_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_DNA13_Corr] = D4wt_Localized_Cumulative_Model_Correlation(THBUS_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_NWUS_Corr] = D4wt_Localized_Cumulative_Model_Correlation(THBUS_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[THBUS_GAP_P4_Corr] = D4wt_Localized_Cumulative_Model_Correlation(THBUS_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);




BigMat = [ME16_MIT_Corr; ME16_GYPSUM_Corr; ME16_USSL_Corr; ME16_DNA13_Corr; ME16_NWUS_Corr; ME16_GAP_P4_Corr; MITP08_GYPSUM_Corr;...
    MITP08_USSL_Corr; MITP08_DNA13_Corr; MITP08_NWUS_Corr; MITP08_GAP_P4_Corr; GYPSUM_USSL_Corr; GYPSUM_DNA13_Corr; GYPSUM_NWUS_Corr;...
     GYPSUM_GAP_P4_Corr; USSL_DNA13_Corr; USSL_NWUS_Corr; USSL_GAP_P4_Corr;...
    DNA13_NWUS_Corr; DNA13_GAP_P4_Corr; NWUS_GAP_P4_Corr; THBUS_ME16_Corr; THBUS_MIT_Corr; ...
    THBUS_GYPSUM_Corr; THBUS_USSL_Corr; THBUS_DNA13_Corr; THBUS_NWUS_Corr; THBUS_GAP_P4_Corr];

%%%% Now do Scale-By-Scale

[ME16_MIT_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MainModel,MIT_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_GYPSUM_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MainModel,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[ME16_USSL_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MainModel,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_DNA13_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MainModel,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_NWUS_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MainModel,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[ME16_GAP_P4_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MainModel,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_GYPSUM_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MIT_Face1_reformat,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[MITP08_USSL_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MIT_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_DNA13_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MIT_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_NWUS_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MIT_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[MITP08_GAP_P4_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MIT_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_USSL_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(GYPSUM_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_DNA13_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(GYPSUM_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_NWUS_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(GYPSUM_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[GYPSUM_GAP_P4_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(GYPSUM_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[USSL_DNA13_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(USSL_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[USSL_NWUS_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(USSL_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[USSL_GAP_P4_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(USSL_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[DNA13_NWUS_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(DNA13_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[DNA13_GAP_P4_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(DNA13_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[NWUS_GAP_P4_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(NWUS_Face1_reformat,GAP_P4_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
%
[THBUS_ME16_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(MainModel,THBUS_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_MIT_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(THBUS_Face1_reformat,MIT_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_GYPSUM_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(THBUS_Face1_reformat,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[THBUS_USSL_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(THBUS_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_DNA13_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(THBUS_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_NWUS_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(THBUS_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[THBUS_GAP_P4_ScaleByScaleCorr] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(THBUS_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);




BigScaleByScaleMat = [ME16_MIT_ScaleByScaleCorr; ME16_GYPSUM_ScaleByScaleCorr; ME16_USSL_ScaleByScaleCorr; ME16_DNA13_ScaleByScaleCorr; ME16_NWUS_ScaleByScaleCorr; ME16_GAP_P4_ScaleByScaleCorr; MITP08_GYPSUM_ScaleByScaleCorr;...
    MITP08_USSL_ScaleByScaleCorr; MITP08_DNA13_ScaleByScaleCorr; MITP08_NWUS_ScaleByScaleCorr; MITP08_GAP_P4_ScaleByScaleCorr; GYPSUM_USSL_ScaleByScaleCorr; GYPSUM_DNA13_ScaleByScaleCorr; GYPSUM_NWUS_ScaleByScaleCorr;...
     GYPSUM_GAP_P4_ScaleByScaleCorr; USSL_DNA13_ScaleByScaleCorr; USSL_NWUS_ScaleByScaleCorr; USSL_GAP_P4_ScaleByScaleCorr;...
    DNA13_NWUS_ScaleByScaleCorr; DNA13_GAP_P4_ScaleByScaleCorr; NWUS_GAP_P4_ScaleByScaleCorr; THBUS_ME16_ScaleByScaleCorr; THBUS_MIT_ScaleByScaleCorr; ...
    THBUS_GYPSUM_ScaleByScaleCorr; THBUS_USSL_ScaleByScaleCorr; THBUS_DNA13_ScaleByScaleCorr; THBUS_NWUS_ScaleByScaleCorr; THBUS_GAP_P4_ScaleByScaleCorr];


Model_RMS_Ratios


[MIT_RMS] = Compute_NonCumulative_ModelRMS_FunctionVer(MIT_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[ME16_RMS] = Compute_NonCumulative_ModelRMS_FunctionVer(MainModel,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[USSL_RMS] = Compute_NonCumulative_ModelRMS_FunctionVer(USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[DNA13_RMS] = Compute_NonCumulative_ModelRMS_FunctionVer(DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[NWUS_RMS] = Compute_NonCumulative_ModelRMS_FunctionVer(NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax)
[GYPSUM_RMS] = Compute_NonCumulative_ModelRMS_FunctionVer(GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[GAP_P4_RMS] = Compute_NonCumulative_ModelRMS_FunctionVer(GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[THBUS_RMS] = Compute_NonCumulative_ModelRMS_FunctionVer(THBUS_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)



[MIT_CumRMS] = Compute_Cumulative_ModelRMS_FunctionVer(MIT_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[ME16_CumRMS] = Compute_Cumulative_ModelRMS_FunctionVer(MainModel,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[USSL_CumRMS] = Compute_Cumulative_ModelRMS_FunctionVer(USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[DNA13_CumRMS] = Compute_Cumulative_ModelRMS_FunctionVer(DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[NWUS_CumRMS] = Compute_Cumulative_ModelRMS_FunctionVer(NWUS_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[GYPSUM_CumRMS] = Compute_Cumulative_ModelRMS_FunctionVer(GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[GAP_P4_CumRMS] = Compute_Cumulative_ModelRMS_FunctionVer(GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[THBUS_CumRMS] = Compute_Cumulative_ModelRMS_FunctionVer(THBUS_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)

YMatrix1 = [MIT_RMS; ME16_RMS; USSL_RMS; DNA13_RMS; NWUS_RMS; GYPSUM_RMS; GAP_P4_RMS; THBUS_RMS;];
YMatrix2 = [MIT_CumRMS; ME16_CumRMS; USSL_CumRMS; DNA13_CumRMS; NWUS_CumRMS; GYPSUM_CumRMS; GAP_P4_CumRMS; THBUS_CumRMS;];

createfigure1(BigMat, Frac_RMS_List, Scale_List, YMatrix2,'Cumulative Correlation')

createfigure1(BigScaleByScaleMat, Frac_RMS_List, Scale_List, YMatrix1,'Scale-By-Scale Correlation')
