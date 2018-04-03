function [RMS] = D4wt_Localized_Model_RMS(Model1,Geographic_Lon,Geographic_Lat,N,Jmax)
% Calculates RMS of a tomographic model on the cubed sphere, over a
% specific geographic region...
% (2^N * 2^N) arrays.
Scale_List = Jmax+1:-1:1;
Grids = load('/home/anant/Software/rem3d/rem3d/files/Csph_Grid_N8_Jmax5_EulerConfig1.mat');
CumCorrList =[];
% Locations of wavelet, i.e. scaling functions and such. 
vwlev1=cube2scale(N,[Jmax Jmax]+1,1,1);

% First threshold over the geographic region. 
[Model1_WaveletOut,Model1_PixelOut] = Threshold_Model_Geographic(N,Grids,Geographic_Lon,Geographic_Lat,Model1,Jmax);


% for llll = 1:length(Scale_List)
%     Desired_Scales = Scale_List(1:llll)
%     Reconstructed_Model1Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
%     Reconstructed_Model2Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
%         
%     for ii = 1:length(Desired_Scales)
%     curr_scale = Desired_Scales(ii);
%     Scale_check = find(vwlev1 == curr_scale);
%     Reconstructed_Model1Wavelet_Map_ScaleInspected(Scale_check) = Model1_WaveletOut(Scale_check);
%     Reconstructed_Model2Wavelet_Map_ScaleInspected(Scale_check) = Model2_WaveletOut(Scale_check);
%     end
% 
%     Reconstructed_Model1VeloMap = angularD4WT(Reconstructed_Model1Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
%     Reconstructed_Model2VeloMap = angularD4WT(Reconstructed_Model2Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
% 
%     
    % Now only extract it over the desired geographic region.
   
    [in,on] = inpolygon( Grids.lon,Grids.lat,Geographic_Lon,Geographic_Lat);   % Logical Matrix
    inon = in | on;     % Combine ?in? And ?on?
    idx = find(inon(:));  
    Geog_Indices = idx;

    Model1_torms = [];
    
    for i = 1:length(Geog_Indices)
        Model1_torms(i) = Model1_PixelOut(Geog_Indices(i));
    end
   
    %%%
    RMS = sqrt(mean(Model1_torms.^2));
    
    
    



end

