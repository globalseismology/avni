function [CorrList] = D4wt_ScaleByScaleLocalized_Cumulative_Model_Correlation(Model1,Model2,Geographic_Lon,Geographic_Lat,N,Jmax)
% Calculates Correlation of Different Tomographic Models over a specified
% Geographic Region. These models are provided as Faces on a cubed sphere.
% i.e. (2^N * 2^N) arrays.
Scale_List = Jmax+1:-1:1;
Grids = load('/home/anant/Software/rem3d/rem3d/files/Csph_Grid_N8_Jmax5_EulerConfig1.mat');
CorrList =[];
% Locations of wavelet, i.e. scaling functions and such. 
vwlev1=cube2scale(N,[Jmax Jmax]+1,1,1);

% First threshold over the geographic region. 
[Model1_WaveletOut,Model1_PixelOut] = Threshold_Model_Geographic(N,Grids,Geographic_Lon,Geographic_Lat,Model1,Jmax);
[Model2_WaveletOut,Model2_PixelOut] = Threshold_Model_Geographic(N,Grids,Geographic_Lon,Geographic_Lat,Model2,Jmax);


for llll = 1:length(Scale_List)
    Desired_Scales = Scale_List(llll)
    Reconstructed_Model1Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
    Reconstructed_Model2Wavelet_Map_ScaleInspected = zeros(2^N,2^N);
        
    
    curr_scale = Desired_Scales;
    Scale_check = find(vwlev1 == curr_scale);
    Reconstructed_Model1Wavelet_Map_ScaleInspected(Scale_check) = Model1_WaveletOut(Scale_check);
    Reconstructed_Model2Wavelet_Map_ScaleInspected(Scale_check) = Model2_WaveletOut(Scale_check);
    

    Reconstructed_Model1VeloMap = angularD4WT(Reconstructed_Model1Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);
    Reconstructed_Model2VeloMap = angularD4WT(Reconstructed_Model2Wavelet_Map_ScaleInspected,[Jmax Jmax],[1 1],'inverse',1);

    
    % Now only extract it over the desired geographic region.
   
    [in,on] = inpolygon( Grids.lon,Grids.lat,Geographic_Lon,Geographic_Lat);   % Logical Matrix
    inon = in | on;     % Combine ?in? And ?on?
    idx = find(inon(:));  
    Geog_Indices = idx;

    Model1_tocorr = [];
    Model2_tocorr = [];
    for i = 1:length(Geog_Indices)
        Model1_tocorr(i) = Reconstructed_Model1VeloMap(Geog_Indices(i));
        Model2_tocorr(i) = Reconstructed_Model2VeloMap(Geog_Indices(i));
    end
   
    %%%
    CumCorr = corrcoef(Model1_tocorr,Model2_tocorr);
    CorrList(llll) = CumCorr(1,2);
    
    
    
end


end

