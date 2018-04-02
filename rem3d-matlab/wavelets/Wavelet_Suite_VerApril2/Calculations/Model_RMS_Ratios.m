%% Make Bar Plot of RMS Ratio
%Get_Fractional_RMS(Model1,Model2,Geographic_Lon,Geographic_Lat,N,Jmax)

[ME16_MIT_RMSRatio] = Get_Fractional_RMS(MainModel,MIT_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_GYPSUM_RMSRatio] = Get_Fractional_RMS(MainModel,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[ME16_USSL_RMSRatio] = Get_Fractional_RMS(MainModel,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_DNA13_RMSRatio] = Get_Fractional_RMS(MainModel,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[ME16_NWUS_RMSRatio] = Get_Fractional_RMS(MainModel,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[ME16_GAP_P4_RMSRatio] = Get_Fractional_RMS(MainModel,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_GYPSUM_RMSRatio] = Get_Fractional_RMS(MIT_Face1_reformat,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[MITP08_USSL_RMSRatio] = Get_Fractional_RMS(MIT_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_DNA13_RMSRatio] = Get_Fractional_RMS(MIT_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[MITP08_NWUS_RMSRatio] = Get_Fractional_RMS(MIT_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[MITP08_GAP_P4_RMSRatio] = Get_Fractional_RMS(MIT_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_USSL_RMSRatio] = Get_Fractional_RMS(GYPSUM_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_DNA13_RMSRatio] = Get_Fractional_RMS(GYPSUM_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[GYPSUM_NWUS_RMSRatio] = Get_Fractional_RMS(GYPSUM_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[GYPSUM_GAP_P4_RMSRatio] = Get_Fractional_RMS(GYPSUM_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[USSL_DNA13_RMSRatio] = Get_Fractional_RMS(USSL_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[USSL_NWUS_RMSRatio] = Get_Fractional_RMS(USSL_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[USSL_GAP_P4_RMSRatio] = Get_Fractional_RMS(USSL_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[DNA13_NWUS_RMSRatio] = Get_Fractional_RMS(DNA13_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[DNA13_GAP_P4_RMSRatio] = Get_Fractional_RMS(DNA13_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[NWUS_GAP_P4_RMSRatio] = Get_Fractional_RMS(NWUS_Face1_reformat,GAP_P4_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[THBUS_ME16_RMSRatio] = Get_Fractional_RMS(MainModel,THBUS_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_MIT_RMSRatio] = Get_Fractional_RMS(THBUS_Face1_reformat,MIT_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_GYPSUM_RMSRatio] = Get_Fractional_RMS(THBUS_Face1_reformat,GYPSUM_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax)
[THBUS_USSL_RMSRatio] = Get_Fractional_RMS(THBUS_Face1_reformat,USSL_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_DNA13_RMSRatio] = Get_Fractional_RMS(THBUS_Face1_reformat,DNA13_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);
[THBUS_NWUS_RMSRatio] = Get_Fractional_RMS(THBUS_Face1_reformat,NWUS_Face1_reformat,NWUSModelLon,NWUSModelLat,N,Jmax);
[THBUS_GAP_P4_RMSRatio] = Get_Fractional_RMS(THBUS_Face1_reformat,GAP_P4_Face1_reformat,Continental_US_Lon,Continental_US_Lat,N,Jmax);





%%%FIX THE NAMES, THE NAMES ARE WRONG...????

 Frac_RMS_List = [ME16_MIT_RMSRatio; ME16_GYPSUM_RMSRatio; ME16_USSL_RMSRatio; ME16_DNA13_RMSRatio; ME16_NWUS_RMSRatio; ME16_GAP_P4_RMSRatio; MITP08_GYPSUM_RMSRatio;...
    MITP08_USSL_RMSRatio; MITP08_DNA13_RMSRatio; MITP08_NWUS_RMSRatio; MITP08_GAP_P4_RMSRatio; GYPSUM_USSL_RMSRatio; GYPSUM_DNA13_RMSRatio; GYPSUM_NWUS_RMSRatio;...
     GYPSUM_GAP_P4_RMSRatio; USSL_DNA13_RMSRatio; USSL_NWUS_RMSRatio; USSL_GAP_P4_RMSRatio;...
    DNA13_NWUS_RMSRatio; DNA13_GAP_P4_RMSRatio; NWUS_GAP_P4_RMSRatio; THBUS_ME16_RMSRatio; THBUS_MIT_RMSRatio; THBUS_GYPSUM_RMSRatio; THBUS_USSL_RMSRatio; THBUS_DNA13_RMSRatio; THBUS_NWUS_RMSRatio; THBUS_GAP_P4_RMSRatio;];


% %c = categorical(x)
% h = barh(Frac_RMS_List)
% 
% h.FaceColor=[0.7 0.7 0.7]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1] is green
% xlabel('RMS Fraction of Model overlap in North America','FontSize',14)
% grid on
% set(gca,'YTickLabel',[])