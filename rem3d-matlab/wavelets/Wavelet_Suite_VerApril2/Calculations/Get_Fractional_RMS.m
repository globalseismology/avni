function [RMS_Ratio] = Get_Fractional_RMS(Model1,Model2,Geographic_Lon,Geographic_Lat,N,Jmax)
% Gets RMSs for each model and divides.
RMS_1 = D4wt_Localized_Model_RMS(Model1,Geographic_Lon,Geographic_Lat,N,Jmax);
RMS_2 = D4wt_Localized_Model_RMS(Model2,Geographic_Lon,Geographic_Lat,N,Jmax);
RMS_Ratio = RMS_1/RMS_2;
end

