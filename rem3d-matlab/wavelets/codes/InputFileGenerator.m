%Input File Generation
%
% Anant Hariharan, 7-8-17


%%%%%%%%%%%Get the repeating list of lats and lons 

VeloFile = 'PRI5w';

fid = fopen(VeloFile, 'r');
A = [textscan(fid, '%f %f %f %f','HeaderLines',1),1]; fclose(fid);
% lat
Model.lon    = A{1};

% londata
Model.lat   = A{2};
% radius
Model.rad    = A{3};

%dv anomaly!
Model.dv= A{4};

%Number of dv values in parametrization
Model.nblobs = length(Model.dv);

repeating_lons = Model.lon(1:99846);
repeating_lats = Model.lat(1:99846);
%%%%%%%%%%%%



% Directory with epix files (first column in Lat, second is Lon)
epix_dir = 'ME16_Vp_epix'; 


Real_DV_Anomalies = load('ME16_Vp_epix_model_at_PRI_locations.mat');
velocity_anomalies = Real_DV_Anomalies.v_vals;

% Directory with epix files (first column in Lat, second is Lon)
epix_dir = 'ME16_Vp_epix'; 

% Work happens below here!
depths = load('GN_Depths.txt')


fileID = fopen('V1ME16_dvvals.txt','w');
fprintf(fileID,'%s\n', 'lon       lat         r     dv(%)');
counter = 0
for j = 1:length(depths)
    dep = depths(j)/1000
    rad = 6371 - dep
    for jj = 1:length(repeating_lons)
        counter = counter + 1;
        %%%%WRITE TO FULL FILE
        if counter <= length(velocity_anomalies)
        fprintf(fileID,'%12f %12f %12f %12f\n', repeating_lons(jj), repeating_lats(jj), rad, velocity_anomalies(counter));
        else
        fprintf(fileID,'%12f %12f %12f %12f\n', repeating_lons(jj), repeating_lats(jj), rad, velocity_anomalies(length(velocity_anomalies)));    
        end            
    end
end
