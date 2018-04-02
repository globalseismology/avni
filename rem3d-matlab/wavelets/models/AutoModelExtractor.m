function [V] = AutoModelExtractor(filename,ime)

%Automates reading/extration of nc tomo models, pre-interpolation

% clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of interpolant to save
%ime = 'NWUS11-P_Matlab_Interpolant'; 
%filename = 'NWUS11-P_percent.nc';
dvp = ncread(filename,'dvp');
dvp(isnan(dvp))=0;
lon = ncread(filename,'longitude');
lat = ncread(filename,'latitude');
depth = ncread(filename,'depth');

% maxlat =  50.0;
% minlat = 24.0;
% minlon = -125.0;
% maxlon =  -67.0;

counter = 0
for k = 1:length(depth)
    for i = 1:length(lon)
        for j = 1:length(lat)
            counter = counter + 1;
            mdl(counter,1) = lon(i); 
            mdl(counter,2) = lat(j); 
            mdl(counter,3) = depth(k); 
            mdl(counter,4) = dvp(i,j,k);
        end
    end
    k
end
mdl=double(mdl);
% Now, conver from spherical coordinates to X,Y,Z (cartesian)
[x,y,z] = sph2cart(pi/180*mdl(:,1),pi/180*mdl(:,2),6371-mdl(:,3)); 
    
V = scatteredInterpolant(x(:),y(:),z(:),mdl(:,4),'natural');
clear mdl; 
save(ime,'V'); 


end

