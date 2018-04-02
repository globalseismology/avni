function [] = plot_ncmodels(filename)
%Orig_Model_Plotter
% Plots nc tomographic models. 

figure()

lon = ncread(filename,'longitude');
lat = ncread(filename,'latitude');
depth = ncread(filename,'depth');
depth_interest = 100;
depthind = find(depth==depth_interest);
dvp = ncread(filename,'dvp');
dvp(isnan(dvp))=0;
%%


lonmat = repmat(lon,[1,length(lat)]);
latmat = repmat(lat',[length(lon),1]);
maxlat =  max(lat);
minlat = min(lat);
minlon = min(lon);
maxlon =  max(lon);




contourf(lonmat,latmat,dvp(:,:,depthind),20)
colorbar

grid on
box on
xlabel('Longitude')
ylabel('Latitude')
title(['Depth =' num2str(depth_interest) ' km' ] )



%%%%%%%
biglon = [-180:0.2:180]';
biglat = [-90:0.2:90]';

biglonmat = repmat(biglon,[1,length(biglat)]);
biglatmat = repmat(biglat',[length(biglon),1]);


end

