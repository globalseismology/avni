%% Read Scott's model!

THBUS = load('BurdickLekic2017.mat');

ime = 'THBUS_Matlab_Interpolant'
dvp = vmeanall;
%dvp(isnan(dvp))=0;
lon = lons;
lat = lats;
depth = 6371-radii;



counter = 0
for k = 1:length(depth)

            counter = counter + 1;
            mdl(counter,1) = lon(counter); 
            mdl(counter,2) = lat(counter); 
            mdl(counter,3) = depth(counter); 
            mdl(counter,4) = dvp(counter);

end
mdl=double(mdl);
% Now, conver from spherical coordinates to X,Y,Z (cartesian)
[x,y,z] = sph2cart(pi/180*mdl(:,1),pi/180*mdl(:,2),6371-mdl(:,3)); 
    
V = scatteredInterpolant(x(:),y(:),z(:),mdl(:,4),'natural');
clear mdl; 
save(ime,'V'); 
