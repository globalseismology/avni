function [x,y,z,megalon,megalat] = Generate_lat_lon_CubedSphere(lfin,eo,write,alfa,bita,gama,sc)

%Generate arbitary group of repeating lat/lon points on the cubed sphere
%using Frederik Simons' subroutines

%Anant Hariharan, Last Updated 13-7-17. 

%----- Explanation for compatibility with GN_Model input model family ---
%99846 lines per depth for lfin = 7. 
% 99846 = 6 * (129)^2
% 98304 = 6 * (2^(2*7)) X: DON'T BE FOOLED BY THIS- UNRELATED!
%

% FIXED: NOTE: KEEP EO default set to 1, when generating wavelet file a la
% PRI5W
% We have an odd number of points!?

% readGNmodel expects 6 * N * N depths. 
% N = 2^lfin + 1: N is the resolution parameter.
% lfin = 7 => N = 129 and it all works out. 
% ---- Explanation ---


defval('lfin',7);
defval('eo',1);
defval('write',0);
defval('alfa',[]);
defval('bita',[]);
defval('gama',[]);
defval('sc',0);

[x,y,z]=cube2sphere(lfin,alfa,bita,gama,eo,sc);
megalon = [];
megalat = [];

% Each run takes care of one chunk on the cubed sphere. 

for in=1:6
   
   % Transform these to longitude and latitude
   [phi,piminth,r]=cart2sph(x(:,:,in),y(:,:,in),z(:,:,in));
   lon=phi(:)*180/pi; lat=piminth(:)*180/pi;
   megalon = [megalon; lon]; 
   megalat = [megalat; lat];
   
end


if write == 1
    % Write to a file of lats and lons that can be used to create the files,
    % given a directory of epix files. See InputFileGenerator.m

    fileID = fopen(['LatLon_CubedSphere' 'N=' num2str(lfin) '.txt'],'w');
    fprintf(fileID,'%s\n', 'lon       lat  ');
    for iii = 1:length(megalon)
        fprintf(fileID,'%12f %12f\n', megalon(iii), megalat(iii));
    end
end
end
