%%Generate arbitary group of repeating lat/lon points on the cubed sphere. 

%99846 lines per depth.
% 99846 = 6 * (129)^2
% 98304 = 6 * (2^(2*7))
%

%FIXED: NOTE: KEEP EO default set to 1. We have an odd number of points!?

% readGNmodel expects 6 * N * N depths. 
% N = 2^lfin + 1: N is the resolution parameter.
% lfin = 7 => N = 129 and it all works out. 

lfin = 7;

[x,y,z]=cube2sphere(lfin);
megalon = [];
megalat = [];

% Each run takes care of one face on the cubed sphere. 

for in=1:6
   
   % Transform these to longitude and latitude
   [phi,piminth,r]=cart2sph(x(:,:,in),y(:,:,in),z(:,:,in));
   lon=phi(:)*180/pi; lat=piminth(:)*180/pi;
   megalon = [megalon; lon]; 
   megalat = [megalat; lat];
   
end

% Write to a file of lats and lons that can be used to create the files,
% given a directory of epix files. See InputFileGenerator.m

fileID = fopen(['Repeating_LatLon_CubedSphere' 'N=' num2str(lfin) '.txt'],'w');
fprintf(fileID,'%s\n', 'lon       lat  ');
for iii = 1:length(megalon)
    fprintf(fileID,'%12f %12f\n', megalon(iii), megalat(iii));
end

