%%%% Write to an output file of the lon,lat,r values.
%Given N. Given Euler Angles. 



lfin = 7; %basically, N, but lets follow Simons' convention here. 
alfa = [];
bita = [];
gama = [];
eo =1;
sc =0;
depths = [0 100 200]
r = 6371 - depths;

[x,y,z,megalon,megalat] = Generate_lat_lon_CubedSphere(lfin,eo,0);

% Write to a file of lats and lons that can be used to create the files,
% given a directory of epix files. See InputFileGenerator.m

fileID = fopen(['LatLon_CubedSphere' 'N=' num2str(lfin) '.txt'],'w');
fprintf(fileID,'%s\n', 'lon       lat  ');
for jjj = 1:length(r)
for iii = 1:length(megalon)
    fprintf(fileID,'%12f %12f %12f %12f\n', megalon(iii), megalat(iii), r(jjj));
end
end