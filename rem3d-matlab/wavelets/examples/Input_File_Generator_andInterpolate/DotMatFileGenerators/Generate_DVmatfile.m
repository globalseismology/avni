%Generate DV .mat file. 
clear all
%GridFile
%%ALL YOU NEED!
N = 7;
Jmax = 4;
eo = 0;
config = 1;
SuperChunkMe = 0; %Output a superchunk grid?
%%
if config ==1
alfa = 0.2089-0.175;
bita= 0.9205+0.25;
gama= 1.2409-0.05;
end
%%
[vwlev,vwlevs]=cube2scale(N,[Jmax+1 Jmax+1],1);

if SuperChunkMe == 1
Name = ['/home/anant/mydbs/Grid_Database/SC_Grid_N' num2str(N) '_Jmax' num2str(Jmax) '_EulerConfig1.mat'];
else
Name = ['/home/anant/mydbs/Grid_Database/Grid_N' num2str(N) '_Jmax' num2str(Jmax) '_EulerConfig1.mat'];    
end

ModelName = 'ME16_Vp'; 
Interpolant = [ModelName '_Matlab_Interpolant.mat']
%Interpolant = 'ME16_Vp_Matlab_Interpolant.mat';
Grid = load(Name);
LoadMe = load(Interpolant);
V = LoadMe.V;


% Depths = [  25
%    100
%    200
%    300
%    400
%    500
%    600
%    700
%    800
%    900
%   1000
%   1100
%   1200
%   1300
%   1400
%   1500
%   1600
%   1700
%   1800
%   1900
%   2000
%   2100
%   2200
%   2300
%   2400
%   2500
%   2600
%   2700
%   2800
%   2890];


Depths = [  25
   100
   200
   300
   400
   500
   600
   700
   800
   900
   1000
   1100
   1200
   1300
   1400
   1500
   1600
   1700
   1800
   1900
   2000
   2100
   2200
   2300
   2400
   2500
   2600
   2700
   2800];


r = 6371 - Depths;
Wavelets = [];
z = [];
lon = [];
lat = [];
ScaleArr = []
for i = 1:length(r)
 TempScaleArr = vwlevs(:);
 latlength = 6*2^(2*N);
 one_array = ones(1,latlength);
 Depth_arr = r(i) * one_array;
 z = [z Depth_arr];
 lon = [lon; Grid.lon];
 lat = [lat; Grid.lat];
 ScaleArr = [ScaleArr; TempScaleArr];
end
 z = z';
 z_arr = z;


%%%% Make it better
%%%%

    % Now, convert from lon, lat, r to x, y, z
    [x,y,z] = sph2cart(pi/180*lon,pi/180*lat,z); 
    %
    vp_vals = zeros(1,length(x)); vs_vals = zeros(size(vp_vals)); 

    step_size = 10000; 
    %
    for j = 1:step_size:length(x)-step_size+1
    100*j/(length(x)-step_size+1)
    if(mod(j,step_size)==1), display(['Working on ' num2str(j) 'st point...']); end; 
    v_vals(j:j+step_size-1) = V(x(j:j+step_size-1),y(j:j+step_size-1),z(j:j+step_size-1)); 
    end 

    %duct tape fix. 
    for j = length(v_vals):length(x)
        100*j/length(x)
        v_vals(j) = V(x(j),y(j),z(j));
    end
    
    
  %%%%
  

Out_Struc.MetaN = N;
Out_Struc.MetaJmax = Jmax;
Out_Struc.Metaeo = eo;
Out_Struc.MetaEulerVals(1) = alfa;
Out_Struc.MetaEulerVals(2) = bita;
Out_Struc.MetaEulerVals(3) = gama;
Out_Struc.MetaEulerNames{1} = 'alfa';
Out_Struc.MetaEulerNames{2} = 'bita';
Out_Struc.MetaEulerNames{3} = 'gama';
Out_Struc.model = v_vals;
Out_Struc.depth = z_arr;
Out_Struc.ScaleIndex = ScaleArr;
Out_Struc.lat = lat;
Out_Struc.lon = lon;
%%%%%

for i = 1:length(r)
    Incr_v = Out_Struc.model(1+(i-1)*6*2^(2*N):(i)*6*2^(2*N));
    Csph_v = reshape(Incr_v,[2^N 2^N 6]);
    Wv_Coeffs = angularD4WT(Csph_v,[Jmax Jmax],[1 1],'forward',1);
    Wavelets = [Wavelets; Wv_Coeffs(:)];
end
Wavelets = Wavelets';
Out_Struc.wvcoeffs = Wavelets;

if SuperChunkMe == 1
FileName = ['/home/anant/Software/rem3d/rem3d/files/' 'SC_' ModelName '.N' num2str(N) '.Jmax' num2str(Jmax) '.EulerConfig' num2str(config) '.mat'];
else
FileName = ['/home/anant/Software/rem3d/rem3d/files/' 'Csph_' ModelName '.N' num2str(N) '.Jmax' num2str(Jmax) '.EulerConfig' num2str(config) '.mat'];
end

save(FileName,'-struct','Out_Struc');

    