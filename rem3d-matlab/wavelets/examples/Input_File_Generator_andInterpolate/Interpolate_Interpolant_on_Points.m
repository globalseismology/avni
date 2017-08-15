% Ved Lekic for Anant Hariharan, July 5th, 2017
% 
% Last edited by Anant Hariharan, 7-27-2017 to fix the step size issue in
% interpolation. 
%
% Based off the second half of Ved Lekic's script. Separated to clarify function, and
% updated to add a little more flexibility.


Interpolant = 'ME16_Vp_Matlab_Interpolant.mat';
Pts_FileName = 'LatLon_CubedSphereN=7.txt';
LoadMe = load(Interpolant);
V = LoadMe.V;
ModelName = 'ME16_Vp'; 

% Now, load in the values at which we want to interpolate
PRI = importdata(Pts_FileName); 

% Now, convert from lon, lat, r to x, y, z
[x,y,z] = sph2cart(pi/180*PRI.data(:,1),pi/180*PRI.data(:,2),PRI.data(:,3)); 
%
vp_vals = zeros(1,length(x)); vs_vals = zeros(size(vp_vals)); 

step_size = 10000; 
%
for j = 1:step_size:length(x)-step_size+1
if(mod(j,step_size)==1), display(['Working on ' num2str(j) 'st point...']); end; 
v_vals(j:j+step_size-1) = V(x(j:j+step_size-1),y(j:j+step_size-1),z(j:j+step_size-1)); 
end 

%duct tape fix. 
for j = length(v_vals):length(x)
    v_vals(j) = V(x(j),y(j),z(j));
end

save([ModelName '_model_at_PRI_locations'],'v_vals'); 

%%%%%%%%%



