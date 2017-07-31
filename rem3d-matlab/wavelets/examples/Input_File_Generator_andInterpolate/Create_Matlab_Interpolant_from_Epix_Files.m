% Create a MATLAB interpolant object from a set of epix files. 
% Ved Lekic for Anant Hariharan, July 5th, 2017
% 
% Last edited by Anant Hariharan, 7-27-2017

clear; close all; 

% Directory with epix files (first column in Lat, second is Lon)
epix_dir = 'ME16_Vp_epix'; 

% Name of interpolant to save
ime = 'ME16_Vp_Matlab_Interpolant'; 

% Work happens below here!
lista = dir([epix_dir '/*.epix']); 

% Read them in, one at a time
for j = 1:length(lista)
    a = importdata([epix_dir '/' lista(j).name]); 
    % Find the depth as number between the periods
    indx1 = strfind(lista(j).name,'.'); 
    dep = str2num(lista(j).name(indx1(1)+1:indx1(2)-1));
    n = size(a.data,1); 
    indx = (j-1)*n + 1:j*n; 
    % In my structure, first column is lon, second is lat, third is dep,
    % and fourth is the velocity perturbation.
    mdl(indx,1) = a.data(:,2); 
    mdl(indx,2) = a.data(:,1); 
    mdl(indx,3) = dep; 
    mdl(indx,4) = a.data(:,4); 
    
    % display progress
    100*j/length(lista)
end

% Now, conver from spherical coordinates to X,Y,Z (cartesian)
[x,y,z] = sph2cart(pi/180*mdl(:,1),pi/180*mdl(:,2),6371-mdl(:,3)); 
    
V = scatteredInterpolant(x(:),y(:),z(:),mdl(:,4),'natural');
clear mdl; 
save(ime,'V'); 



