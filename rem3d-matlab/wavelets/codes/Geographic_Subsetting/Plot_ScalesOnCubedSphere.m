function [ right_indices ] = Plot_ScalesOnCubedSphere( N,Jmax,J_Interest )
% Makes a map view plot of the location of scale N points on the cubed
% sphere. 

%Inputs:

%J_Interest is a vector of all the scales you're interested in plotting. 
% N is the resolution parameter in all its ambiguity. 
% J is the max you want to go in the expansion. 

[vwlev, vwlevs] = cube2scale(N,[Jmax Jmax],1);
figure
for i = 1:length(J_Interest)
    right_indices = (vwlevs == J_Interest(i))   
    make2d_map_plot(right_indices,N,Jmax)
    hold on
end
title(['scales of interest:' num2str(J_Interest) ' N = ' num2str(N) ' Full J =' num2str(Jmax)])
