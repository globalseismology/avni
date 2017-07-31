function [ xpts,ypts ] = Get_Interior_PolygonPts_on_CubedSphere( N,xbounds,ybounds,plot )
%Get x and y vectors that correspond to points within a polygon, where the
%polygon is defined on the cubed sphere.
defval('N',7)
defval('plot',0)
% Anant Hariharan, 7-24-17

biglinex = [];
for i = 1:2^N
    line = 1:2^N;
    biglinex = [biglinex line];
end

bigliney = [];
for i = 1:2^N
    line = i*ones(1,2^N);
    bigliney = [bigliney line];
end

in = inpolygon(biglinex,bigliney,xbounds,ybounds);
%figure

xpts = biglinex(in);
ypts = bigliney(in);

if plot ==1
plot(biglinex(in),bigliney(in),'r+') % points inside
hold on
plot(biglinex(~in),bigliney(~in),'bo') % points outside
end
end

