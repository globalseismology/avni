%Outputs .mat file for grid. This contains all the points and their
%lat,lon,and index. Along with Metadata. 
clear all

eo = 0;
N = [7 8 9];
Jmax = [4 5 6];
config = 1; %Which euler angle config do you want to use?
SuperChunkMe = 1 %Output a superchunk grid?

if config ==1
alfa = 0.2089-0.175;
bita= 0.9205+0.25;
gama= 1.2409-0.05;
end


for iii = 1:length(N)

lfin = N(iii);


if SuperChunkMe == 1
[x,y,z]=cube2sphere(lfin,alfa,bita,gama,eo,SuperChunkMe);    
else    
[x,y,z]=cube2sphere(lfin,alfa,bita,gama,eo);
end
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

[vwlev,vwlevs]=cube2scale(N(iii),[Jmax(iii)+1 Jmax(iii)+1],1);

Face_Interval = length(megalon)/6;

Out_Struc.lat = megalat;
Out_Struc.lon = megalon;
for i=1:length(megalon)
    if i <= Face_Interval
        Out_Struc.face(i) = 1;
    elseif i <= Face_Interval*2
        Out_Struc.face(i) = 2;
    elseif i <= Face_Interval*3
        Out_Struc.face(i) = 3;
    elseif i <= Face_Interval*4
        Out_Struc.face(i) = 4;
    elseif i <= Face_Interval*5
        Out_Struc.face(i) = 5;
    elseif i <= Face_Interval*6
        Out_Struc.face(i) = 6;
    end
    
    Out_Struc.ScaleIndex(i) = vwlevs(i);
    
end
Out_Struc.face = Out_Struc.face';
Out_Struc.MetaN = N(iii);
Out_Struc.MetaJmax = Jmax(iii);
Out_Struc.Metaeo = eo;
Out_Struc.MetaEulerVals(1) = alfa;
Out_Struc.MetaEulerVals(2) = bita;
Out_Struc.MetaEulerVals(3) = gama;
Out_Struc.MetaEulerNames{1} = 'alfa';
Out_Struc.MetaEulerNames{2} = 'bita';
Out_Struc.MetaEulerNames{3} = 'gama';

if SuperChunkMe == 1
FileName = ['SC_Grid_' 'N' num2str(N(iii)) '_Jmax' num2str(Jmax(iii)) '_EulerConfig' num2str(config)];
else
FileName = ['Grid_' 'N' num2str(N(iii)) '_Jmax' num2str(Jmax(iii)) '_EulerConfig' num2str(config)];
end
save(FileName,'-struct','Out_Struc');


end