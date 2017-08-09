%Test out the capabilities of superchunks!!!


lfin = 7;
eo = 0;
write = 0
alfa = [];
bita = []
gama = []
sc = 0;


[x,y,z,megalonCSPH,megalatCSPH] = Generate_lat_lon_CubedSphere(lfin,eo,write,alfa,bita,gama,sc);
[x,y,z,megalonSC,megalatSC] = Generate_lat_lon_CubedSphere(lfin,eo,write,alfa,bita,gama,1);

figure

subplot(1,2,1)
scatter(megalonCSPH(1:length(megalonCSPH)/6),megalatCSPH(1:length(megalonCSPH)/6),5,'r')
hold on
scatter(megalonCSPH(1+length(megalonCSPH)/6:2*length(megalonCSPH)/6),megalatCSPH(1+length(megalatCSPH)/6:2*length(megalatCSPH)/6),5,'b')
scatter(megalonCSPH(1+2*length(megalonCSPH)/6:3*length(megalonCSPH)/6),megalatCSPH(1+2*length(megalatCSPH)/6:3*length(megalatCSPH)/6),5,'g')
scatter(megalonCSPH(1+3*length(megalonCSPH)/6:4*length(megalonCSPH)/6),megalatCSPH(1+3*length(megalatCSPH)/6:4*length(megalatCSPH)/6),5,'c')
scatter(megalonCSPH(1+4*length(megalonCSPH)/6:5*length(megalonCSPH)/6),megalatCSPH(1+4*length(megalatCSPH)/6:5*length(megalatCSPH)/6),5,'k')
scatter(megalonCSPH(1+5*length(megalonCSPH)/6:6*length(megalonCSPH)/6),megalatCSPH(1+5*length(megalatCSPH)/6:6*length(megalatCSPH)/6),5,'f')

hold on
%scatter(megalonCSPH,megalatCSPH,'blue')
title('Original Construction: Cubed Sphere')
%xlim([-180,180])
%ylim([-90,-50])
subplot(1,2,2)

scatter(megalonSC(1:length(megalonCSPH)/6),megalatSC(1:length(megalonCSPH)/6),5,'r')
hold on
scatter(megalonSC(1+length(megalonCSPH)/6:2*length(megalonCSPH)/6),megalatSC(1+length(megalatCSPH)/6:2*length(megalatCSPH)/6),5,'b')
scatter(megalonSC(1+2*length(megalonCSPH)/6:3*length(megalonCSPH)/6),megalatSC(1+2*length(megalatCSPH)/6:3*length(megalatCSPH)/6),5,'g')
scatter(megalonSC(1+3*length(megalonCSPH)/6:4*length(megalonCSPH)/6),megalatSC(1+3*length(megalatCSPH)/6:4*length(megalatCSPH)/6),5,'c')
scatter(megalonSC(1+4*length(megalonCSPH)/6:5*length(megalonCSPH)/6),megalatSC(1+4*length(megalatCSPH)/6:5*length(megalatCSPH)/6),5,'k')
scatter(megalonSC(1+5*length(megalonCSPH)/6:6*length(megalonCSPH)/6),megalatSC(1+5*length(megalatCSPH)/6:6*length(megalatCSPH)/6),5,'f')

title('Second Construction: Superchunked Sphere')
%xlim([-180,180])
%ylim([-90,-50])