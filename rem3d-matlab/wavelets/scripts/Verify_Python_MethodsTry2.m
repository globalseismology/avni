
ModelFile = load('ME16_Vp.N7.Jmax4.EulerConfig1.mat');
WaveletBasis = load('/home/anant/Software/rem3d/rem3d/files/DVmaps.N7.J4.D4.mat');
GridFile = load('/home/anant/Software/rem3d/rem3d-matlab/wavelets/examples/Input_File_Generator_andInterpolate/DotMatFileGenerators/Grid_Database/Grid_N7_Jmax4_EulerConfig1.mat');



%%%Verify Python Methods Try 2: 
Jmax = 4; N = 7;
MIT_Results_90_Km = load('/home/anant/Software/rem3d/rem3d-matlab/wavelets/examples/MIT_Wave_Results/loris5_GN_0090_7_4_D4_1_1.mat')
Face1Coeffs = MIT_Results_90_Km.vw(:,:,1);
[vwlev,vwlevs]=cube2scale(N,[5 5],1)

%%%%%%%%%%%%%%
All_Scale4Indices = find(vwlev == 4);
%%%%%%%%%%%%%%

Out_Map = 0;
for i = 1:length(All_Scale4Indices)
    WvCoeff = Face1Coeffs(All_Scale4Indices(i));
    Curr_Map = WvCoeff*full(WaveletBasis.Me((All_Scale4Indices(i))).dv);
    Out_Map = Out_Map + Curr_Map;
end


plotoncube(Out_Map,'2D')

figure

plotoncube(angularD4WT(Face1Coeffs,[Jmax Jmax],[1 1],'inverse',1),'2D')
