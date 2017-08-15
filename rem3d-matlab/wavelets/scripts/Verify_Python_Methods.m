Primal_Test = 1;
N = 7;
J = 4;
[vwlev,vwlevs]=cube2scale(N,[5 5],1)
All_Scale4Indices = find(vwlev == 4);
%%%%%%%%%%%%%%%%Reconstruct Face Entirely



if Primal_Test ==1
N = 7
Indices = 2^(2*N)
Out_Map = 0;
for i = 1:length(All_Scale4Indices)
    %WvCoeff = ModelFile.wvcoeffs(i);
    WvCoeff = ModelFile.wvcoeffs(All_Scale4Indices(i));
    Curr_Map = WvCoeff*full(WaveletBasis.Me(All_Scale4Indices(i)).dv);
    Out_Map = Out_Map + Curr_Map;
end

%contourf(Out_Map)
plotoncube(Out_Map,'2D')
else

ModelFile = load('ME16_Vp.N7.Jmax4.EulerConfig1.mat');
WaveletBasis = load('/home/anant/Software/rem3d/rem3d/files/DVmaps.N7.J4.D4.mat');
GridFile = load('/home/anant/Software/rem3d/rem3d-matlab/wavelets/examples/Input_File_Generator_andInterpolate/DotMatFileGenerators/Grid_Database/Grid_N7_Jmax4_EulerConfig1.mat');


%%%
Dep_Interest = 100;
Radius_Interest = 6371-Dep_Interest;
%%%
Dep_Indices = find(ModelFile.depth==Radius_Interest);
Scale_Interest= 4;
%%%%
Indices_Interest = find(GridFile.ScaleIndex == Scale_Interest)

Face1Map = 0;
for i = 1:length(Indices_Interest)
    Curr_Index = Indices_Interest(i);
    Curr_Face = GridFile.face(Curr_Index);
    if Curr_Face ==1
        CurrMap = ModelFile.wvcoeffs(Curr_Index)* WaveletBasis.Me(Curr_Index).dv;
        Face1Map = Face1Map + CurrMap;
    end
end


contourf(Face1Map)
colorbar

end



