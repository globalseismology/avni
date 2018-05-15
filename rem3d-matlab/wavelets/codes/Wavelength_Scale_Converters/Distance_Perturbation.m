%%Get average region of perturbation from different scale wavelets. 

N = 7
Jmax = 4
Scale_Interest= 4;


depth = 0;
re = 6371*1000; %m
NumBox = 2^(2*N) ;  %assuming odd gridding.  %number

PerFaceSurface_Area = (1/6)*((4)*pi*((re-depth*1000)^2)); %m^2

BoxArea = PerFaceSurface_Area/NumBox; %m^2

Length_Dimension = (BoxArea)^(0.5); %m %This corresponds to the difference between two samples.
Length_Dimension  = Length_Dimension/1000 %km

[vwlev,vwlevs] = cube2scale(N,[Jmax+1 Jmax+1],1);
TransformMe = zeros(2^N,2^N,6);

list_lengths = [];
    CurrScaleIndices = find(vwlev == Scale_Interest);
    
        for j = 1:length(CurrScaleIndices)
            TransformMe = zeros(2^N,2^N);
            TransformMe(CurrScaleIndices(j)) = 1;
            recoverMe = angularD4WT(TransformMe,[Jmax Jmax],[1 1],'inverse',1);
            length_indices = find(recoverMe ~= 0);
            length_square = sqrt(length(length_indices));
            list_lengths = [list_lengths length_square];
            100*j/length(CurrScaleIndices)
        end
     Avg_Scal_Length = mean(list_lengths)
     
     
     
  Total_Scal_Length_km = Avg_Scal_Length*Length_Dimension
  Total_Scal_Length_kmRad = Total_Scal_Length_km/2
  GcScalLegthDeg = km2deg(Total_Scal_Length_km)
  GcScalLengthDegRad = GcScalLegthDeg/2