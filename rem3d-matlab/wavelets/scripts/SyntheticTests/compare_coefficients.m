%%COMPARE WAVELET COEFFICIENTS FOR DIFFERENT INPUT MODELS. 
%%Currently set to compare coefficients for Lines. 

Get_Edgy = 1;
N = 9; Jmax = 5; face = 1;
[vwlev,vwlevs] = cube2scale(N,[Jmax Jmax],1);
analyse_truncated_line = 1;
%-------

Desired_Coefficient = 1000;

%%LINES - Get Indices Here%%
Desired_Indexbigline_dim1 = [64:84];
Desired_Indexbigline_dim2 = [64*ones(21,1)];


TransformMe_BigLine = zeros(2^N,2^N,6);
for i = 1:length(Desired_Indexbigline_dim1)
    %disp(i)
    TransformMe_BigLine(Desired_Indexbigline_dim1(i),Desired_Indexbigline_dim2(i),face) = Desired_Coefficient;
end

%%LINES - Get Indices Here%%
Desired_Indexsmallline_dim1 = [64:74];
Desired_Indexsmallline_dim2 = [64*ones(11,1)];
%%%%%

TransformMe_smallLine = zeros(2^N,2^N,6);
for i = 1:length(Desired_Indexsmallline_dim1)
    %disp(i)
    TransformMe_smallLine(Desired_Indexsmallline_dim1(i),Desired_Indexsmallline_dim2(i),face) = Desired_Coefficient;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vws_Big = angularD4WT(TransformMe_BigLine,[Jmax Jmax],[1 1],'forward',1);

if analyse_truncated_line ==1
[ xpts,ypts ] = Get_Interior_PolygonPts_on_CubedSphere( N,[64 64 74 74],[64 64 64 64] );
%Detect Wavelets where we want them to be. 
if Get_Edgy == 1
    [ wavelet_locs ] = get_wavelet_pos( xpts,ypts,N,Jmax,face,0,1);
else
    [ wavelet_locs ] = get_wavelet_pos( xpts,ypts,N,Jmax,face,0,0);

end

%Now zero out wavelets at those locations. 
[vws_Big ] = zero_wavelets(vws_Big,wavelet_locs,'db4',N,1,Jmax );
end

vws_Small = angularD4WT(TransformMe_smallLine,[Jmax Jmax],[1 1],'forward',1);


Big_Indices = find(vws_Big ~= 0);
Big_Wavelets = vws_Big(Big_Indices);

Small_Indices = find(vws_Small ~= 0);
Small_Wavelets = vws_Small(Small_Indices);

Scal1Indices = find(vwlevs == 1);
Scal2Indices = find(vwlevs == 2);
Scal3Indices = find(vwlevs == 3);
Scal4Indices = find(vwlevs == 4);




subplot(2,2,[1,2])

for i = 1:length(Big_Indices)
    

if ismember(Big_Indices(i),Scal1Indices)
plot(Big_Indices(i),Big_Wavelets(i),'r+')
hold on
elseif ismember(Big_Indices(i),Scal2Indices)
plot(Big_Indices(i),Big_Wavelets(i),'b+')
hold on
elseif ismember(Big_Indices(i),Scal3Indices)
plot(Big_Indices(i),Big_Wavelets(i),'g+')
hold on
elseif ismember(Big_Indices(i),Scal4Indices)
plot(Big_Indices(i),Big_Wavelets(i),'c+')
hold on
end

end



hold on
for i = 1:length(Small_Indices)
if ismember(Small_Indices(i),Scal1Indices)
plot(Small_Indices(i),Small_Wavelets(i),'ro')
elseif ismember(Small_Indices(i),Scal2Indices)
plot(Small_Indices(i),Small_Wavelets(i),'bo')
elseif ismember(Small_Indices(i),Scal3Indices)
plot(Small_Indices(i),Small_Wavelets(i),'go')
elseif ismember(Small_Indices(i),Scal4Indices)
plot(Small_Indices(i),Small_Wavelets(i),'co')
end

end



title('Wavelet Coefficient Points for Half Lines: Predicted from Synthetic and Inverted from a direct Half-Line\n + is zeroed out')
xlabel('Index on cubed Sphere')
ylabel('Coefficient Value')


%%%%% CHECK %%%%%

Recover_Direct = angularD4WT(vws_Small,[Jmax Jmax],[1 1],'inverse',1);
Recover_Circuitours = angularD4WT(vws_Big,[Jmax Jmax],[1 1],'inverse',1);

subplot(2,2,3)
l=imagefnan([1 1],[2^N 2^N],Recover_Direct(:,:,face));
title('Direct Recovery from Inverse')
subplot(2,2,4)
j=imagefnan([1 1],[2^N 2^N],Recover_Circuitours(:,:,face));
title('Recovery from Truncated Large Line')