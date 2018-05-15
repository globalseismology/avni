function [ Model ] = ReadEpixArr(EpixFile)
%Anant Hariharan, ah824-at-cornell.edu

%Ignore top two rows- should just be comments!

%=================== input model dv =======================================
%
fid = fopen(EpixFile, 'r');
A = [textscan(fid, '%f %f %f %f','HeaderLines',2),1]; fclose(fid);
% lat
Model.lat    = A{1};

% londata
Model.lon   = A{2};
% pixsize
Model.pixsize    = A{3};

%dv anomaly!
Model.dv= A{4};

%Number of dv values in parametrization
Model.nblobs = length(Model.dv);

