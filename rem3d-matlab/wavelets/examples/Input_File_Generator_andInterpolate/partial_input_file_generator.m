%input file generator

dvs = load('MIT_P08_epix_model_at_PRI_locations.mat')

header1 = '   '
fid=fopen('mitdvsonly.txt','w');
fprintf(fid, [ header1 '\n']);
fprintf(fid, '%f \n', [v]');
fclose(fid);