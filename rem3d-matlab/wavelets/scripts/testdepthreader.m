%%%%%TESTS!
d=load(fullfile(getenv('IFILES'),'EARTHMODELS','MONTELLI','GNdepths'));

indix = 20

depth_interest = d(indix)/1000
depth_interest = d(indix)/100 + 1
awwyeah = readGNmodel(depth_interest)