function createfigure1(cdata1, yvector1, X1, YMatrix1,labelsidetitle)
%CREATEFIGURE1(cdata1, yvector1, X1, YMatrix1)
%  CDATA1:  image cdata
%  YVECTOR1:  bar yvector
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.394957983193277 0.334659090909091 0.530042016806724]);
hold(axes1,'on');

% Create image
image(cdata1,'Parent',axes1,'CDataMapping','scaled');

% Create title
title({'Intermodel Correlation'});

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[0.5 6.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[0.5 28.5]);
box(axes1,'on');
axis(axes1,'ij');
% Set the remaining axes properties
set(axes1,'Layer','top','XTick',[1 2 3 4 5 6],'XTickLabel',...
    {'6','5','4','3','2','1'},'YTick',...
    [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28],...
    'YTickLabel',...
    {'MITP08-ME16','GYPSUM-ME16','USSL-ME16','DNA13-ME16','NWUS-ME16','GAP\_P4-ME16','GYPSUM-MITP08','USSL-MITP08','DNA13-MITP08','NWUS-MITP08','GAP\_P4-MITP08','USSL-GYPSUM','DNA13-GYPSUM','NWUS-GYPSUM','GAP\_P4-GYPSUM','DNA13-USSL','NWUS-USSL','GAP\_P4-USSL','NWUS-DNA13','GAP\_P4-DNA13','GAP\_P4-NWUS', 'ME16-PTHB\_NA16','MITP08-PTHB\_NA16','GYPSUM-PTHB\_NA16', 'USSL-PTHB\_NA16','DNA13-PTHB\_NA16','NWUS-PTHB\_NA16','GAP\_P4-PTHB\_NA16'});
% Create colorbar
hhhh = colorbar('peer',axes1,'Position',...
    [0.576642335766424 0.395261845386534 0.024330900243309 0.529925187032419]);


ylabel(hhhh, labelsidetitle)
colormap(flipud(inferno))

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.47625 0.394014962593516 0.0857937956204375 0.531172069825436]);
hold(axes2,'on');

% Create bar
barh(flipud(yvector1),'BarWidth',0.9,'Parent',axes2,'FaceColor',[0.7 0.7 0.7]);

%barh(ones(length(yvector1),1),'Parent',axes2,'FaceColor',[1 1 1])
plot([1 1],[0 100],'--','Color','k','Parent',axes2,'LineWidth',0.75)
% Create title
title({'RMS Ratios '});

% Uncomment the following line to preserve the X-limits of the axes
xlim(axes2,[0 7.5]);
% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes2,[0.5 28.5]);
box(axes2,'on');
grid(axes2,'on');
% Set the remaining axes properties
set(axes2,'YTickLabel','');
% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.13 0.0952380952380953 0.334659090909091 0.233333333333333]);
hold(axes3,'on');

% Create multiple lines using matrix input to semilogy
semilogy1 = semilogy(X1,YMatrix1,'Parent',axes3,'LineWidth',2);
set(semilogy1(1),'DisplayName','MITP08');
set(semilogy1(2),'DisplayName','ME16');
set(semilogy1(3),'DisplayName','USSL-2014');
set(semilogy1(4),'DisplayName','DNA13');
set(semilogy1(5),'DisplayName','NWUS');
set(semilogy1(6),'DisplayName','GYPSUM');
set(semilogy1(7),'DisplayName','GAP\_P4');
set(semilogy1(8),'DisplayName','PTHB\_NA16','Color',[1 0.7 0.9]);

% Create ylabel
ylabel('RMS Amplitude','FontSize',12);

% Create xlabel
xlabel('Wavelet Scale','FontSize',12);

% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes3,[0.01 1]);
box(axes3,'on');
grid(axes3,'on');
% Set the remaining axes properties
set(axes3,'XDir','reverse','XTick',[1 2 3 4 5 6],'XTickLabel',...
    {'1','2','3','4','5','6'},'YMinorTick','on','YScale','log');
% Create legend
legend1 = legend(axes3,'show');
set(legend1,...
    'Position',[0.47625 0.0942307692307692 0.1 0.234615384615385]);

