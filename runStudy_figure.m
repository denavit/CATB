%clear all; close all; clc;

iCase = 4;
compare = 'Pu-to-Pnx';
%compare = 'Pu-to-Pnca';
%compare = 'Pu-to-Pnz';
%compare = 'Pnca-to-Pnx';
%compare = 'Pnca-to-Pnz';

%% Load Data
load('Computed_Data.mat');
numShapes  = length(ShapeData_Wide_Flange);
numLengths = length(ShapeData_Wide_Flange(1).L);
L_over_d   = ShapeData_Wide_Flange(1).L/ShapeData_Wide_Flange(1).d;

%% Retrieve Data
W           = [ShapeData_Wide_Flange(:).W]';
h_over_tw   = [ShapeData_Wide_Flange(:).h_over_tw]';

%% Figure
fs = figureStyle('Thesis');
fs.FontSize = 8;

if strcmp(compare,'Pnca-to-Pnx') || strcmp(compare,'Pnca-to-Pnz')
    hf = fs.figure(3.25,2.25);
    ha = fs.axes([0.14 0.18 0.83 0.79]);
else
    hf = fs.figure(3.5,2.25);
    ha = fs.axes([0.14 0.18 0.83 0.79]);
end

xmin = min(L_over_d);
xmax = max(L_over_d);
%plot([xmin xmax],[1 1],'k-','LineWidth',1);
xlim([xmin xmax]);

set(ha,'XTick',xmin:5:xmax)

c = h_over_tw;
cmax = 60; %max(c);
cmin = 10; %min(c);
caxis([cmin cmax])
num_colors = 100;
colors = gray(num_colors);
colormap(colors)
%colorbar('SouthOutside')

ymin = Inf;
ymin_shape = '';
for i = 1:numShapes
    if W(i) > 150
        continue
    end

    switch compare
        case 'Pu-to-Pnx'
            y = ShapeData_Wide_Flange(i).Pr_beta(iCase,:)./(0.9*ShapeData_Wide_Flange(i).Pnx(:)');
        case 'Pu-to-Pnz'
            y = ShapeData_Wide_Flange(i).Pr_beta(iCase,:)./(0.9*ShapeData_Wide_Flange(i).Pnz(:)');
        case 'Pu-to-Pnca'
            y = ShapeData_Wide_Flange(i).Pr_beta(iCase,:)./(0.9*ShapeData_Wide_Flange(i).Pnca(:)');
        case 'Pnca-to-Pnx'
            y = ShapeData_Wide_Flange(i).Pnca(:)./ShapeData_Wide_Flange(i).Pnx(:);
        case 'Pnca-to-Pnz'
            y = ShapeData_Wide_Flange(i).Pnca(:)./ShapeData_Wide_Flange(i).Pnz(:);
        otherwise
            error('Error')
    end
    
    icolor = ceil(num_colors*(c(i)-cmin)/(cmax-cmin));
    if icolor < 1 || icolor > num_colors
        fprintf('Hello\n')
    end
    icolor = max(icolor,1);
    icolor = min(icolor,num_colors);    
    plot(L_over_d,y,'-','LineWidth',0.5,'Color',colors(icolor,:))
    
    if min(y) < ymin
        ymin = min(y);
        ymin_shape = ShapeData_Wide_Flange(i).label;
    end
end
fprintf('ymin = %g (%s)\n',ymin,ymin_shape)

switch compare
    case 'Pu-to-Pnx'
        ylabel('\it{P_{u,braced}/\phiP_{nx}}')
        ylim([0 3])
    case 'Pu-to-Pnz'
        ylabel('\it{P_{u,braced}/\phiP_{nz}}')
        ylim([0 3])
    case 'Pu-to-Pnca'
        ylabel('\it{P_{u,braced}/\phiP_{nca}}')
        ylim([0 6])
    case 'Pnca-to-Pnx'
        ylabel('\it{\phiP_{nca}/\phiP_{nx}}')
        ylim([0 3])
    case 'Pnca-to-Pnz'
        ylabel('\it{\phiP_{nca}/\phiP_{nz}}')
        ylim([0 1])        
    otherwise
        error('Error')
end
xlabel('Normalized Unbraced Length, \it{L_b/d}')

temp = ylim;
ylim([0 temp(2)])

filename = sprintf('FigureXX_%i_%s',iCase,compare);
print(hf,filename,'-dsvg');


%% Colorbar Figure
% hf = fs.figure(3.5,2.5);
% ha = fs.axes;
% 
% c = h_over_tw;
% cmax = 60; %max(c);
% cmin = 10; %min(c);
% caxis([cmin cmax])
% num_colors = 100;
% colors = gray(num_colors);
% colormap(colors)
% colorbar('SouthOutside')
% set(ha,'Visible','off')
% 
% filename = sprintf('FigureXX_Colorbar');
% print(hf,filename,'-dsvg');