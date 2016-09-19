clear all; close all; clc;

load('Computed_Data.mat')

numShapes  = length(ShapeData_Wide_Flange);
numLengths = length(ShapeData_Wide_Flange(1).L);
numCases   = size(ShapeData_Wide_Flange(1).Pr_beta,1);
shape_labels = {ShapeData_Wide_Flange(:).label}';

% Write to Excel
excel_filename = 'output.xls';

% Main Data
sheet = 'Main';
status = xlswrite(excel_filename,shape_labels,sheet,'A3');

% Pnx
sheet = 'Pnx';
data1 = reshape([ShapeData_Wide_Flange(:).Pnx],numLengths,numShapes)';
data2 = reshape([ShapeData_Wide_Flange(:).beta_Tb_Pnx],numLengths,numShapes)';
data = [data1 data2];
status = xlswrite(excel_filename,shape_labels,sheet,'A3');
status = xlswrite(excel_filename,data,sheet,'B3');

% Pny
sheet = 'Pny';
data = reshape([ShapeData_Wide_Flange(:).Pny],numLengths,numShapes)';
status = xlswrite(excel_filename,shape_labels,sheet,'A3');
status = xlswrite(excel_filename,data,sheet,'B3');

% Pnz
sheet = 'Pnz';
data1 = reshape([ShapeData_Wide_Flange(:).Pnz],numLengths,numShapes)';
data2 = reshape([ShapeData_Wide_Flange(:).beta_Tb_Pnz],numLengths,numShapes)';
data = [data1 data2];
status = xlswrite(excel_filename,shape_labels,sheet,'A3');
status = xlswrite(excel_filename,data,sheet,'B3');

% Pnca
sheet = 'Pnca';
data = reshape([ShapeData_Wide_Flange(:).Pnca],numLengths,numShapes)';
status = xlswrite(excel_filename,shape_labels,sheet,'A3');
status = xlswrite(excel_filename,data,sheet,'B3');

% Case Data
for iCase = 1:numCases
    sheet = sprintf('Case %i',iCase);
    data = nan(numShapes,numLengths);
    for iShape = 1:numShapes
        data(iShape,:) = ShapeData_Wide_Flange(iShape).Pr_beta(iCase,:);
    end
    status = xlswrite(excel_filename,shape_labels,sheet,'A3');
    status = xlswrite(excel_filename,data,sheet,'B3');
end