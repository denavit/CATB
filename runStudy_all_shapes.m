clear all; close all; clc;

Fy = 50;
L_over_d = 10:10:100;
code = 'AISC2016';

%% Run Analysis
load('ShapeData_Wide_Flange.mat')

% Open file and write header line
fid = fopen('beta_for_all_shapes_results.csv','w');
fprintf(fid,'shape,h/t');
for i = 1:length(L_over_d)
    fprintf(fid,',L/d = %i',L_over_d(i));
end
fprintf(fid,'\n');

% Compute for each shape
for iShape = 1:length(ShapeData_Wide_Flange)
    
    shapeName = ShapeData_Wide_Flange(iShape).label;
    wf = wf_caftb(shapeName,Fy,code);
    
    fprintf(fid,'%s,%g',shapeName,wf.h_over_tw);
    for iL = 1:length(L_over_d)
        L       = L_over_d(iL)*wf.d;
        phiPnx  = wf.phi_c*wf.Pnx(L,1.0);
        phiPnca = wf.phi_c*wf.Pnca(L,1.0);
        if phiPnx > phiPnca
            beta_Tb = wf.beta_Tb(phiPnx,L);
        else
            beta_Tb = 0;
        end
        fprintf(fid,',%g',beta_Tb);
    end
    fprintf(fid,'\n');
    
end

% Close file
fclose(fid);

