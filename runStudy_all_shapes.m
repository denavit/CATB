clear all; close all; clc;

Fy = 50;
L_over_d = 5:5:50;

% Case 1 - Steel Deck - Infinite Strength
cases(1).beta           = 30;
cases(1).check_strength = false;
% Case 2 - Steel Deck - Typical Strength
cases(2).beta           = 30;
cases(2).check_strength = true;
cases(2).phiMn_deck     = 1;
cases(2).phiTn_conn     = 0.4;
% Case 3 - Composite Deck - Infinite Strength
cases(3).beta           = 400;
cases(3).check_strength = false;
% Case 4 - Composite Deck - Typical Strength
cases(4).beta           = 400;
cases(4).check_strength = true;
cases(4).phiMn_deck     = 5;
cases(4).phiTn_conn     = 10;

code = 'AISC2016';

%% Run calculations
load('ShapeData_Wide_Flange.mat')
numCases   = length(cases);
numShapes  = length(ShapeData_Wide_Flange);
numLengths = length(L_over_d);

tic
fprintf('Running calculations for all shapes...\n');
hwait = waitbar(0,'Initilizing');
for iShape = 1:numShapes
    shapeName = ShapeData_Wide_Flange(iShape).label;
    hwait = waitbar((iShape-1)/numShapes,hwait,...
        sprintf('Running shape %i of %i (%s)',iShape,numShapes,shapeName));
    
    % Create CATB object
    wf = wf_catb(shapeName,Fy,code);
    phiMn_web = 0.9*Fy*wf.tw^2/4;
    
    % Initilize data
    ShapeData_Wide_Flange(iShape).L    = L_over_d*ShapeData_Wide_Flange(iShape).d;
    ShapeData_Wide_Flange(iShape).Pnx  = nan(1,numLengths);
    ShapeData_Wide_Flange(iShape).Pny  = nan(1,numLengths);
    ShapeData_Wide_Flange(iShape).Pnz  = nan(1,numLengths);
    ShapeData_Wide_Flange(iShape).Pnca = nan(1,numLengths);
    ShapeData_Wide_Flange(iShape).beta_Tb_Pnx = nan(1,numLengths);
    ShapeData_Wide_Flange(iShape).beta_Tb_Pnz = nan(1,numLengths);
    ShapeData_Wide_Flange(iShape).Pr_beta = nan(numCases,numLengths);
    
    % Run for each length
    for iL = 1:numLengths
        L = ShapeData_Wide_Flange(iShape).L(iL);
        
        % Nominal Strengths
        ShapeData_Wide_Flange(iShape).Pnx(iL)  = wf.Pnx(L,1.0);
        ShapeData_Wide_Flange(iShape).Pny(iL)  = wf.Pny(L,1.0);
        ShapeData_Wide_Flange(iShape).Pnz(iL)  = wf.Pnz(L,1.0);
        ShapeData_Wide_Flange(iShape).Pnca(iL) = wf.Pnca(L,1.0);
        
        % beta for Pnx
        if ShapeData_Wide_Flange(iShape).Pnx(iL) > ShapeData_Wide_Flange(iShape).Pnca(iL)
            ShapeData_Wide_Flange(iShape).beta_Tb_Pnx(iL) = ...
                wf.beta_Tb(wf.phi_c*ShapeData_Wide_Flange(iShape).Pnx(iL),L);
        else
            ShapeData_Wide_Flange(iShape).beta_Tb_Pnx(iL) = 0;
        end
        
        % beta for Pnz
        if ShapeData_Wide_Flange(iShape).Pnz(iL) > ShapeData_Wide_Flange(iShape).Pnca(iL)
            ShapeData_Wide_Flange(iShape).beta_Tb_Pnz(iL) = ...
                wf.beta_Tb(wf.phi_c*ShapeData_Wide_Flange(iShape).Pnz(iL),L);
        else
            ShapeData_Wide_Flange(iShape).beta_Tb_Pnz(iL) = 0;
        end
        
        % Pr given beta = 0
        ShapeData_Wide_Flange(iShape).Pr_beta0(iL) = ...
            wf.Pr_given_betaTb(L,0);
        
        % Pr given beta cases
        for iCase = 1:numCases           
            beta  = cases(iCase).beta;
            if cases(iCase).check_strength
                phiMn_conn = cases(iCase).phiTn_conn*wf.bf/3*(1/12);
                phiMn = min([cases(iCase).phiMn_deck phiMn_web phiMn_conn]);
                ShapeData_Wide_Flange(iShape).Pr_beta(iCase,iL) = ...
                    wf.Pr_given_betaTb_and_phiMn(L,beta,phiMn);
            else
                ShapeData_Wide_Flange(iShape).Pr_beta(iCase,iL) = ...
                    wf.Pr_given_betaTb(L,beta);
            end
        end
    end
end
fprintf('  calculations complete\n');
close(hwait);

%% Save data
save('Computed_Data','ShapeData_Wide_Flange','-v7.3');
fprintf('  results saved\n');
toc