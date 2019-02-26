clear all; close all; clc

shape_name = 'W18x35';
Fy = 50;
K  = 1;
wf = wf_caftb(shape_name,Fy,'AISC2016');
Py = wf.Py;
bf = wf.bf;

phiMn_web = 0.9*Fy*wf.tw^2/4;
phiMn_deck_A = 5;
phiMn_deck_B = 1;
phiPn_conn_A = 10;
phiPn_conn_B = 0.4;

phiMn_A = min([phiMn_deck_A phiPn_conn_A*(bf/3)*(1/12) phiMn_web]);
phiMn_B = min([phiMn_deck_B phiPn_conn_B*(bf/3)*(1/12) phiMn_web]);

beta_A = 400;
beta_B = 30;

fs = figureStyle('Thesis');

%% Specified Length
Li = 24*12;

% Calculations
PnxLi  = wf.Pnx(K,Li);
PnyLi  = wf.Pny(K,Li);
PnzLi  = wf.Pnz(K,Li);
PncaLi = wf.Pnca(K,Li);

wf.beta_T(250,Li);
wf.beta_Tb(250,Li);
[Pnys,tau,x] = wf.Pnys(250,Li);

P = linspace(0,Py,100);

beta_T  = nan(size(P));
beta_Tb = nan(size(P));
Pnys    = nan(size(P));
tau     = nan(size(P));
x       = nan(size(P));

for i = 1:length(P)
    beta_T(i)  = wf.beta_T(P(i),Li);
    beta_Tb(i) = wf.beta_Tb(P(i),Li);
    [iPnys,itau,ix] = wf.Pnys(P(i),Li);
    Pnys(i) = iPnys;
    tau(i) = itau;
    x(i) = ix;
end


% Figure 3
hf = fs.figure(3.25,2.0);
ha = fs.axes();
plot([0 Py],[0 Py],'k--','LineWidth',1)
plot(P,Pnys,'k-','LineWidth',2)
xlabel('Applied Axial Load, P (kips)')
ylabel('P_{ny}* (kips)')
%axis equal
xlim([0 Py])
ylim([0 1.1*PnyLi])
saveas(hf,'Figure_03');
print(hf,'-dsvg','Figure_03');


% Figure 4
hf = fs.figure(3.25,3.25);
ha = fs.axes([0.12 0.12 0.84 0.84]);
plot(PnxLi*[1 1],[0 max(beta_Tb)],'k--','LineWidth',2)
plot(PnyLi*[1 1],[0 max(beta_Tb)],'k--','LineWidth',2)
plot(PnzLi*[1 1],[0 max(beta_Tb)],'k--','LineWidth',2)
plot(PncaLi*[1 1],[0 max(beta_Tb)],'k--','LineWidth',2)
plot(P,beta_T,'k--','LineWidth',2)
plot(P,beta_Tb,'k-','LineWidth',2)
xlabel('Axial Load (kips)')
ylabel('Required Stiffness (kip-in./rad/in.)')
xlim([0 450])
ylim([0 70])
saveas(hf,'Figure_04');
print(hf,'-dsvg','Figure_04');


%% Varying Length
L = linspace(0,50*12,100);

% Calculations
Pnx             = nan(size(L));
Pny             = nan(size(L));
Pnz             = nan(size(L));
Pnca            = nan(size(L));
Pnca_0          = nan(size(L));
Pnca_A          = nan(size(L));
Pnca_B          = nan(size(L));
beta_Tb_Pnx     = nan(size(L));
beta_Tb_Pnz     = nan(size(L));

Pnx(1)          = wf.Pno_lb;
Pny(1)          = Pnx(1);
Pnz(1)          = Pnx(1);
Pnca(1)         = Pnx(1);
Pnca_0(1)       = Pnx(1);
Pnca_A(1)       = Pnx(1);
Pnca_B(1)       = Pnx(1);

for i = 2:length(L)
    Pnx(i)       = wf.Pnx(L(i),K);
    Pny(i)       = wf.Pny(L(i),K);
    Pnz(i)       = wf.Pnz(L(i),K);
    Pnca(i)      = wf.Pnca(L(i),K);
    Pnca_0(i)    = wf.Pr_given_betaTb(L(i),0);
    Pnca_A(i)    = wf.Pr_given_betaTb_and_phiMn(L(i),beta_A,phiMn_A);
    Pnca_B(i)    = wf.Pr_given_betaTb_and_phiMn(L(i),beta_B,phiMn_B);
    
    if Pnx(i) > Pnca(i)
        beta_Tb_Pnx(i) = wf.beta_Tb(Pnx(i),L(i));
    else
        beta_Tb_Pnx(i) = 0;
    end
      
    if Pnz(i) > Pnca(i)
        beta_Tb_Pnz(i) = wf.beta_Tb(Pnz(i),L(i));
    else
        beta_Tb_Pnz(i) = 0;
    end
end


% Figure 2
hf = fs.figure(3.25,3.25);
ha = fs.axes([0.15 0.15 0.80 0.80]);
plot([0 Li/12 Li/12],[PnxLi PnxLi 0],'k--','LineWidth',1)
plot([0 Li/12],[PnyLi PnyLi],'k--','LineWidth',1)
plot([0 Li/12],[PnzLi PnzLi],'k--','LineWidth',1)
plot([0 Li/12],[PncaLi PncaLi],'k--','LineWidth',1)
plot(L/12,Pnx,'k-','LineWidth',1)
plot(L/12,Pny,'k-','LineWidth',1)
plot(L/12,Pnz,'k-','LineWidth',1)
plot(L/12,Pnca,'k-','LineWidth',1)
plot(Li/12,PnxLi,'ok','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',4)
plot(Li/12,PnyLi,'ok','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',4)
plot(Li/12,PnzLi,'ok','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',4)
plot(Li/12,PncaLi,'ok','LineWidth',1,'MarkerFaceColor',[1 1 1],'MarkerSize',4)
xlabel('Length (ft)')
ylabel('Axial Capacity (kips)')
saveas(hf,'Figure_02');
print(hf,'-dsvg','Figure_02');


% Figure 7
hf = fs.figure(3.25,2.0);
ha = fs.axes([0.15 0.18 0.80 0.77]);
plot(L/12,0.9*Pnca,'k--','LineWidth',2)
plot(L/12,Pnca_0,'k-','LineWidth',2)
plot(L/12,Pnca_A,'k-','LineWidth',2)
plot(L/12,Pnca_B,'k-','LineWidth',2)
xlabel('Length (ft)')
ylabel('Axial Capacity (kips)')
saveas(hf,'Figure_07');
print(hf,'-dsvg','Figure_07');


% Figure 8
hf = fs.figure(3.25,2.0);
ha = fs.axes([0.15 0.18 0.80 0.77]);
plot(L/12,beta_Tb_Pnx,'k-','LineWidth',2)
plot(L/12,beta_Tb_Pnz,'k-','LineWidth',2)
xlabel('Length (ft)')
ylabel('Required Stiffness (kip-in./rad/in.)    .')
saveas(hf,'Figure_08');
print(hf,'-dsvg','Figure_08');
