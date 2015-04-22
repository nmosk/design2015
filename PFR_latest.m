clear;
clc;
close all;

% test
global T MR r R Peb0 Ptot

T = 600+273; % kelvin
MR = 8;   % steam/feed
r = 1.987; % cal/mol K
R = 8.314 * 10^-5; % m3 bar/K mol
Ptot = 2.02;

Peb0 = 1/(1+MR)* Ptot ;
Feb=6; % initial guess for our inlet
Fsteam=Feb*MR*Ptot;

rho = Ptot/R/T;

global k1 kn1 k2 k3

%rate constants
k1 = 1.177 * 10^8 * exp(-21708/r/T) ;
kn1 = 20.965 * exp(7804/r/T);
k2 = 9.206 * 10^12 * exp(-45675/r/T) ;
k3 = 4.724 * 10^7 * exp(-18857/r/T) ;

MW_st = 104; % g/mol
Pst = 100 * 10^6 / MW_st * 1000; % mol/yr
hr_per_yr = 8400; %hr/yr
s_per_hr = 3600;

N = 100 ;
Ceb = zeros(N,1);
Cst = Ceb; % mol/m^3
Ch2 = Ceb; % mol/m^3
Cb=Ceb;
Ct=Ceb;
X = Ceb;
s = Ceb;
Fin = Ceb;
Fin_tot = Ceb;
q = Ceb;
V = Ceb;


Ceb(1) = (1/(1+MR)) * rho;

delta_Tau = 0.00005;

for n = 2:10000
    % CALCULATING CONCENTRATIONS OF EB , ST , H2 ---------------
    Ceb(n) = Ceb(n-1) + delta_Tau * (kn1*R*R*T*T*Cst(n-1)*Ch2(n-1) - k1*R*T*Ceb(n-1) - k2*R*T*Ceb(n-1) - k3*R*R*T*T*Ceb(n-1)*Ch2(n-1));
    Cst(n) = Cst(n-1) + delta_Tau * (k1*R*T*Ceb(n-1) - kn1*R*R*T*T*Cst(n-1)*Ch2(n-1));
    Ch2(n) = Ch2(n-1) + delta_Tau * (k1*R*T*Ceb(n-1) - k3*R*R*T*T*Ceb(n-1)*Ch2(n-1) - kn1*R*R*T*T*Cst(n-1)*Ch2(n-1));
    Cb (n) = Cb (n-1) + delta_Tau * (k2*R*T*Ceb(n-1));
    Ct (n) = Ct (n-1) + delta_Tau * (k3*R*R*T*T*Ceb(n-1)*Ch2(n-1));
    
    % calculating conversion of eB
    X(n) = (Ceb(1) - Ceb(n))./Ceb(1);
    % calculating selectivity of styrene
    s(n) = Cst(n)./(Ceb(1) - Ceb(n));
    
    % inlet flowrate
    Fin(n) = Pst./X(n)./s(n) / hr_per_yr;   % mol/hr
    
    % total flowrate
    Fin_tot(n) = Fin(n) * (1 + MR);            % mol/hr
    
    % volumetric flowrate
    q(n) = Fin_tot(n) / rho; % m3/hr
    % volume of reactor in m^3
    V(n) = delta_Tau / s_per_hr * n * q(n);
    %%
    % CALCULATING SELECTIVITIES OF ALL RXNS -------------------
    
    % For the matrix:
    % selectivity of rxn 1 -> 1
    % selectivity of rxn 2 -> 2 (benzene and ethylene)
    % selectivity of rxn 3 -> 3 (toluene and methane)
    
    s_123(n,:) = [s(n) Cb(n)./(Ceb(1) - Ceb(n)) Ct(n)./(Ceb(1) - Ceb(n))];
end

figure
plot(X,s)
xlabel('Conversion, X')
ylabel('Selectivity, S')
title('Selectivity T = 600C, MR = 8')

figure
plot(X,V)
xlabel('Conversion, X')
ylabel('Volume, V [m^3]')
title('volume, T = 600C, MR = 8')

%  Cb=Cb(1:200:end)   ;
%   Ceb=Ceb(1:200:end)   ;
%     Ch2=Ch2(1:200:end) ;
%      Cst=Cst(1:200:end) ;
%       Ct=Ct(1:200:end)   ;
%  Fin=Fin(1:200:end)  ;
%   Fin_tot=Fin_tot(1:200:end)   ;
%  q=q(1:200:end)   ;
%
%  s=s(1:200:end)   ;
%   s_123=s_123(1:200:end,:)   ;
%
%  V=V(1:200:end)   ;
%  X=X(1:200:end)  ;



%%
% ---------------------------------------------------------

% CALCULATING ECONOMIC POTENTIAL --------------------------
% Inputs selectivity, molar ratio, and conversion
% Outputs economic potential (EP), economic potential if all byproducts
% sold as fuel  (EP_fuel), and working capital (WC)


[EP,EPF,WC, priceeb, pricest, priceby] = EP_PFR_FINAL_vector(s_123,MR,X);

% ---------------------------------------------------------

% ECONOMIC ANALYSIS ---------------------------------------
% Inputs V, WC, EP, X
% Outputs ROI_BT, TCI, NPVs, etc
tic
[H_E,ROI_BT, reac, V_ft, D_fact, WC_CF, PO_CF , TCI, H, D, FC,TI, SU, WC, Profit_BT, Profit_AT, C_F, Cashflow_d, Bond_Fin, D_CF, NPV_0, NPV_proj,NPV_percent,Depreciation] = conceptual_econ_new(V, WC, EP,X);

time=toc;
% ---------------------------------------------------------
%%
% FIGURE 1: s vs X
% Plots all the different MRs in one plot
figure
plot(X,s,'b-','LineWidth',2)
xlabel('Conversion, X','FontSize',26)
ylabel('Selectivity, S','FontSize',26)
axis([0 1 0 1])
title('PFR; S vs X;  T = 793k')
legend('MR=8','FontSize',20)
set(gca,'FontSize',26)
line([0.75 0.75], ylim,'Color','r','LineStyle','--','LineWidth',2);


% FIGURE 2: V vs X
% Plots volume against conversion
figure
plot(X,V,'b-','LineWidth',2)
xlabel('Conversion, X','FontSize',26)
ylabel('Volume (m^3)','FontSize',26)
% axis([0 1 0 10])
title('PFR; V vs X;  T = 600C')
legend('MR=8','FontSize',20)
line([0.75 0.75], ylim,'Color','r','LineStyle','--','LineWidth',2);

set(gca,'FontSize',26)

%FIGURE 3: EP vs V
%Plots all the economic potentials against corresponding volumes
figure
	plot(V,EP,'b-','LineWidth',2)
        axis([0 2 0 3*10^7])
        ylabel('EP  [$/yr]','FontSize',26)
        xlabel('Volume [m^3]','FontSize',26)
        set(gca,'FontSize',26)




% FIGURE 4: NPV_proj vs X and TCI vs X
% Plots NPV_proj versus conversion and TCi versus conversion
figure
plot(X,NPV_proj,'g-',X,TCI,'b-','LineWidth',3);
legend('Net Present Value','Total Capital Investment') % label left y-axis
ylabel('Dollars [$]')
xlabel('Conversion, X','FontSize',26)
axis([0.4 0.9 0 1*10^8])
set(gca,'FontSize',26)
line([0.75 0.75], ylim,'Color','r','LineStyle','--','LineWidth',2);


% FIGURE 5: NPV_percent vs X and ROI_BT vs X
% Plots NPV_percent versus conversion and ROI_BT versus conversion



figure
plot(X,NPV_percent,'g-',X,ROI_BT*100,'b-','LineWidth',3);
ylabel('Percent per Year [%/yr]')
legend('Net Present Value' ,'Return on Investment') % label left y-axis
xlabel('Conversion, X','FontSize',26)
axis([0 0.9 0 40])
set(gca,'FontSize',26)
line([0.75 0.75], ylim,'Color','r','LineStyle','--','LineWidth',2);


figure
	plot(s,EP,'b-','LineWidth',2)
        axis([0.4732 1 0 3*10^7])
        ylabel('EP  [$/yr]','FontSize',26)
        xlabel('Selectivity, S','FontSize',26)
        line([0.75 0.75], ylim,'Color','r','LineStyle','--','LineWidth',2);

        set(gca,'FontSize',26)
