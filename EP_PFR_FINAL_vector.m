% DESCRIPTION -------------------------------------------
% Economic Potential calculation

% Is used in the mainscript of the PFR code

% NOTE: most values are outputted in arrays
% -------------------------------------------------------

% ABBREVIATIONS of SPECIES ------------------------------
% ethylbenzene -> eB 
% styrene -> St 
% benzene -> B 
% toulene -> T 
% ethylene -> eN 
% hydrogen -> H2 
% methane -> Me 
% water-> H20
% -------------------------------------------------------

% NEW ORDER OF SPECIES IN ARRAYS ------------------------
% eb  -> 1
% St  -> 2
% B  ->  3
% T  ->  4
% eN  -> 5
% H2  -> 6
% Me  -> 7
% H20 -> 8
% -------------------------------------------------------

function [EP,EPF,WC, cost_eb, value_sty, value_by] = EP_PFR_FINAL_vector(s_123,MR,X)

% PRICES OF SPECIES -------------------------------------
% prices in dollars -- eB, St, B, T, en, H2 [$/kg]
price_eb = 1.1;
price_st = 1.37;
price_b = 0.86;
price_t = 0.97;
price_en = 1.2;
price_h = 2;
price_fuel = 3; % [$/MMBTU]
price_water = 0.06/1000; % [$/kg]

price=[price_eb price_st price_b price_t price_en price_h price_fuel price_water]; %price_fuel in [$/MMBTU], rest in [$/kg]

%fuel prices [$/kg fuel]
fuel_b = 0.12;
fuel_t = 0.12;
fuel_en = 0.14;
fuel_h = 0.40;
fuel_me = 0.16; 

price_fuel=[0 0 fuel_b fuel_t fuel_en fuel_h fuel_me 0]; %[$/kg fuel]
% -------------------------------------------------------

%m MOLECULAR WEIGHTS [g/mol]
% -------------------------------------------------------
MW_eb = 106.17;
MW_st = 104.18;
MW_b = 78.11;
MW_h = 2.016;
MW_en = 28.05;
MW_t = 92.138;
MW_me = 16.04;

MW = [106.17 104.18 78.1118 92.1384 28.0532 2.016 16.0425 18.016]; %[g/mol]
% -------------------------------------------------------

% OTHER -------------------------------------------------
time = 8400; % [hrs/yr]

% Product of Styrene
P_st = 100 * 10^6; % [kg/yr]

% Converting product of styrene from kg/yr to kmol/yr
n_st = P_st/MW(2); % [kmol/yr]
% -------------------------------------------------------

%%

% CALCULATING ECONOMIC POTENTIAL ------------------------
% Calculates economic potential if all byproducts sold as chemical species
% NOTE %%%***********************************************************************************
% selectivity spans 0 to 1, s_123(i,1) is selectivity of styrene reaction, s_123(i,2) is
% selectivity of benz reaction
% *******************************************************************************************

% Initializes EP matrix
EP=zeros(length(s_123),1);

    % Separation cost model
        sep1 = 1.9e6 ; % cost of column 1
        sep2 = 1.6e6; %  cost of column 2
        sep3 = 1.6e6; %  cost of column 3
%- 6.36e3 - 3.3e3 - sep1 - sep2 - sep3
    EP = (price_st*P_st - (n_st*MW_eb./s_123(:,1)).*price_eb - (MR*n_st*MW(8)*price_water) + (n_st*MW_en*s_123(:,2)./s_123(:,1)).*fuel_en ...
        + (n_st*MW_b*s_123(:,2)./s_123(:,1)).*price_b + (n_st*MW_t)*(s_123(:,3))./s_123(:,1)*price_t ...
        + (n_st*MW_me)*(s_123(:,3))./s_123(:,1)*fuel_me + (n_st*MW_h)*(1-s_123(:,2))./s_123(:,1).*price_h);
    
% FIGURE 1: EP vs s_123
% Plots economic potential versus selectivity of styrene (s_123(:,1))
figure
hold on
title('Annual economic potential (chemical by product) vs selectivity')
ylabel('EP, $/yr')
xlabel('s_1')
plot(s_123(:,1),EP,'LineWidth',2)
axis([0 1 0 7*10^8])
set(gca,'FontSize',26)
% -------------------------------------------------------


% CALCULATING ECONOMIC POTENTIAL of FUEL ----------------
% Calculates economic potential if all byproducts sold as fuel
% Calculates working capital (WC) for 2 months
% Initializes EP_fuel matrix

%EP_fuel=zeros(length(s_123),1);
%WC=EP_fuel; %working cost
% for i=1:length(s_123)
% 
%     EP_fuel(i) = (price_st*P_st - (n_st*MW_eb./s_123(i,1)).*price_eb - (MR*n_st*MW(8)*price_water) + (n_st*MW_en*s_123(i,2)./s_123(i,1)).*fuel_en ...
%         + (n_st*MW_b*s_123(i,2)./s_123(i,1)).*fuel_b + (n_st*MW_t)*(s_123(i,3))./s_123(i,1)*fuel_t ...
%         + (n_st*MW_me)*(s_123(i,3))./s_123(i,1)*fuel_me + (n_st*MW_h)*(1-s_123(i,2))./s_123(i,1).*fuel_h-(5.*10.^5/X(i)));
%     WC(i)= ((n_st*MW_eb./s_123(i,1)).*price_eb - (MR*n_st*MW(8)*price_water))/12*2 ; % cost for 2 months
% end

%     EPF = (price_st*P_st - (n_st*MW_eb./s_123(:,1)).*price_eb - (MR*n_st*MW(8)*price_water) +(n_st*MW_en*s_123(:,2)./s_123(:,1)).*fuel_en ...
%         + (n_st*MW_b*s_123(:,2)./s_123(:,1)).*fuel_b + (n_st*MW_t)*(s_123(:,3))./s_123(:,1)*fuel_t ...
%         + (n_st*MW_me)*(s_123(:,3))./s_123(:,1)*fuel_me + (n_st*MW_h)*(1-s_123(:,2))./s_123(:,1).*fuel_h);


% CALCULATING ECONOMIC POTENTIAL of FUEL ----------------
% Calculates economic potential if all byproducts sold as fuel
% Calculates working capital (WC) for 2 months
% Initializes EP_fuel matrix

   cost_eb= (n_st*MW_eb./s_123(:,1)).*price_eb ;
   value_sty= price_st*P_st;
   value_by=((n_st*MW_en*s_123(:,2)./s_123(:,1)).*fuel_en ...
        + (n_st*MW_b*s_123(:,2)./s_123(:,1)).*price_b + (n_st*MW_t)*(s_123(:,3))./s_123(:,1)*price_t ...
        + (n_st*MW_me)*(s_123(:,3))./s_123(:,1)*fuel_me + (n_st*MW_h)*(1-s_123(:,2))./s_123(:,1).*fuel_h)- (MR*n_st*MW(8)*price_water) ;
   
    EPF= value_sty - cost_eb + value_by ;
    
    % Working Capital for 2 months
    WC= ((n_st*MW_eb./s_123(:,1)).*price_eb - (MR*n_st*MW(8)*price_water))/12*2 ; % cost for 2 months

% FIGURE 2: EP_fuel vs s_123
% Plots EP_fuel versus selectivity of styrene (s_123(:,1))
figure
hold on
title('Annual economic potential (fuel by product) vs selectivity')
ylabel('EPF, $/yr')
xlabel('s_1')
plot(s_123(:,1),EPF,'LineWidth',2)
axis([0 1 0 7*10^8])
set(gca,'FontSize',26)

            
end


