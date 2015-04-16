%Economic potential level 2 calculation
% function of selectivity
%1 -> eB ; 2 -> St ; 3 -> B ; 4 -> T ; 5 -> eN ; 6 -> H2 ; 7 -> Me ; 8 ->
%H20

function [EP,EPF,WC, cost_eb, value_sty, value_by] = EP_PFR_FINAL_vector(s_123,MR,X)
%prices in dollars -- eB, St, B, T, en, H2 [$/kg]

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

%molecular weights [g/mol]
MW_eb = 106.17;
MW_st = 104.18;
MW_b = 78.11;
MW_h = 2.016;
MW_en = 28.05;
MW_t = 92.138;
MW_me = 16.04;

MW = [106.17 104.18 78.1118 92.1384 28.0532 2.016 16.0425 18.016]; %[g/mol]

time = 8400; % [hrs/yr]

P_st = 100 * 10^6; % [kg/yr]
n_st = P_st/MW(2); % [kmol/yr]

% selectivity spans 0 to 1, s_123(i,1) is selectivity of styrene reaction, s_123(i,2) is
% selectivity of benz reaction

EP=zeros(length(s_123),1);
%for i=1:length(s_123)

    EP = (price_st*P_st - (n_st*MW_eb./s_123(:,1)).*price_eb - (MR*n_st*MW(8)*price_water) +(n_st*MW_en*s_123(:,2)./s_123(:,1)).*price_en ...
        + (n_st*MW_b*s_123(:,2)./s_123(:,1)).*price_b + (n_st*MW_t)*(s_123(:,3))./s_123(:,1)*price_t ...
        + (n_st*MW_me)*(s_123(:,3))./s_123(:,1)*fuel_me + (n_st*MW_h)*(1-s_123(:,2))./s_123(:,1).*price_h-(5.*10.^5./X));
%end
figure
hold on
title('Annual economic potential (chemical by product) vs selectivity')
ylabel('EP, $/yr')
xlabel('s_1')
plot(s_123(:,1),EP,'LineWidth',2)
axis([0 1 0 7*10^8])
set(gca,'FontSize',26)

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
   cost_eb= (n_st*MW_eb./s_123(:,1)).*price_eb ;
   value_sty= price_st*P_st;
   value_by=((n_st*MW_en*s_123(:,2)./s_123(:,1)).*fuel_en ...
        + (n_st*MW_b*s_123(:,2)./s_123(:,1)).*fuel_b + (n_st*MW_t)*(s_123(:,3))./s_123(:,1)*fuel_t ...
        + (n_st*MW_me)*(s_123(:,3))./s_123(:,1)*fuel_me + (n_st*MW_h)*(1-s_123(:,2))./s_123(:,1).*fuel_h)- (MR*n_st*MW(8)*price_water)+(n_st*MW_en*s_123(:,2)./s_123(:,1)).*fuel_en ;
   
    EPF= value_sty - cost_eb + value_by -(5.*10.^5./X);

    WC= ((n_st*MW_eb./s_123(:,1)).*price_eb - (MR*n_st*MW(8)*price_water))/12*2 ; % cost for 2 months

    
figure
hold on
title('Annual economic potential (fuel by product) vs selectivity')
ylabel('EPF, $/yr')
xlabel('s_1')
plot(s_123(:,1),EPF,'LineWidth',2)
axis([0 1 0 7*10^8])
set(gca,'FontSize',26)

            
end


