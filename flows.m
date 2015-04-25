
function [flow_eb,flow_h2,flow_t,flow_me,flow_b,flow_en,val_eb,val_t,val_b,val_fuel] = flows(s1,s2)

s2;

global MW_eb MW_st MW_b MW_h MW_en MW_t MW_me time P_st n_st MW_w
%molecular weights [g/mol]
MW_eb = 106.17;
MW_st = 104.18;
MW_b = 78.11;
MW_h = 2.016;
MW_en = 28.05;
MW_t = 92.138;
MW_me = 16.04;
time = 8400; % [hrs/yr]
P_st = 100*10^6; % [kg/yr]
n_st = P_st/MW_st; % [kmol/yr]
MW_w = 18; 

% selectivity spans 0 to 1, s1 is selectivity of styrene reaction, s2 is
% selectivity of benz reaction

disp('Flow rates of Level 2 streams [kg/yr]:')    

% level 2
flow_eb = (n_st*MW_eb./s1);
flow_h2 = (n_st*MW_h)*(1-s2)./s1;
flow_t = (n_st*MW_t)*(1-s1-s2)./s1;
flow_me = (n_st*MW_me)*(1-s1-s2)./s1;
flow_b = (n_st*MW_b*s2./s1);
flow_en = (n_st*MW_en*s2./s1);

% value of product per year 
price_eb = 1.1;
price_st = 1.37;
price_b = 0.86;
price_t = 0.97;
price_en = 1.2;
price_h = 2;
price_fuel = 3; % [$/MMBTU]
price_water = 0.06/1000; % [$/kg]

% calcs the value of the product in a year for conceptual design
val_eb=flow_eb*price_eb;
val_t=flow_t*price_t;
val_b=flow_b*price_b;
val_fuel=(flow_me+flow_h2+flow_en).*price_fuel;
end



