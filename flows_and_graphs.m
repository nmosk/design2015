%ULTIMATE GRAPH MAKER SCRIPT 2015 
clear;
clc;
close all;
PFR_latest

global MW_eb MW_st MW_b MW_h MW_en MW_t MW_me P_st n_st time MW_w
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
MR = 8;

xx = [0.01:0.01:0.97]';
ss = xs_interp(xx,X,s_123);

[flow_eb,flow_h2,flow_t,flow_me,flow_b,flow_en] = flows(ss(:,1),ss(:,2));

recycle_eb = n_st*MW_eb./ss(:,1).*(1-xx)./xx;
n_st = n_st*ones(size(recycle_eb));
flow_steam = MR.*(recycle_eb + flow_eb).*MW_w./MW_eb;
flow_out_reactor = recycle_eb+flow_h2+flow_me+flow_b+flow_en+flow_t+P_st+flow_steam;
moles_out_reactor = [(recycle_eb./MW_eb + flow_h2./MW_h + flow_me./MW_me...
    + flow_b./MW_b + flow_en./MW_en + flow_t./MW_t + n_st...
    + flow_steam./MW_w)];

y_eb = recycle_eb/MW_eb./moles_out_reactor;
y_h = flow_h2/MW_h./moles_out_reactor;
y_me = flow_me/MW_me./moles_out_reactor;
y_b = flow_b/MW_b./moles_out_reactor;
y_en = flow_en/MW_en./moles_out_reactor;
y_t = flow_t/MW_t./moles_out_reactor;
y_st = n_st./moles_out_reactor;
y_w = flow_steam/MW_w./moles_out_reactor;


figure(1)
plot(xx,y_eb,xx,y_h,'-.',xx,y_b,'-..',xx,y_st,'--',xx,y_t,'-..k',xx,y_w,'LineWidth',3)
legend('y_e_B','y_H_2','y_B/y_e_N','y_S_t','y_T/y_M_e','y_H_2_O')
xlabel('conversion, x')
ylabel('mole fraction, y')
title('from reactor to separation')
set(gca,'FontSize',26)

figure(2)
plot(xx,recycle_eb./time,'LineWidth',3)
xlabel('conversion, x')
ylabel('recycled eb [kg/hr]')
set(gca,'FontSize',26)

figure(3)
plot(xx,flow_eb./time)
xlabel('conversion, x')
ylabel('Fresh feed EB, [kg/hr]')
set(gca,'FontSize',26)

figure(4)
plot(xx,flow_out_reactor./time,'LineWidth',3)
xlabel('conversion, x')
ylabel('Flow to separation [kg/hr]')

figure(5)
plot(xx,(flow_eb+recycle_eb+flow_steam)./time,'LineWidth',3)
xlabel('conversion, x')
ylabel('Feed to reactor [kg/hr]')
set(gca,'FontSize',26)

set_x = 75

format short

inlet_reactor_kg_hr = (flow_eb(set_x)+ recycle_eb(set_x) + flow_steam(set_x))./time
outlet_reactor_kg_hr = flow_out_reactor(set_x)./time
recycle_eb_kg_hr = recycle_eb(set_x)./time
moles_out_reactor_hr = moles_out_reactor(set_x)./time
mole_frac_eb = y_eb(set_x)
mole_frac_st = y_st(set_x)
mole_frac_water = y_w(set_x)
mole_frac_h2 = y_h(set_x)
mole_frac_tol = y_t(set_x)
mole_frac_benz = y_b(set_x)

%distillation
disp('Distillation info:')
basis_molar_dist = (mole_frac_st + mole_frac_eb + mole_frac_benz + mole_frac_tol);
moles_hr_into_distillation = basis_molar_dist* moles_out_reactor(set_x)./time

dist_frac_tol = mole_frac_tol/basis_molar_dist
dist_frac_eb = mole_frac_eb/basis_molar_dist
dist_frac_st = mole_frac_st/basis_molar_dist
dist_frac_benz = mole_frac_benz/basis_molar_dist

