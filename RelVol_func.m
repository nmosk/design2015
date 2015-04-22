% relative volatilities (relative to styrene!)
% use with vap_P.m

function [RelVol] = RelVol_func(T)

T = T-273.15; % [C]
RelVol_be_to = vap_P(T,'be')/vap_P(T,'to')
RelVol_to_eb = vap_P(T,'to')/vap_P(T,'eb')
RelVol_eb_st = vap_P(T,'eb')/vap_P(T,'st')
RelVol = [RelVol_be_to; RelVol_to_eb; RelVol_eb_st]

TT = 10:200;
figure
plot(TT,vap_P(TT,'eb')./vap_P(TT,'st'),TT,vap_P(TT,'to')./vap_P(TT,'eb'),TT,vap_P(TT,'be')./vap_P(TT,'to'))
ylabel('relative volatilities')
xlabel('Temperature [c]')
legend('eb-st','to-eb','be-to')

end