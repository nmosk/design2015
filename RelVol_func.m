% relative volatilities (relative to styrene!)
% use with vap_P.m

function [RelVol] = RelVol_func(T)

T = T-273.15; % [C]
RelVol_be = vap_P(T,'be')/vap_P(T,'st')
RelVol_to = vap_P(T,'to')/vap_P(T,'st')
RelVol_eb = vap_P(T,'eb')/vap_P(T,'st')
RelVol = [RelVol_be; RelVol_to; RelVol_eb]

TT = 10:200;
plot(TT,vap_P(TT,'be')./vap_P(TT,'st'),TT,vap_P(TT,'to')./vap_P(TT,'st'),TT,vap_P(TT,'eb')./vap_P(TT,'st'))
ylabel('relative volatilities (to styrene)')
xlabel('Temperature [c]')
legend('benzene','toluene','ethylbenzene')

end