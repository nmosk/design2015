



TT = 10:200;
plot(TT,vap_P(TT,'eb')./vap_P(TT,'st'),TT,vap_P(TT,'to')./vap_P(TT,'eb'),TT,vap_P(TT,'be')./vap_P(TT,'to'))
ylabel('relative volatilities')
xlabel('Temperature [c]')
legend('eb-st','to-eb','be-to')
