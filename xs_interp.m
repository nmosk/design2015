function [ss_outputs] = xs_interp(x_want, X, s_123)

%x_want = x values you want to find s values at


ss_outputs=zeros(length(x_want),3);
for j=1:length(x_want)
    for i=1:3
    ss_outputs(j,i)=spline(X,s_123(:,i),x_want(j)); %a is cycle_tinit values
    end
end
ss_outputs; %s values at the x_want values

end