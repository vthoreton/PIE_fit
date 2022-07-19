function[res]= fitPIE(p)
global fitting_data %fitting_param
% weighting mass 34:
weight_34=1;
% getting the data from fitting_data
T=fitting_data{1};
R0=fitting_data{2};
f18=fitting_data{3};
f32=fitting_data{4};
f34=fitting_data{5};
f36=fitting_data{6};
f18_RT=fitting_data{7};
f36_RT=fitting_data{8};
S=fitting_data{9};
molar_flow_rate=fitting_data{10};
low=fitting_data{11};
high=fitting_data{12};
%f36c=(1-p).^2*f18_RT^2.*exp(-2.*R0*S/molar_flow_rate)./(1-2.*p)+f36_RT.*exp(-R0.*S./(molar_flow_rate.*p))-(1-p).^2*f18_RT^2.*exp(-R0.*S./(molar_flow_rate.*p))./(1-2.*p);
f36c=(1-p).^2*f18_RT^2.*exp(-2.*R0*S/molar_flow_rate)./(1-2.*p)+f36_RT.*exp(-R0.*S./(molar_flow_rate.*p))-(1-p).^2*f18_RT^2.*exp(-R0.*S./(molar_flow_rate.*p))./(1-2.*p);
f34c=2*(f18-f36c);
f32c=1-f34c-f36c;
res = sum((minus(f36(low:high),f36c(low:high))).^2)+weight_34*sum((minus(f34(low:high),f34c(low:high))).^2)+sum((minus(f32(low:high),f32c(low:high))).^2);
end