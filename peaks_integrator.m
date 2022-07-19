function [integral,peaks_coordinate,back_step,integ_length] = peaks_integrator(vector,peaks_coordinate,integration_method,RB_type_override)
back_step=integration_method{1}; % set value or (lower) limit value (depending of the integration mode)
integ_length=integration_method{2}; % set value or (higher) limit value (depending of the integration mode)
precise_extrema=integration_method{3}; % precise maxima on or off
LB_type=integration_method{4}; %set or deriv.
RB_type=integration_method{5}; % length (RB-LB) set, mid, dmid, deriv. or auto-shrinked
background_type=integration_method{6}; % constant, linear or constant_R
peak_symbol=integration_method{7};
LB_drv_trshld=integration_method{8};
RB_drv_trshld=integration_method{9};

if strcmp(RB_type_override,'no')==1
elseif strcmp(RB_type_override,'mid')==1 | strcmp(RB_type_override,'dmid')==1  | strcmp(RB_type_override,'auto-shrinked')==1
 RB_type= RB_type_override;   
end


back_step_max=round(mean(back_step));
integ_length_max=round(mean(integ_length));
for kk = 1:size(peaks_coordinate,2)

%% locate precisely the extrema 
if strcmp(precise_extrema,'yes')==1
    if peak_symbol==1
fine_extrema=max(vector(peaks_coordinate(kk)-20:peaks_coordinate(kk)+20,2));
fine_extrema_position=find(vector(peaks_coordinate(kk)-20:peaks_coordinate(kk)+20,2)==fine_extrema);
elseif peak_symbol==-1
fine_extrema=min(vector(peaks_coordinate(kk)-20:peaks_coordinate(kk)+20,2));
fine_extrema_position=find(vector(peaks_coordinate(kk)-20:peaks_coordinate(kk)+20,2)==fine_extrema);  
end
peaks_coordinate(kk)=peaks_coordinate(kk)-21+fine_extrema_position;
end

%% locate the left bound
if strcmp(LB_type,'deriv.')==1
for ii=1:back_step_max
local_slope=(vector(peaks_coordinate(kk)-ii,2)-vector(peaks_coordinate(kk)-ii-1,2))/(vector(peaks_coordinate(kk)-ii,1)-vector(peaks_coordinate(kk)-ii-1,1));    
if peak_symbol==-1
local_slope=-local_slope;   
end
if local_slope <=LB_drv_trshld 
back_step(kk)=ii;
break
end
end
elseif strcmp(LB_type,'set')==1 | strcmp(LB_type,'set (ref)')==1
 back_step(kk)=back_step_max;   
end

%% determine the integration length
if strcmp(RB_type,'mid')==1
integ_length(kk)=back_step(kk);
elseif strcmp(RB_type,'dmid')==1
integ_length(kk)=2*back_step(kk);
elseif strcmp(RB_type,'set')==1 | strcmp(RB_type,'set (ref)')==1
integ_length(kk)=integ_length_max;
elseif strcmp(RB_type,'deriv.')==1 | strcmp(RB_type,'auto-shrinked')==1
    for ii=5:integ_length_max-back_step(kk)
local_slope=(vector(peaks_coordinate(kk)+ii,2)-vector(peaks_coordinate(kk)+ii-1,2))/(vector(peaks_coordinate(kk)+ii,1)-vector(peaks_coordinate(kk)+ii-1,1));    
if peak_symbol==-1
local_slope=-local_slope;   
end
if local_slope >=-RB_drv_trshld 
integ_length(kk)=ii+back_step(kk);
break
else integ_length(kk)=integ_length_max;
end
    end
if strcmp(RB_type,'auto-shrinked')==1    
 integ_length(kk)=round(integ_length(kk)/2+back_step(kk));   
end   
end

%% integration!
    left_base_coordinate=peaks_coordinate(kk)-back_step(kk);
    right_base_coordinate=peaks_coordinate(kk)-back_step(kk)+integ_length(kk);
 integ=trapz(vector(left_base_coordinate:right_base_coordinate,1),vector(left_base_coordinate:right_base_coordinate,2));
% %if trapz is not available, replace the previous line by:
%  integ=0;
% for nn=left_base_coordinate:1:right_base_coordinate-1
% integ=integ+(vector(nn,2)+vector(nn+1,2))/2*(vector(nn+1,1)-vector(nn,1));
% end

if strcmp(background_type,'linear')==1
 background=(vector(left_base_coordinate,2)+vector(right_base_coordinate,2))/2*(vector(right_base_coordinate,1)-vector(left_base_coordinate,1));
elseif strcmp(background_type,'constant')==1 %| strcmp(background_type,'auto_mid')==1
     background=vector(left_base_coordinate,2)*(vector(right_base_coordinate,1)-vector(left_base_coordinate,1));
     elseif strcmp(background_type,'constant_R')==1
     background=vector(right_base_coordinate,2)*(vector(right_base_coordinate,1)-vector(left_base_coordinate,1));
end

integral(kk)=integ-background;

end
end