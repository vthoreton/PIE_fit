function [positions] = peaks_finder(vector,resolution,modifier,low_treshold,high_treshold)
section_min=1;
index=1;
positions=[];
while section_min < size(vector,1)
         if  section_min+resolution <=size(vector,1)
        section_max=section_min+resolution;
    else
        section_max=size(vector,1);
         end              
local_min_log=min(log10(vector(section_min:section_max,1)));
treshold_level_L=10^(low_treshold+local_min_log);
treshold_level_H=10^(high_treshold+local_min_log);         
        max_section=max(vector(section_min:section_max,1));
if max_section > treshold_level_L
        if  section_min+modifier*resolution <=size(vector,1)        
        section_max=section_min+modifier*resolution;
    else
        section_max=size(vector,1);
        end
        if  section_min-modifier*resolution >=1
        section_min=section_min-modifier*resolution;    
        else
            section_min=1;
        end        
        end  
local_max=max(vector(section_min:section_max,1));
    if local_max > treshold_level_H
    local_max_position=find(vector(section_min:section_max,1)==local_max);
     positions(index)=local_max_position+section_min-1;
        index=index+1;    
    end
    section_min=section_max+1;
end
end