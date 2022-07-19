function savloarch_log(content,path,sl) %function [conditions]= save_load_parameters(content,path,sl)
% This function can save,load variables (string and doubles) used as
% experimental parameters OR archive variables. Variables must be declared as global in the
% calling program.
%
% It is called with a line such as savloarch_log save_load_parameters(content,path,'save')
% content must be a 2 columns cell array defined as:
%  content={'sample' 'raw_data_filename';
%         'data_source' 'data_source';
%         'time_unit' 'time_unit';
%         'm'  'm_mg' ;
%         'M' 'M' ;
%         'rho' 'rho';};
% the 1st column being the logged names of the parameters, the 2nd column
% being the corresponding varialbe.
% The path is a string indicating the location and name of the lof file such as: C:\Users\User\Desktop\Dataset\experimental_conditions_log.csv 
% sl is 'save' or 'load'
%
% NB: This script is retro-compatible with former log files.

eval(sprintf('global %s',strjoin(content(:,2))))

%% LOADING
if strcmp(sl,'load')==1
parameters_to_load=content;
loaded_parameters=readcell(path);%preload to get the size
if size(loaded_parameters,1)==2
%% old format in lines!
loaded_parameters=readtable(path);
loaded_names=loaded_parameters.Properties.VariableNames'; %transposed
loaded_values=table2cell(loaded_parameters)';%transposed
elseif size(loaded_parameters,2)==2
%% new format in columns
loaded_parameters=readcell(path);
loaded_names=loaded_parameters(:,1);
loaded_values=loaded_parameters(:,2);
end
    for i = 1:size(parameters_to_load,1)
        position=find(double(ismember(loaded_names(:,1),parameters_to_load(i,1))));
        if isempty(position)==0
            if isnan(str2double(string(loaded_values(position,1))))==1
                eval(sprintf('%s=''%s'';',char(parameters_to_load(i,2)),char(loaded_values(position))));
            else
                eval(sprintf('%s=%g;',char(parameters_to_load(i,2)),str2double(string(loaded_values(position)))));
            end
        elseif isempty(position)==1
            eval(sprintf('%s=''NA'';',char(parameters_to_load(i,2))));
        end
    end
  
%% SAVING OR ARCHIVING 
elseif strcmp(sl,'save')==1 | strcmp(sl,'archive')==1
parameters_to_save=content;
export_fit_string='parameters_values={';
        export_names_string='parameters_names={';
        for i=1:size(parameters_to_save,1)
            if exist(string(parameters_to_save(i,2)))==1 & eval(sprintf('isempty(%s)==0',string(parameters_to_save(i,2))))
                export_fit_string=strcat(export_fit_string,string(parameters_to_save(i,2)),',');
                export_names_string=strcat(export_names_string,' ''',string(parameters_to_save(i,1)),'''');
            else
                export_fit_string=strcat(export_fit_string,' ''NA'' ,');
                export_names_string=strcat(export_names_string,' ''',string(parameters_to_save(i,1)),'''');
            end
        end
        export_fit_string=char(export_fit_string);
        export_names_string=char(export_names_string);
        export_fit_string=export_fit_string(1:end-1);
        export_names_string=export_names_string(1:end-1);
        export_fit_string=strcat(export_fit_string,'}'';');
        export_names_string=strcat(export_names_string,'''}'';');
        % evaluate the generated code:
        eval(export_fit_string);
        eval(export_names_string);
        
    if strcmp(sl,'save')==1
    parameters_log = [parameters_names,parameters_values];
    parameters_log=string(parameters_log);   
    writematrix(parameters_log,path);
    
    elseif strcmp(sl,'archive')==1
        if exist(path,'file') == 2
            archive_file=readcell(path);
            if size(parameters_names(:,1),1)==size(archive_file(:,1),1) & prod(cellfun(@strcmp, parameters_names(:,1), archive_file(:,1)))==1 %checks that it has the same size AND that contains the same parameters                
                archive_file=[archive_file,parameters_values];  % append a new column
            else
                [filepath,filename,ext] = fileparts(path);
                                new_filename=[filename,'_new',ext];
                new_path = fullfile(filepath,new_filename);
                 warning(sprintf('The dimension or content of the new log line does not match the original log file. A new log file has been generated: ''%s''.',new_filename));
                savloarch_log(parameters_to_save,new_path,'archive');                       
            end            
        else
    archive_file = [parameters_names,parameters_values]; % create a new log file
                end
        % write the file:
     archive_file=string(archive_file);   
    writematrix(archive_file,path); 
        
    end
    
end
end