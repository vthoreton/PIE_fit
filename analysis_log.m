function analysis_log(content,path)
% This function saves tables  Variables must be declared as global in the
% calling program.
%
% It is called with a line such as savloarch_log save_load_parameters(content,path,'save')
% content must be a 2 columns cell array defined as:
          
% the 1st column being the logged names of the parameters, the 2nd column
% being the corresponding variable.
% The path is a string indicating the location and name of the lof file such as: C:\Users\User\Desktop\Dataset\experimental_conditions_log.csv 

eval(sprintf('global %s',strjoin(content(:,2))))

export_fit_string='export_fit=cat(2,';
        export_fit_columns_string='export_fit_columns={';
        for i=1:size(content,1)
            if exist(string(content(i,2)))==1 & eval(sprintf('isempty(%s)==0',string(content(i,2)))) % & eval(sprintf('size(%s,1)==size(all_temperatures,1)',string(content(i,1))))
                export_fit_string=strcat(export_fit_string,string(content(i,2)),',');
                export_fit_columns_string=strcat(export_fit_columns_string,'''',string(content(i,1)),'''',',');
            end
        end
        export_fit_string=char(export_fit_string);
        export_fit_columns_string=char( export_fit_columns_string);
        export_fit_string=export_fit_string(1:end-1);
        export_fit_columns_string=export_fit_columns_string(1:end-1);
        export_fit_string=strcat(export_fit_string,');');
        export_fit_columns_string=strcat(export_fit_columns_string,'};');
        eval(export_fit_string);
        eval(export_fit_columns_string);
        export_fit=array2table(export_fit);
        export_fit.Properties.VariableNames = export_fit_columns;
        writetable(export_fit,path,'Delimiter',',');
end