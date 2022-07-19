    % PIE_fit v1.2 beta 19/04/2021
% by  Vincent Thoréton (email address for questions, comments and requests: vincent.thoreton.smn.uio.no (UiO) or  work.vincent@thoreton.net (permanent) )
% This script opens raw data from a folder named after the measurement temperature (ascii files exported by Quadera), plots and processes the data in a semi-automatic way.
% It has been written for a research and educational purposes. The code
% does is explained to a great extent. You may copy it or take inspiration
% from it to make your own routines. Suggestions are very welcome!
% Acknoledgements are also welcome if you are using this script for
% publication :-)


%% TO IMPROVE:
% reconstruct_data_file: deal with incomplete lines (requires Tord's data)
% summarising figure :done

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   List of associated scripts: %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   -reconstruct_data_file.m (homogeneise the number of columns of the data file)
%
%   -peaks_finder.m & peak_finder1.m (finding peaks) (v2 at beta stage v1 soon deprecated)
%
%   -peaks_integrator.m (peaks integration)
%
%   -fitPIE.m (fitting of p)
%
%   -savloarch_log.m (to load, save and archive log files)
%
%   -analysis_log.m (to export vectors in columns)
%
%   -rgb.m (third party function)

% clear the workspace, figures etc:
clear all;
clearvars;
close all force;
% declare the global variables that may be shared with subroutines and
% functions:    (NB: all variables that will be exported using the
% exporting functions MUST be declared as global)
global  integrals_log averaged_integrals stdev_log error_log fitting_data experimental_conditions_log_path masses_index timestamp sample_name data_source time_unit EA_eV_R0 EA_eV_Rads EA_eV_Rinc EA_R0 EA_Rads EA_Rinc stdev_R0 stdev_Rads stdev_Rinc M m_mg rho SSA pO2 MFC_T bed_length pulse_volume flow_rate followed_masses ref_mass f32_RT f34_RT f36_RT labelled_ref use_experimental_32 ID default_back_step default_integ_length background_type RB_type_override  LB_drv_trshld RB_drv_trshld default_LB_drv_trshld default_RB_drv_trshld all_temperatures all_TTT all_f18 all_frac32 all_frac34 all_frac36 stdev_34 stdev_36 err_34 err_36 all_k all_R0 temperatures TTT f18 frac32 frac34 frac36 k  R0 Rads Rinc p_calc Arrhenius_fit temper TTTT f18c frac32c frac34c frac36c R0c Radsc Rincc pc nu



%say hi
disp('Hello World!');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% ********************************************************************* %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% ***************** %% Choice of parameters and options %% ************ %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% ********************************************************************* %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%


%% setting experimental parameters:
%data_source='Quadera';%Quadera (Oslo) or Quadstar (Twente)
override_xp_cond='no'; % override saved data (enter it again) or not 

%% parameters of the peak finder:
% Optimised for data from Oslo (old values):
% % resolution=100;
% % modifier=4;
% % low_treshold=0.2;
% % high_treshold=0.5;
finder_version='v2';% v1 or v2; v1 is based on mean value of the whole vector, v2 is based on the local minimum. v2 can deal way better with background fluctuation

resolution=30;
modifier=3;
low_treshold=1;
high_treshold=2;

% Optimised for data from Oslo (old values):
% % resolution=50;
% % modifier=3;
% % low_treshold=0.5;
% % high_treshold=1.5;

% Optimised for data from Oslo with small peaks:
% % resolution=50;
% % modifier=3;
% % low_treshold=0.1;
% % high_treshold=0.5;

% Optimised for data from Oslo (old values):
% % resolution=100;
% % modifier=4;
% % low_treshold=0.2;
% % high_treshold=0.5;

% Optimised for data from Twente:
% % resolution=30;
% % modifier=3;
% % low_treshold=1;
% % high_treshold=2;


%% parameters for the integration:
default_back_step=50;
default_integ_length=100;
background_type='constant'; % constant, linear or constant_R
RB_type_override='no';% no, mid, dmid or auto-shrinked
default_LB_drv_trshld=1e-5;% default derivative treshold for the left bound
default_RB_drv_trshld=1e-5;% default derivative treshold for the right bound
% integration parameters are specific to each mass:
masses_settings=[2  4  14 16 17 18 19 20 21 22 28     30 32     34      36      40      44   46   48; %mass
                 0  0  0   0  0  0  0  0  0  0  0      0  0     0        0       0       0    0    0; %background
                 1  1  1   1  1  1  1  1  1  1  1      1  1     1        1       1.00	1.12 1.12 1.12; %yield N2: 1.3315 yield Ar: 1.28 ; set to 1 for simplicity
                 1  1  1   1  1  1  1  1  1  1 -1      1 -1     1        1       1       1    1    1; % positive (1) or negative (-1) peak
                 0  0  0   0  1  0  0  0  0  0  1      1  1     1        1       1       1    1    1;% precise maxima 1 for on, 0 for off
                -1 -1 -1  -1 -1 -1 -1 -1 -1 -1  1      1 -1    -1       -1       1       1    1    1;%LB (lower bound) location ; can be set to any number but 1,0 and -1 are reserved values: 1 for auto  (-1 is set to the same as the reference peak) 0 is simply set stupidly
                -1 -1 -1  -1 -1 -1 -1 -1 -1 -1  1      1 -1    -1       -1       1       1    1    1; % RB (right bound) location ; can be set to any number but 1,0 and -1 are reserved values: 1 for auto (-1 is set to the same as the reference peak) 0 is simply set stupidly
                -1 -1 -1  -1 1e-2 -1 -1 -1 -1 -1 -1     -1 -1     1       -1    1e-2      -1   -1   -1; % LB derivative treshold ; can be set to any number but -1 is reserved (set to the default value)
                -1 -1 -1  -1 -1 -1 -1 -1 -1 -1 -1     -1 -1     1     	-1    1e-5      -1   -1   -1]; % RB derivative treshold (absolute value) ; can be set to any number but -1 is reserved (set to the default value)
            
  masses_colors    ={2          4          14          16          17          18          19          20          21          22          28          30          32          34         36           40          44          46          48          ; % mass
                     'Silver'   'Maroon'   'Teal'      'Navy'      'RoyalBlue' 'DarkRed'   'Tomato'    'Gold'      'Gray'      'Black'     'Lime'     'Indigo'     'Blue'    'Purple'    'Red'        'Green'      'Cyan'      'Fuschia'   'Salmon'    };% color in plots
           
            
            
            
 % for futur feveloppement (not in use)):           
 masses_settings_v2={4          14          16          17          18          19          20          21          22          28          30          32          34         36           40          44          46          48          ; % mass
                     0          0           0           0           0           0           0           0           0           0           0           0           0          0            0           0           0           0           ; % MS background
                     1          1           1           1           1           1           1           1           1           1.3315      1           1           1          1            1.28        1.12        1.12        1.12        ; % MS yield
                     1          1           1           1          -1           1          -1           1           1           1           1           1           1          1            1           1           1           1           ; % positive (1) or negative (-1) peak
                     1          1           1           1          -1           1           0           0           1           1           1           1           1          1            1           1           1           1           ; % precise maxima 1 for on, 0 for off
                     1          1           1           1           1           1           1           1           1           1           1           1           1          1            1           1           1           1           ; % LB (lower bound) location ; can be set to any number; 1,0 and -1 are reserved values: 1 for auto, -1 sets it to the same as the reference peak, 0 is set to the default value 
                     1          1           1           1           1           'mid'       1          -1          -1           1           1           1           1          1            1           1           1           1           ; % RB (upper bound) location ; can be set to any number or a string; 1,0 and -1 are reserved values: 1 for auto, -1 sets it to the same as the reference peak, 0 is set to the default value. Alternatively, the string can be set to 'mid', 'dmid' 'auto', '???' or '???')
                     1          1           1           1           1           1           1           1           1           1           1           1           1          1            1           1           1           1           ; % background type (0 for constant (left bound), l for linear (from left to right bound), 2 for constant (right bound)
                     'Maroon'   'Teal'      'Navy'      'RoyalBlue' 'DarkRed'   'Tomato'    'Gold'      'Gray'      'Black'     'Lime'     'Indigo'     'Blue'    'Purple'    'Red'        'Green'      'Cyan'      'Fuschia'   'Salmon'    };% color in plots
% %             
% "Optimised" for data from Twente:           
% % masses_settings=[4 14 16 17 18 19 20 21 22 28     30 32 34 36 40   44   46   48; %mass
% %                  0 0  0  0  0  0  0  0  0  0      0  0  0  0  0    0    0    0; %background
% %                  1 1  1  1  1  1  1  1  1  1.3315 1  1  1  1  1.28 1.12 1.12 1.12; %yield
% %                  1 1  1  1  1  1  1  1  1  1      1 -1  1  1  1    1    1    1; % positive (1) or negative (-1) peak
% %                  1 1  1  1  1  1  1  1  1  1      1  0  0  1  1    1    1    1;% precise maxima 1 for on, 0 for off
% %                  1 1  1  1  1  1  1  1  1  6      1 -1 -1 -1  1    1    1    1;%LB (lower bound) location ; can be set to any number but 1,0 and -1 are reserved values: 1 for auto  (-1 is set to the same as the reference peak) 0 is simply set stupidly
% %                  1 1  1  1  1  1  1  1  1  22     1 -1 -1 -1  1    1    1    1]; % RB (right bound) location ; can be set to any number but 1,0 and -1 are reserved values: 1 for auto (-1 is set to the same as the reference peak) 0 is simply set stupidly
             
             
             
%% plotting parameters:
plot_found_peaks='no'; % useful to check visually at the 1st run that all the peaks were found
plot_individual_peaks='no'; % no by default
peaks_YScale='log'; % log or lin
plot_ref_34_36='no';
plot_integrals='yes';
plot_all_integrals='no' ;
integrals_YScale='lin'; % log or lin
left_gap_peaks_display=10;
right_gap_peaks_display=30;

%% Fitting parameters:
argon_36_correction='yes';
argon36_corrective_factor=0.3336/99.6035;
options_fit=optimset('Display','notify','TolX',1e-11,'TolFun',1e-11,'MaxIter',1000);% default optimisation of the fminsearch function used for fitting
force_reference='no';%GPA or thermo or no (log)

%% Export parameters:
logging_folder='fitting'; % set it back to 'processed_data' if you want to keep compatibility with v1 beta
fig_res=300; % figures resolution
figure_format='png'; % png or jpg
export_fig='no'; % additionally exports matlab figures (fig format)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a few constants %
%%%%%%%%%%%%%%%%%%%%%%%%%%
R=8.314462;% ideal gas constant (J/mol/K)
FarC=96485.33212; % Faraday constant (C/mol)

%%%%%%%%%%%%%%%%%%%%
%% Importing data %%
%%%%%%%%%%%%%%%%%%%%
timestamp=datestr(now,'mmmm-dd-yyyy HHMMSS'); %timestamp used for naming each log folder
%choose the raw data file (directly exported from the MS software)
[raw_data_filename,raw_data_root_folder] = uigetfile('*.asc','Choose any temperature file from the working folder');
saving_path=fullfile(raw_data_root_folder,logging_folder);
saving_directory=fullfile(saving_path,sprintf('%s',timestamp));%set the path for the saving directory
script_folder = pwd;
eval(sprintf('cd (''%s'');',raw_data_root_folder));
[~,sample_name,~]=fileparts(pwd);
temperatures_index = dir('*.asc');
temperatures_index=struct2cell(temperatures_index);
for temp_index=1:1:size(temperatures_index,2)
    file=char(temperatures_index(1,temp_index));
    temperature=char(temperatures_index(1,temp_index));
   [filepath,temperature,ext] = fileparts(temperature);    
    temperature=str2double(temperature);
    date=char(temperatures_index(3,temp_index));
    transitional_temperatures_index{1,temp_index}=file;
    transitional_temperatures_index{2,temp_index}=temperature;
    transitional_temperatures_index{3,temp_index}=date;
end
temperatures_index=transitional_temperatures_index;
temperatures_index=sortrows(temperatures_index',2)';% sorts the temperatures in ascending order
temperatures=str2double(string(temperatures_index(2,:)))';
eval(sprintf('cd ''%s'';',script_folder))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load or enter the experimental conditions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  old_parameters_list={'sample' 'raw_data_filename'; % for compatibility with former log version
        'data_source' 'data_source';
        'time_unit' 'time_unit';
        'm'  'm_mg' ;
        'M' 'M' ;
        'rho' 'rho';
        'SSA' 'SSA' ;
        'pO2' 'pO2';
        'MFC_T'  'MFC_T';
        'bed_length' 'bed_length';
        'flow_rate' 'flow_rate';
        'pulse_vol' 'pulse_volume';
        'followed_masses' 'followed_masses';
        'ref_mass' 'ref_mass';
        'f34' 'f34_RT';
        'f36' 'f36_RT';
        'reference' 'labelled_ref';
        'xp_32' 'use_experimental_32';
        'inner_diam' 'ID'};
    
    new_parameters_list={'sample'  'sample_name'; % current version (more readable with 1 parameters per row)
             'data_source' 'data_source';
             'time_unit' 'time_unit';
             'm [mg]' 'm_mg';
             'M [g/mol]' 'M';
             'nu' 'nu';
             'rho [g/cc]' 'rho';
             'SSA [m2/g]' 'SSA';
             'pO2 [bar]' 'pO2';
             'MFC_T [degC]' 'MFC_T';
             'bed_length [mm]' 'bed_length';
            'flow_rate [mL/min]' 'flow_rate';
             'pulse_vol [uL]' 'pulse_volume';
            'followed_masses' 'followed_masses';
            'ref_mass' 'ref_mass';
            'f34' 'f34_RT';
            'f36' 'f36_RT';
              'reference' 'labelled_ref';
             'xp_32' 'use_experimental_32' ;
             'inner_diam. [mm]' 'ID' };
         
         

         
         


experimental_conditions_log_path = fullfile(saving_path,'experimental_conditions_log.csv');
if exist(experimental_conditions_log_path,'file') == 2  %load the parameters

    if strcmp(override_xp_cond,'yes')==0
    experimental_conditions=readtable(experimental_conditions_log_path);
    experimental_conditions_content=experimental_conditions.Properties.VariableNames;%get the names of variables from the log
    experimental_conditions=table2cell(experimental_conditions);%converts the table to array before extracting the values
    
    if size(experimental_conditions,1)==1
%% old format in 2 lines!
    parameters_to_load=old_parameters_list;    
    elseif size(experimental_conditions,2)==2    
        parameters_to_load=new_parameters_list;
    end
    parameters_to_save=new_parameters_list;    
 savloarch_log(parameters_to_load,experimental_conditions_log_path,'load');    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% OR input parameters in a dialog %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %defaults parameters
    %% sample
    data_source='Quadera';
    time_unit='ms';
    m_mg=24.2;%sample mass (mg)
    M=136.74; %g/mol
    nu=3;
    rho=4.01;%g/cm3
    bed_length=5;% mm
    SSA=0.28;%m2/g
    %% gas
    flow_rate=50;%ml/min
    pO2=0.02;%atm
    pulse_volume=500;%uL
    followed_masses='40 34 36';
    ref_mass=40; %change to 4 if helium is used etc
    MFC_T=27;% temperature of the MFC (degC)
    % isotopologue fractions at RT:
    f36_RT=0.935;
    f34_RT=0.0582;
    labelled_ref='GPA'; % thermo or GPA
    use_experimental_32='no';
    %% reactor
    ID=2;%inner diameter (mm)
    
    % input of the experimental parameters
    default_input =string({data_source time_unit m_mg M nu rho bed_length SSA flow_rate pO2 pulse_volume followed_masses ref_mass MFC_T f34_RT f36_RT labelled_ref use_experimental_32 ID});
    experimental_parameters = inputdlg({'Data source [Quadera (e.g. from Oslo) or Quadstar(e.g. from Twente)]','Time unit [ms or s]','m (mg)','M (g/mol)','Oxygen stiochiometry \\nu','rho (g/cm³)','bed length (mm)','specific surface area (m²/g)','flow rate (mL/min)','pO2 (bar)','pulse volume (uL)','followed masses','reference mass','MFC temperature (°C)','f34 (GPA)','f36 (GPA)','ref. calibration [GPA or thermo]','use exp. 32 [uses mass 32 if it was recorded]','reactor inner diameter (mm)'},...
        'Input of the experimental parameters',...
        [1 35],default_input);
    if size(experimental_parameters)~=0
        data_source=experimental_parameters(1);
        time_unit=experimental_parameters(2);
        m_mg=str2double(experimental_parameters(3));
        M=str2double(experimental_parameters(4));
        nu=str2double(experimental_parameters(5));
        rho=str2double(experimental_parameters(6));
        bed_length=str2double(experimental_parameters(7));
        SSA=str2double(experimental_parameters(8));
        flow_rate=str2double(experimental_parameters(9));
        pO2=str2double(experimental_parameters(10));
        pulse_volume=str2double(experimental_parameters(11));
        followed_masses=experimental_parameters(12);
        ref_mass=str2double(experimental_parameters(13));
        MFC_T=str2double(experimental_parameters(14));
        f34_RT=str2double(experimental_parameters(15));
        f36_RT=str2double(experimental_parameters(16));
        labelled_ref=string(experimental_parameters(17));
        use_experimental_32=experimental_parameters(18);
        ID=str2double(experimental_parameters(19));        
    elseif size(experimental_parameters)==0  % if cancel is pressed
        choice = questdlg('Continue with default values?',...
            'Warning', ...
            'yes','no','no');
        switch choice
            case 'no'
                return
            case 'yes'
        end
    end    
    %save  the experimental conditions log if not there yet:
    mkdir(fullfile(raw_data_root_folder,logging_folder));%create the logging folder
     parameters_to_save=new_parameters_list;
end % end of input of parameters (saved of dialog)
   
% for compatibility with previous versions (prior the introduction of the data source):
if strcmp(data_source,'NA')==1
    data_source='Quadera';
end
if strcmp(time_unit,'NA')==1
    time_unit='ms';
end

if strcmp(nu,'NA')==1
  nu = inputdlg({'Oxygen stoichiometry (nu)'},...
        'Input of the experimental parameters');
    nu=str2double(nu);
end

    savloarch_log(parameters_to_save,experimental_conditions_log_path,'save');%save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Preliminary calculations  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=m_mg/1E3;
f32_RT=1-f36_RT-f34_RT;
f18_GPA=f36_RT+f34_RT/2;% for comparison with the thermodynamic equilibrium approximation
molar_flow_rate=(2*pO2*10^5/(R*(273.15+MFC_T)))*flow_rate*10^-6/60;
c_0=3*rho*10^6/M;
S=m*SSA;
residence_time=pulse_volume/1000/flow_rate*60; % s
pulse_length=pulse_volume/(pi*(ID/2)^2); % mm
bed_porosity=m/(bed_length/10*pi()*(ID/20)^2)/rho;
bed_density=1-bed_porosity;
n_O_pulse=2*pO2*10^5*pulse_volume*10^-9/R/(MFC_T+273.15);
n_O_bed=nu*m/M;

followed_masses=split(followed_masses)';
possible_masses_num=str2double(followed_masses);
fig_res=sprintf('-r%d',fig_res);

if strcmp(force_reference,'GPA')==1
    labelled_ref='GPA';
elseif strcmp(force_reference,'thermo')==1
    labelled_ref='thermo';
elseif strcmp(force_reference,'no')==1
end
if strcmp(labelled_ref,'thermo')==1
    content_annotation_ref={'Labelled gas reference:','thermodynamic approximation'};
elseif strcmp(labelled_ref,'GPA')==1
    content_annotation_ref={'Labelled gas reference:','GPA measurement'};
end
if strcmp(use_experimental_32,'yes')==1
content_annotation_ref={'Labelled gas reference:','internal reference (mass 32)'};
end

message=msgbox('Processing data, please wait...','Notification');


%% SCAN THE FILES BY INCREASING TEMPERATURE
for temp_index=1:1:size(temperatures_index,2)
    temperature=str2double(string(temperatures_index(2,temp_index)));
    path=strcat(raw_data_root_folder,char( temperatures_index(1,temp_index)));
    file=reconstruct_data_file(sprintf('%s',path));
    
     if strcmp(data_source,'Quadera')==1
    line_num_of_masses=find(strcmp(file(:,1),'Time'))-1;    
    % locate the masses of interest and create a specific array for each of them
    masses_index=[];%list of the recorded mass and corresponding column
            %% SCAN THE CURRENT FILE (TEMPERATURE) FOR AVAILABLES MASSES
    for i = 1:size(followed_masses,2)
        col=find(strcmp(file(line_num_of_masses,:),followed_masses{1,i}));
        if strcmp(file(line_num_of_masses,col),followed_masses{1,i}) == 1
            masses_index=[masses_index [possible_masses_num(1,i);col]];
        end
    end
    %% SCAN THE CURRENT TEMPERATURE BY MASS AND CREATE VECTORS  BEFORE INTEGRATION
    for mm = 1:size(masses_index,2)%creates a 2 columns array (time and intensity) for each mass (time and intensity)
        %eval(sprintf('clear raw_%g',masses_index(1,m)));
        eval(sprintf('raw_%g=[];',masses_index(1,mm)));
        if strcmp(time_unit,'ms')==1        
        eval(sprintf('raw_%g(:,1)=str2double(file(line_num_of_masses+3:end-1,masses_index(2,mm)-1))/1000;',masses_index(1,mm)));
        elseif strcmp(time_unit,'s')==1
        eval(sprintf('raw_%g(:,1)=str2double(file(line_num_of_masses+3:end-1,masses_index(2,mm)-1));',masses_index(1,mm)));
        end
        eval(sprintf('raw_%g(:,2)=str2double(file(line_num_of_masses+3:end-1,masses_index(2,mm)));',masses_index(1,mm)));
        col=find(masses_settings(1,:)==masses_index(1,mm));
        eval(sprintf('clear corr_%g',masses_index(1,mm)));
        eval(sprintf('corr_%g(:,1)=raw_%g(:,1);',masses_index(1,mm),masses_index(1,mm)));
        eval(sprintf('corr_%g(:,2)=(raw_%g(:,2)-masses_settings(2,col))/masses_settings(3,col);',masses_index(1,mm),masses_index(1,mm)));
    end
    
  elseif strcmp(data_source,'Quadstar')==1 
  [experiment_name,experiment_date,experiment_time,raw_data_headers,masses_index,raw_data] = import_QS32(file);
   for mm = 1:size(masses_index,2)%creates a 2 columns array (time and intensity) for each mass (time and intensity)
        eval(sprintf('raw_%g=[];',masses_index(1,mm)));
  eval(sprintf('raw_%g(:,1)=raw_data(:,1);',masses_index(1,mm)));
        eval(sprintf('raw_%g(:,2)=raw_data(:,%g);',masses_index(1,mm),masses_index(2,mm)));
        col=find(masses_settings(1,:)==masses_index(1,mm));
                eval(sprintf('corr_%g=[];',masses_index(1,mm)));
 eval(sprintf('corr_%g(:,1)=raw_%g(:,1);',masses_index(1,mm),masses_index(1,mm)));
        eval(sprintf('corr_%g(:,2)=(raw_%g(:,2)-masses_settings(2,col))/masses_settings(3,col);',masses_index(1,mm),masses_index(1,mm)));
   end
   end   
        eval(sprintf('reference_vector=corr_%g;',ref_mass)) 
        
        %%%%%%%%%%%%%%%%%%%%%% Debug same time scale ???????Qu'est ce que
        %%%%%%%%%%%%%%%%%%%%%% c'est que ca ???????????       
% % %             for m = 1:size(masses_index,2)
% % %          col=find(masses_settings(1,:)==masses_index(1,m));
% % %         eval(sprintf('corr_%g(:,1)=reference_vector(:,1);',masses_index(1,m)));
% % %     end        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    %% SEARCHING FOR PEAKS
    if strcmp(finder_version,'v1')==1
    peaks_coordinate=peaks_finder1(reference_vector(:,2),resolution,modifier,low_treshold,high_treshold);%',ref_mass))%using reference mass for finding peaks ; the peak integrator will determine the position of individual peaks more accurately
    elseif strcmp(finder_version,'v2')==1
    peaks_coordinate=peaks_finder(reference_vector(:,2),resolution,modifier,low_treshold,high_treshold);%',ref_mass))%using reference mass for finding peaks ; the peak integrator will determine the position of individual peaks more accurately
    end
    
    %% DISPLAYING FOUND PEAKS     
    if strcmp(plot_found_peaks,'yes')==1   
        content_annotation={sprintf('Peak(s) found: %g',size(peaks_coordinate,2))};
        eval(sprintf('plot_corr_log_%gC=figure(''Name'',''Pulses at %g %cC'');',temperature,temperature,char(176)))
        for var=1:size(masses_index,2)            
            col=find(masses_settings(1,:)==masses_index(1,var));
            color=rgb(string(masses_colors(2,col)));
            r=color(1);
            g=color(2);
            b=color(3);            
        eval(sprintf('plot(corr_%g(:,1),corr_%g(:,2),''.'',''color'',[%g %g %g],''Linewidth'',2)',masses_index(1,var),masses_index(1,var),r,g,b))      
        hold on;
        end 
        %eval(sprintf('set(gca,''XScale'',''lin'',''YScale'',''%s'')',peaks_YScale))%
        set(gca, 'YScale', 'lin')
        legend_string='legend(''i_{';
        for var=1:size(masses_index,2)-1
        legend_string=strcat(legend_string,string(masses_index(1,var)),'}'',''i_{');
        end
        legend_string=strcat(legend_string,string(masses_index(1,end)),'}'',''Location'',''northeast'');');
        eval(legend_string);             
        xlabel('time (s)')
        ylabel('Current (A)')
        axis([reference_vector(peaks_coordinate(1),1)-left_gap_peaks_display reference_vector(peaks_coordinate(end),1)+right_gap_peaks_display 0 max(reference_vector(:,2))*11/8])
        eval(sprintf('title(''Peaks at %g %cC'')',temperature,char(176)))
        annotation('textbox',[0.4,0.9,0,0],'String',content_annotation,'FitBoxToText','on')
        
          %indicates peaks with arrows:
                        for n = 1:size(peaks_coordinate,2)
            eval(sprintf('arrow%g  = annotation(''arrow'');',n))
            eval(sprintf('arrow%g.Parent = gca;',n) )
            eval(sprintf('arrow%g.Position= [reference_vector(peaks_coordinate(n),1), max(reference_vector(:,2))+0.01, 0, -0.0001];',n));
            end    
        
        
              
     if strcmp(plot_individual_peaks,'yes')==1         
                     for peak_i = 1:size(peaks_coordinate,2)                
                eval(sprintf('peak%g_%gC=figure(''Name'',''Pulse %g at %g %cC'');',peak_i,temperature,peak_i,temperature,char(176)))
        for var=1:size(masses_index,2)            
            col=find(masses_settings(1,:)==masses_index(1,var));
            color=rgb(string(masses_colors(2,col)));
            r=color(1);
            g=color(2);
            b=color(3);            
        eval(sprintf('plot(corr_%g(:,1),corr_%g(:,2),''.'',''color'',[%g %g %g],''Linewidth'',2)',masses_index(1,var),masses_index(1,var),r,g,b))      
        hold on;
        end 
        eval(sprintf('set(gca,''XScale'',''lin'',''YScale'',''%s'')',peaks_YScale))%set(gca, 'YScale', 'log')
        legend_string='legend(''i_{';
        for var=1:size(masses_index,2)-1
        legend_string=strcat(legend_string,string(masses_index(1,var)),'}'',''i_{');
        end
        legend_string=strcat(legend_string,string(masses_index(1,end)),'}'',''Location'',''northeast'');');
        eval(legend_string);             
        xlabel('time (s)')
        ylabel('Current (A)')
        %axis([reference_vector(peaks_coordinate(1),1)-left_gap_peaks_display reference_vector(peaks_coordinate(end),1)+right_gap_peaks_display 0 max(reference_vector(:,2))*11/8])
        axis([reference_vector(peaks_coordinate(peak_i),1)-left_gap_peaks_display reference_vector(peaks_coordinate(peak_i),1)+right_gap_peaks_display  0 max(reference_vector(:,2))*11/8])%',masses_index(1,var),masses_index(1,var)))
        eval(sprintf('title(''Peak %g at %g %cC'')',peak_i,temperature,char(176)))
        eval(sprintf('set(gca,''XScale'',''lin'',''YScale'',''%s'')',peaks_YScale))
        end
     end
    
        
     
       
        if strcmp(plot_ref_34_36,'yes')==1
        content_annotation={sprintf('Peak(s) found: %g',size(peaks_coordinate,2))};
        eval(sprintf('plot_corr_lin_%gC=figure(''Name'',''Pulses at %g %cC'');',temperature,temperature,char(176)))
        hold on;
        plot(reference_vector(:,1),reference_vector(:,2),'k.','Linewidth',2)
        plot(corr_34(:,1),corr_34(:,2),'.','color',[1.0  0.0  1.0],'Linewidth',2)
        plot(corr_36(:,1),corr_36(:,2),'r.','Linewidth',2)   
        eval(sprintf('legend(''i_{34}'',''i_{36}'',''i_{%g}'',''Location'',''east'')',ref_mass))
        xlabel('time (s)')
        ylabel('Current (A)')
        axis([reference_vector(peaks_coordinate(1),1)-left_gap_peaks_display reference_vector(peaks_coordinate(end),1)+right_gap_peaks_display 0 max(reference_vector(:,2))*11/8])
        eval(sprintf('title(''Peaks at %g %cC'')',temperature,char(176)))
        annotation('textbox',[0.4,0.9,0,0],'String',content_annotation,'FitBoxToText','on')
                %% OR:
        % disp_num_peaks = annotation('textbox','String',content_annotation,'FitBoxToText','on');
        % disp_num_peaks.Parent = gca;
        % disp_num_peaks.Position= [corr_40(peaks_coordinate(end),1), max(corr_40(:,2))+0.01, 20, 0];

                %indicates peaks with arrows:
                        for n = 1:size(peaks_coordinate,2)
            eval(sprintf('arrow%g  = annotation(''arrow'');',n))
            eval(sprintf('arrow%g.Parent = gca;',n) )
            eval(sprintf('arrow%g.Position= [reference_vector(peaks_coordinate(n),1), max(reference_vector(:,2))+0.01, 0, -0.0001];',n));
            end                
        end 
    end
    
    
    %% INTEGRATION OF PEAKS
    integrals_log_NB=[];%% NB stands for "New Block"
    
    
    % the integration parameters (back_Step and integ_lenght in particular) must be first
        % determined on the reference mass:
        eval(sprintf('corr_ref_mass=corr_%g;',ref_mass))
        ref_col=find(masses_settings(1,:)==ref_mass);
        ref_peak_symbol=masses_settings(4,ref_col);
        ref_precise_extrema=masses_settings(5,ref_col);
        ref_LB_type=masses_settings(6,ref_col);
        ref_RB_type=masses_settings(7,ref_col);
        ref_LB_drv_trshld=masses_settings(8,ref_col);
        ref_RB_drv_trshld=masses_settings(9,ref_col);
        
        if ref_precise_extrema==1
            ref_precise_extrema='yes';
        elseif ref_precise_extrema==0
            ref_precise_extrema='no';
        end
        
        if ref_LB_drv_trshld==-1
        ref_LB_drv_trshld=default_LB_drv_trshld;
        end
        if ref_RB_drv_trshld==-1
        ref_RB_drv_trshld=default_RB_drv_trshld;
        end
        
        if ref_LB_type==1
            ref_LB_type='deriv.';
            ref_back_step=default_back_step;
        elseif ref_LB_type==0
            ref_LB_type='set';
            ref_back_step=default_back_step;
            elseif ref_LB_type~=0 & ref_LB_type~=1
              default_back_step=ref_LB_type;  
              ref_LB_type='set';
        end
        
        if ref_RB_type==1
            ref_RB_type='deriv.';
            ref_integ_length=default_integ_length;
        elseif ref_RB_type==0
            ref_RB_type='set';
            ref_integ_length=default_integ_length;
            elseif strcmp(ref_RB_type,'mid')==1 | strcmp(ref_RB_type,'dmid')==1
            
             elseif ref_RB_type~=0 & ref_RB_type~=1
              default_integ_length=ref_RB_type; 
              ref_RB_type='set';
        end
        ref_background_type=background_type;
        ref_peaks_coordinate=peaks_coordinate;
        integration_method={default_back_step,default_integ_length,ref_precise_extrema,ref_LB_type,ref_RB_type,ref_background_type,ref_peak_symbol,ref_LB_drv_trshld,ref_RB_drv_trshld};
        [ref_integrals ref_peaks_coordinate ref_back_step ref_integ_length]=peaks_integrator(corr_ref_mass,ref_peaks_coordinate,integration_method,RB_type_override);

    
    
    
    
    
    
    for n = 1:size(masses_index,2)
        masss=masses_index(1,n);
        eval(sprintf('corr_current_mass=corr_%g;',masss))
        col=find(masses_settings(1,:)==masses_index(1,n));
        peak_symbol=masses_settings(4,col);
        precise_extrema=masses_settings(5,col);
        LB_type=masses_settings(6,col);
        RB_type=masses_settings(7,col);
        LB_drv_trshld=masses_settings(8,col);
        RB_drv_trshld=masses_settings(9,col);
        if precise_extrema==1
            precise_extrema='yes';
        elseif precise_extrema==0
            precise_extrema='no';
        end
        if LB_drv_trshld==-1
        LB_drv_trshld=default_LB_drv_trshld;
        end
        if RB_drv_trshld==-1
        RB_drv_trshld=default_RB_drv_trshld;
        end
        
                
        
        if LB_type==1
            LB_type='deriv.';
            back_step=default_back_step;
        elseif LB_type==0
            LB_type='set';
            back_step=default_back_step;
        elseif  LB_type==-1
            LB_type='set (ref)';
            back_step=ref_back_step;
        elseif LB_type~=-1 & LB_type~=0 & LB_type~=1
                  back_step=LB_type;   
        LB_type='set';
        end
        
 
        if RB_type==1
            RB_type='deriv.';
            integ_length=default_integ_length;
        elseif RB_type==0
            RB_type='set';
            integ_length=default_integ_length;
        elseif  RB_type==-1
            RB_type='set (ref)';
            integ_length=ref_integ_length;
        elseif RB_type~=-1 & RB_type~=0 & RB_type~=1
               integ_length=RB_type;    
            RB_type='set';
                
        end
        
                if strcmp(RB_type_override,'no')==1                                 
RB_type_disp=RB_type;
                else                    
        RB_type_disp=RB_type_override;        
                end
        
        integration_method={back_step,integ_length,precise_extrema,LB_type,RB_type,background_type,peak_symbol,LB_drv_trshld,RB_drv_trshld};
        [integrals peaks_coordinate back_step integ_length]=peaks_integrator(corr_current_mass,peaks_coordinate,integration_method,RB_type_override);
        integrals_log_NB=horzcat(integrals_log_NB,integrals');
       sampling_rate=1/((corr_current_mass(end,1)-corr_current_mass(1,1))/(size(corr_current_mass,1)-1));
        
      if strcmp(argon_36_correction,'yes')==1 & exist('corr36')& exist('corr40') & ref_mass==40
col40=find(masses_index(1,:)==40);
col36=find(masses_index(1,:)==36);
integrals_log(:,1+col36)=integrals_log(:,col36)-integrals_log(:,col40)*argon36_corrective_factor;
else
end
        
        if peak_symbol==1
            peak_symbol='positive';
        elseif peak_symbol==-1
            peak_symbol='negative';
        end
        
        
        %% PLOTTING THE INTEGRATION PROCESS FOR VISUAL VERIFICATION
        if strcmp(plot_integrals,'yes')==1
            if strcmp(plot_all_integrals,'yes')==0
                kk_max=1;
            else
                kk_max=size(peaks_coordinate,2);
            end
            for kk = 1:kk_max
                eval(sprintf('plot_corr_log_%g_%gC=figure(''Name'',''Pulse %g of mass %g at %g %cC'');',masss,temperature,kk,masss,temperature,char(176)))
                eval(sprintf('semilogy(corr_%g(:,1),corr_%g(:,2),''r.'',''Linewidth'',2)',masss,masss))
                hold on;
                xlabel('time (s)')
                ylabel('Current (A)')
                eval(sprintf('axis([corr_%g(peaks_coordinate(kk),1)-left_gap_peaks_display corr_%g(peaks_coordinate(kk),1)+right_gap_peaks_display  -inf inf])',masss,masss))
                left_base_coordinate=peaks_coordinate(kk)-back_step(kk); %% kk remove)
                right_base_coordinate=peaks_coordinate(kk)-back_step(kk)+integ_length(kk);%% kk remove)
                eval(sprintf('time_coordinates=corr_%g(left_base_coordinate:right_base_coordinate,1);',masss))
                eval(sprintf('semilogy(time_coordinates,corr_%g(left_base_coordinate:right_base_coordinate,2));',masss))
                
                eval(sprintf('integ_curve=[corr_%g(left_base_coordinate:right_base_coordinate,2)];',masss)) %%%%%%%%%! trick! ; background_curve(end)

                if strcmp(background_type,'linear')==1
                    
                    eval(sprintf('background_curve=linspace(corr_%g(left_base_coordinate,2),corr_%g(right_base_coordinate,2),size(integ_curve,1))'';',masss,masss))
                    %% if the boundary values are the same (rare but it happened!!) the following line just returns [] ... 
                    %eval(sprintf('background_curve=(corr_%g(left_base_coordinate,2):(corr_%g(right_base_coordinate,2)-corr_%g(left_base_coordinate,2))/integ_length(kk):corr_%g(right_base_coordinate,2))'';',masss,masss,masss,masss))
                                
                elseif strcmp(background_type,'constant')==1
                    eval(sprintf('background_curve=zeros(integ_length(kk)+1,1)+corr_%g(left_base_coordinate,2);',masss))
                elseif strcmp(background_type,'constant_R')==1
                    eval(sprintf('background_curve=zeros(integ_length(kk)+1,1)+corr_%g(right_base_coordinate,2);',masss))
                end
                
                content_annotation={sprintf('\\intidt = %g A.s',integrals(kk)),sprintf('Sampling rate = %0.1f Hz',sampling_rate),sprintf('Position: %g',peaks_coordinate(kk)),sprintf('Back-steps: %g',back_step(kk)),sprintf('Integration length: %g',integ_length(kk)),sprintf('Background: %s',background_type),sprintf('Precise extremum: %s',precise_extrema),sprintf('Lower bound: %s',LB_type),sprintf('Integ. length: %s',RB_type_disp),sprintf('Deriv. tresh. (left): %0.0e',LB_drv_trshld),sprintf('Deriv. tresh. (right): %0.0e',RB_drv_trshld)};%\\partial ,sprintf('Type of peak: %s',peak_symbol)
                annotation('textbox',[0.45 0.55 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
                patch([time_coordinates' fliplr(time_coordinates')],[background_curve' fliplr(integ_curve')],'r')
                eval(sprintf('set(gca,''XScale'',''lin'',''YScale'',''%s'')',integrals_YScale))
                eval(sprintf('legend(''i_{%g}'',''Location'',''east'')',masss))
                eval(sprintf('title(''Integration of peak %g of mass %g at %g %cC'')',kk,masss,temperature,char(176)))
                %  end
            end
        end
    end
    
    eval(sprintf('col=find(masses_index(1,:)==%g);',ref_mass))
    norm_integrals=integrals_log_NB(:,1:end)./integrals_log_NB(:,col);
    means_NB=mean(norm_integrals,1);
    if min(size(norm_integrals))==1
    stdev_NB=zeros(1,max(size(norm_integrals)));    
    else
    stdev_NB=std(norm_integrals,1);
    end    
    error_NB=stdev_NB./means_NB;
    averaged_integrals=vertcat(averaged_integrals,means_NB);
    stdev_log=vertcat(stdev_log,stdev_NB);
    error_log=vertcat(error_log,error_NB);
    
    integrals_log_NB=horzcat(zeros(size(peaks_coordinate,2),1)+temperature,integrals_log_NB);
    integrals_log=vertcat(integrals_log,integrals_log_NB);    
end


col=find(masses_index(1,:)==34);
stdev_34=stdev_log(:,col);
err_34=error_log(:,col);

col=find(masses_index(1,:)==36);
stdev_36=stdev_log(:,col);
err_36=error_log(:,col);

delete(message);

if strcmp(plot_integrals,'yes')==1 | strcmp(plot_found_peaks,'yes')
    choice = questdlg('Save the figures?', ...
        'Notification', ...
        'yes','no','yes');
    switch choice
        case 'no'
            close all;
        case 'yes'
                        message=msgbox('Saving figures, please wait...','Notification');
        if exist(saving_directory,'dir') == 7
        else
                mkdir(saving_directory);%create a directory for each fit
                end
           % saving_directory=strcat(saving_path,sprintf('%s',timestamp),'\');%set the path for the saving directory
            
            FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
            for iFig = 2:length(FigList)
                FigHandle = FigList(iFig);
                FigName   = get(FigHandle, 'Name');
                % savefig(FigHandle, fullfile(saving_directory, FigName, '.fig'));
                print('-dpng',fig_res,FigHandle,fullfile(saving_directory,sprintf('%s.png',FigName)));
            end
            delete(message);
            close all;
    end    
end

choice = questdlg('Continue with fitting?', ...
    'Notification', ...
    'yes','no','yes');
switch choice
    case 'no'
        return
    case 'yes'
end

%% Various data handling
col32=find(masses_index(1,:)==32);
col34=find(masses_index(1,:)==34);
col36=find(masses_index(1,:)==36);
room_temperature=min(temperatures);
rowRT=find(temperatures(:,1)==room_temperature);
A34=averaged_integrals(rowRT,col34);

% if strcmp(,'yes')==1
% A36=averaged_integrals(rowRT,col36);
% else
A36=averaged_integrals(rowRT,col36);
%end

if isempty(col32)==0
    A32_exp=averaged_integrals(rowRT,col32);
end

if strcmp(use_experimental_32,'yes')==1 & isempty(col32)==0
    AO2=-averaged_integrals(rowRT,col32);
    A32=1-A36-A34;
else
    A32=A36*(1/f36_RT-1)/(f32_RT+f34_RT)-A34-A36;
    AO2_GPA=A32+A34+A36;
    AO2_thermo=A34+A36+A34^2/(4*A36);
    
    if strcmp(labelled_ref,'GPA')==1
        AO2= AO2_GPA;
    elseif strcmp(labelled_ref,'thermo')==1
        AO2= AO2_thermo;
    end
end
% averaged_integrals2=averaged_integrals;
% averaged_integrals2(1:end,col34)=averaged_integrals2(1:end,col34)+averaged_integrals{1,col34}
%
% averaged_integrals

%for u=1:size(followed_masses,2)
%   eval(sprintf('frac%s=averaged_integrals(:,col%s)/AO2;;',string(followed_masses(1,u)),string(followed_masses(1,u))))
%eval(sprintf('all_frac%g=frac%g;',str2double(split(followed_masses(1,u))),str2double(split(followed_masses(1,u)))));
%eval(sprintf('all_frac%g=frac%g;',str2double(followed_masses(1,u)),str2double(followed_masses(1,u))));
frac34=averaged_integrals(:,col34)/AO2;
frac36=averaged_integrals(:,col36)/AO2;
%end

frac32=1-frac34-frac36;

if isempty(col32)==0
    frac32_exp=averaged_integrals(:,col32)/AO2;
  %  all_frac32./frac32_exp
   all_frac32=frac32_exp;
end

f18=frac34/2+frac36;

f18_2=(averaged_integrals(:,col34)/2+averaged_integrals(:,col36))/AO2;

f18_RT=f18(rowRT);

%disp('f18_RT/f18_GPA')
f18_RT/f18_GPA;

R0=-molar_flow_rate/S*log(f18/f18_RT);
TTT=1000./(temperatures+273.15);

k=R0/c_0;


content_annotation={sprintf('Sample : %s',sample_name),sprintf('Date & time : %s',date),sprintf('m = %g mg , M = %g g/mol',m_mg,M),sprintf('\\rho = %g g/cm³ , SSA = %g m²/g',rho,SSA),sprintf('bed length = %g mm , rel. density = %0.0f %%',bed_length,bed_density*100),sprintf('p(O2) = %g bar , flow = %g mL/min',pO2,flow_rate),sprintf('pulse vol. = %g \\mum , inner diam. = %g mm',pulse_volume,ID),sprintf('res. time = %0.2f s , pulse length = %0.0f mm',residence_time,pulse_length),sprintf(''),sprintf('n_O^{pulse}/n_O^{bed} = %0.3f',n_O_pulse/n_O_bed),sprintf('gas ref. : %s , f^{thermo}_{18}/f^{GPA}_{18} = %0.3f ',labelled_ref,f18_RT/f18_GPA)};
isotopologues_fraction=figure('Name','Isotopologues fraction','NumberTitle','off');
plot(temperatures(1:end),frac32(1:end),'bs','Linewidth',2);
%errorbar(temperatures(2:end),frac32(2:end),err)
hold on;
plot(temperatures(1:end),frac34(1:end),'s','color',[1.0  0.0  1.0],'Linewidth',2)
plot(temperatures(1:end),frac36(1:end),'rs','Linewidth',2)
xlabel(sprintf('Temperature (%cC)',char(176)))
ylabel('Oxygen isotopologues fraction')
title('Experimental summary')
axis([-Inf Inf 0 1])
annotation('textbox',[0.15 0.45 0.3 0.3],'String',content_annotation,'FitBoxToText','on')
legend('f32','f34','f36','Location','eastoutside')


Arrhenius_plot=figure('Name','Provisional Arrhenius plot','NumberTitle','off');
semilogy(TTT, k, 'r--s', 'LineWidth', 2);
xlabel('1000/T (K^{-1})');
ylabel('k (m.s^{-1})');

message = sprintf('Select the points of interest to extract the activation energies (select two diagonally opposite corners of a rectangular box).');
uiwait(msgbox(message));
area_of_interest=ginput(2);
close(Arrhenius_plot);
TTT_min=min(area_of_interest(:,1));
TTT_max=max(area_of_interest(:,1));
k_min=min(area_of_interest(:,2));
k_max=max(area_of_interest(:,2));
TTT_to_delete = TTT < TTT_min | TTT > TTT_max;
all_TTT=TTT;
all_temperatures=temperatures;
all_k=k;
all_R0=R0;
all_f18=f18;

% for u=1:size(followed_masses,2)
%     eval(sprintf('all_frac%g=frac%g;',followed_masses(1,u),followed_masses(1,u)))
all_frac32=frac32;
all_frac34=frac34;
all_frac36=frac36;
%end
TTT(TTT_to_delete) = [];
temperatures(TTT_to_delete) = [];
k(TTT_to_delete) = [];
R0(TTT_to_delete) = [];
f18(TTT_to_delete) = [];
frac32(TTT_to_delete) = [];
frac34(TTT_to_delete) = [];
frac36(TTT_to_delete) = [];

%Getting the activation energy
log_k=log10(k);
[polynom,quality]=polyfit(TTT,log_k,1);
Arrhenius_fit=10.^(polynom(2)+polynom(1)*TTT);
EA=-polynom(1)*R*log(10);
EA_eV=EA/FarC*1000;
PreExp=10^polynom(2);

%Pre-scaling the Arrhenius plot
k_min=log10(k_min);
k_max=log10(k_max);
TTT_min=round(mean(TTT_min*10))/10-0.1;
TTT_max=round(mean(TTT_max*10))/10+0.1;
TTT_n=round(mean((TTT_max-TTT_min)/0.1))+1;

content_annotation={sprintf('E_a: %0.2g eV',EA_eV),sprintf('Pre-exp. factor: %0.3g',PreExp)};
Arrhenius_plot=figure('Name','k','NumberTitle','off');
% setup bottom axis
ax = axes();
hold(ax);
ax.YAxis.Scale = 'log';
xlabel(ax, '1000/T (K^{-1})');
ylabel(ax, 'k (m.s^{-1})');

% setup top axis
ax_top = axes(); % axis to appear at top
hold(ax_top);
ax_top.XAxisLocation = 'top';
ax_top.YAxisLocation = "right";
ax_top.YTick = [];
%ax_top.XDir = 'reverse';
ax_top.Color = 'none';
xlabel(ax_top, 'Temperature (°C)');

% linking axis
linkprop([ax, ax_top],{'Units','Position','ActivePositionProperty'});
ax.Position(4) = ax.Position(4) * .95;
TT = linspace(TTT_min,TTT_max,TTT_n);

% configure limits of bottom axis
ax.XLim = [TT(1) TT(end)];
ax.XTick = TT;
ax.XAxis.TickLength = [0.015, 0.00];
ax.YAxis.TickLength = [0.02, 0.00];
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = linspace(TTT_min,TTT_max,2*TTT_n-1);

% configure limits and labels of top axis
temp_ticks = [100 200 300 400 500 600 700 800 900 1000];
TTT_tick = 1000./(temp_ticks+273.15);
ax_top.XLim = [TT(1) TT(end)];
ax_top.XTick = fliplr(TTT_tick);
ax_top.XTickLabel = compose('%1.0f', fliplr(temp_ticks));
ax_top.XAxis.TickLength = [0.02, 0.01];
ax_top.XAxis.MinorTick = 'on';
plot(ax, all_TTT, all_k, 'b--s', 'LineWidth', 1);
plot(ax, TTT, Arrhenius_fit, 'r-', 'LineWidth', 2);
legend(ax,'data','fit','Location','south')
annotation('textbox',[0.6 0.45 0.3 0.3],'String',content_annotation,'FitBoxToText','on')


fitting_data={temperatures,R0,f18,frac32,frac34,frac36,f18_RT,f36_RT,S,molar_flow_rate,1,size(temperatures,1)};
guess=0.8;
p_calc=fminsearch('fitPIE',guess,options_fit);
guess=ones(size(temperatures,1),1)*p_calc;
p_calc=fminsearch('fitPIE',guess,options_fit);

f36c=(1-p_calc).^2*f18_RT^2.*exp(-2.*R0*S/molar_flow_rate)./(1-2.*p_calc)+f36_RT.*exp(-R0.*S./(molar_flow_rate.*p_calc))-(1-p_calc).^2*f18_RT^2.*exp(-R0.*S./(molar_flow_rate.*p_calc))./(1-2.*p_calc);
f34c=2*(f18-f36c);
f32c=1-f34c-f36c;


twostepsmodelfitsimple=figure('Name','Two steps model fitting','NumberTitle','off');
plot(all_temperatures,all_frac32,'bs','Linewidth',2);
hold on;
plot(temperatures,f32c,'b:','Linewidth',2)
plot(all_temperatures,all_frac34,'s','color',[1.0  0.0  1.0],'Linewidth',2)
plot(temperatures,f34c,':','color',[1.0  0.0  1.0],'Linewidth',2)
plot(all_temperatures,all_frac36,'rs','Linewidth',2)
plot(temperatures,f36c,'r:','Linewidth',2)
xlabel(sprintf('Temperature (%cC)',char(176)))
ylabel('Oxygen isotopologues fraction')
axis([min(all_temperatures)-20 max(all_temperatures)+20 0 1])
title('Fitting of p')
legend('f32','fit','f34','fit','f36','fit','Location','best')


Rads=R0./p_calc;
Rinc=R0./(1-p_calc);

temper=(0:10:1100)';
TTTT=1000./(temper+273.15);

log_R0=log10(R0);
[polynom,quality]=polyfit(TTT,log_R0,1);
[R0c,delta] = polyval(polynom,TTT,quality);
R0c=10.^(polynom(2)+polynom(1)*TTTT);
EA_R0=-polynom(1)*R*log(10);
EA_eV_R0=EA_R0/FarC*1000;
stdev_R0=mean(delta)*R*log(10)/FarC*1000; % to verify

log_Rads=log10(Rads);
[polynom,quality]=polyfit(TTT,log_Rads,1);
[Radsc,delta] = polyval(polynom,TTT,quality);
Radsc=10.^(polynom(2)+polynom(1)*TTTT);
EA_Rads=-polynom(1)*R*log(10);
EA_eV_Rads=EA_Rads/FarC*1000;
stdev_Rads=mean(delta)*R*log(10)/FarC*1000;% to verify

log_Rinc=log10(Rinc);
[polynom,quality]=polyfit(TTT,log_Rinc,1);
[Rincc,delta] = polyval(polynom,TTT,quality);
Rincc=10.^(polynom(2)+polynom(1)*TTTT);
EA_Rinc=-polynom(1)*R*log(10);
EA_eV_Rinc=EA_Rinc/FarC*1000;
stdev_Rinc=mean(delta)*R*log(10)/FarC*1000; % to verify

%pre-scaling the Arrhenius plot
R_min=10^(round(min(log10(R0))-0.5))/1;%*8/10
R_max=10^(round(max(log10(Rinc))+0.5))*1;


content_annotation={sprintf('E_a(\\Re_0) = %0.2g +/- %0.1g eV',EA_eV_R0,stdev_R0),sprintf('E_a(\\Re_{ads}) = %0.2g +/- %0.1g eV',EA_eV_Rads,stdev_Rads),sprintf('E_a(\\Re_{inc}) = %0.2g +/- %0.1g eV',EA_eV_Rinc,stdev_Rinc)};
rates=figure('Name','Rates','NumberTitle','off');
% setup bottom axis
ax = axes();
hold(ax);
ax.YAxis.Scale = 'log';
xlabel(ax, '1000/T (K^{-1})');
ylabel(ax, sprintf('\\Re_0, \\Re_{ads} and \\Re_{inc} (mol.m^{-2}.s^{-1})'));

% setup top axis
ax_top = axes(); % axis to appear at top
hold(ax_top);
ax_top.XAxisLocation = 'top';
ax_top.YAxisLocation = "right";
ax_top.YTick = [];
%ax_top.XDir = 'reverse';
ax_top.Color = 'none';
xlabel(ax_top, 'Temperature (°C)');

% linking axis
linkprop([ax, ax_top],{'Units','Position','ActivePositionProperty'});
ax.Position(4) = ax.Position(4) * .95;
TT = linspace(TTT_min,TTT_max,TTT_n);

% configure limits of bottom axis
ax.XLim = [TT(1) TT(end)];
%%ax.YLim = [R_min R_max];
ax.XTick = TT;
ax.XAxis.TickLength = [0.015, 0.00];
ax.YAxis.TickLength = [0.02, 0.00];
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = linspace(TTT_min,TTT_max,2*TTT_n-1);

% configure limits and labels of top axis
temp_ticks = [100 200 300 400 500 600 700 800 900 1000];
TTT_tick = 1000./(temp_ticks+273.15);
ax_top.XLim = [TT(1) TT(end)];
ax_top.XTick = fliplr(TTT_tick);
ax_top.XTickLabel = compose('%1.0f', fliplr(temp_ticks));
ax_top.XAxis.TickLength = [0.02, 0.01];
ax_top.XAxis.MinorTick = 'on';
plot(ax, TTT, R0, 'ks', 'LineWidth', 1);
plot(ax, TTT, Rads, 'bs', 'LineWidth', 1);
plot(ax, TTT, Rinc, 'gs', 'LineWidth', 1);
plot(ax, TTTT, R0c, 'k:', 'LineWidth', 1);
plot(ax, TTTT, Radsc, 'b:', 'LineWidth', 1);
plot(ax, TTTT, Rincc, 'g:', 'LineWidth', 1);
legend(ax,sprintf('\\Re_0'),sprintf('\\Re_{ads}'),sprintf('\\Re_{inc}'),'Location','east');
annotation('textbox',[0.45 0.55 0.3 0.3],'String',content_annotation,'FitBoxToText','on');
%[0.2 0.05 0.3 0.3]

pc=R0c./Radsc;
f18c=f18_RT.*exp(-(R0c*S)/molar_flow_rate);
% from the excel sheet... with an error: frac36c=(1-pc).^2*f18_RT^2.*exp(-2*R0c*S/molar_flow_rate)./(1-2*pc)+f18_RT.*exp(-R0c*S./(molar_flow_rate*pc))-(1-pc).^2*f36_RT^2.*exp(-R0c*S./(molar_flow_rate*pc))./(1-2*pc);
frac36c=(1-pc).^2*f18_RT^2.*exp(-2*R0c*S/molar_flow_rate)./(1-2*pc)+f36_RT.*exp(-R0c*S./(molar_flow_rate*pc))-(1-pc).^2*f18_RT^2.*exp(-R0c*S./(molar_flow_rate*pc))./(1-2*pc);
frac34c=2*(f18c-frac36c);
frac32c=1-frac34c-frac36c;

twostepsmodelfit=figure('Name','Two steps model fitting','NumberTitle','off');
plot(all_temperatures,all_frac32,'bs','Linewidth',2);
hold on;
plot(temper,frac32c,'b:','Linewidth',1);
plot(all_temperatures,all_frac34,'s','color',[1.0  0.0  1.0],'Linewidth',2)
plot(temper,frac34c,':','color',[1.0  0.0  1.0],'Linewidth',1)
plot(all_temperatures,all_frac36,'rs','Linewidth',2)
plot(temper,frac36c,'r:','Linewidth',1)
xlabel(sprintf('Temperature (%cC)',char(176)))
ylabel('Oxygen isotopologues fraction')
axis([min(all_temperatures)-20 max(all_temperatures)+20 0 1])
title('Simulation based on activation energies')
legend('f32','fit','f34','fit','f36','fit','Location','best')
annotation('textbox',[0.15 0.05 0.3 0.3],'String',content_annotation_ref,'FitBoxToText','on')




choice = questdlg('Save the figures?', ...
    'Notification', ...
    'yes','no','yes');
switch choice
    case 'no'
    case 'yes'
        message=msgbox('Saving figures, please wait...','Notification');
        
        if exist(saving_directory,'dir') == 7
        else
            mkdir(saving_directory);%create a directory for each fit
        end
        if strcmp(figure_format,'png')==1
            print('-dpng',fig_res,twostepsmodelfitsimple,fullfile(saving_directory,'simulation_from_fitting_p.png'));
            print('-dpng',fig_res,twostepsmodelfit,fullfile(saving_directory,'simulation_from_activation energies.png'));
            print('-dpng',fig_res,Arrhenius_plot,fullfile(saving_directory,'Arrhenius_k.png'));
            print('-dpng',fig_res,isotopologues_fraction,fullfile(saving_directory,'summary.png'));
            print('-dpng',fig_res,rates,fullfile(saving_directory,'rates.png'));
        elseif strcmp(figure_format,'jpg')==1
            saveas(isotopologues_fraction,fullfile(saving_directory,'summary.jpg'),'jpg');
            saveas(twostepsmodelfit,fullfile(saving_directory,'simulation_from_activation energies.jpg'),'jpg');
            saveas(twostepsmodelfitsimple,fullfile(saving_directory,'simulation_from_fitting_p.jpg'),'jpg');
            saveas(Arrhenius_plot,fullfile(saving_directory,'Arrhenius_k.jpg'),'jpg');
            saveas(rates,fullfile(saving_directory,'rates.jpg'),'jpg');
        end
        if strcmp(export_fig,'yes')==1
            FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
            for iFig = 2:length(FigList)
                FigHandle = FigList(iFig);
                FigName   = get(FigHandle, 'Name');
                savefig(FigHandle, fullfile(saving_directory, FigName, '.fig'));
                %print('-dpng',fig_res,FigHandle,strcat(saving_directory,sprintf('%s.png',FigName)));
            end
        end
        delete(message);
        %close all;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% exporting the vectors %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
choice = questdlg('Export the fitted values and data?', ...
    'Saving', ...
    'yes','no','yes');
switch choice
    case 'no'
    case 'yes'
        if exist(saving_directory,'dir') == 7
        else
            mkdir(saving_directory);%create a directory for each fit
        end
        
        
        
        %% saving the integration parameters
        integration_log_list={'default_back_step'  'default_back_step';
            'default_integ_length' 'default_integ_length';
            'background_type' 'background_type';
            'RB_type_override' 'RB_type_override';
            'deflt_LB_drv_trshld' 'default_LB_drv_trshld';
            'deflt_RB_drv_trshld' 'default_RB_drv_trshld'};
        integ_log_path=fullfile(saving_directory,'integration_log.csv');
        savloarch_log(integration_log_list,integ_log_path,'save');
        
        mass_set_log_path=fullfile(saving_directory,'masses-settings_log.csv');
        writematrix(masses_settings, mass_set_log_path);
        
        %% saving the integrals
        integrals_log_path=fullfile(saving_directory,'integrals_log.csv');
        writematrix(integrals_log, integrals_log_path);
        
        
         %% saving k and R0 at all temperatures
         %possible vectors to save and names:
         potential_data_to_save={ 'T' 'all_temperatures';
             '1000/T' 'all_TTT';
             'f18' 'all_f18' ;
             'f32' 'all_frac32';
             'f34' 'all_frac34' ;
             'f36' 'all_frac36';
             'stdev_34' 'stdev_34';
             'stdev_36' 'stdev_36';
             'err_34' 'err_34';
             'err_36' 'err_36';
             'k' 'all_k';
             'R0' 'all_R0'};
         export_fit_path=fullfile(saving_directory,'k and R0.csv');
         analysis_log(potential_data_to_save,export_fit_path);
         
       %% saving p fit
         %possible vectors to save and names:
         potential_data_to_save={'T'  'temperatures' ;
             '1000/T' 'TTT' ;
             'f18' 'f18';
             'f32' 'frac32';
             'f34' 'frac34' ;
             'f36' 'frac36';
             'k' 'k';
             'R0' 'R0';
             'Rads' 'Rads';
             'Rinc' 'Rinc';
             'p' 'p_calc';
             'k_calc' 'Arrhenius_fit'};
         export_fit_path=fullfile(saving_directory,'p_fit.csv');
         analysis_log(potential_data_to_save,export_fit_path);
         
  
         
         %% saving the simulation
         %possible vectors to save and names:
         potential_data_to_save={ 'T' 'temper' ;
             '1000/T' 'TTTT' ;
             'f18' 'f18c';
             'f32' 'frac32c';
             'f34' 'frac34c' ;
             'f36' 'frac36c';
             'R0' 'R0c';
             'Rads' 'Radsc';
             'Rinc' 'Rincc';
             'p' 'pc' };
         export_fit_path=fullfile(saving_directory,'simulation.csv');
         analysis_log(potential_data_to_save,export_fit_path);
         
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% exporting fitted parameters %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %saves the experimental values and fitted parameters in the log
        %file in the local log folder; these can be saved to the main
        %database later
        % 1st row: names in the log file, 2nd row: variables' names
        
        followed_masses=strjoin(followed_masses);
        
        potential_parameters_to_save={'timestamp' 'timestamp';
             'sample'  'sample_name';
             'EA(R0) [eV]' 'EA_eV_R0';
             'EA(Rads) [eV]' 'EA_eV_Rads';
             'EA(Rinc) [eV]' 'EA_eV_Rinc';
             'EA(R0) [kJ/mol]' 'EA_R0';
            'EA(Rads) [kJ/mol]' 'EA_Rads' ;
             'EA(Rinc) [kJ/mol]' 'EA_Rinc';
             'stdev(R0) [eV]' 'stdev_R0';
            'stdev(Rads) [eV]' 'stdev_Rads';
            'stdev(Rinc) [eV]' 'stdev_Rinc';
              'm [mg]'  'm_mg';
            'M [g/mol]' 'M' ;
            'rho [g/cm3]' 'rho';
            'SSA [m2/g]' 'SSA' ;
            'pO2 [bar]' 'pO2';
            'MFC_T [degC]'  'MFC_T';
            'bed_length [mm]' 'bed_length';
            'flow_rate [mL/min]' 'flow_rate';
             'pulse_vol [uL]' 'pulse_volume';
            'followed_masses' 'followed_masses';
            'ref_mass' 'ref_mass';
            'f34' 'f34_RT';
            'f36' 'f36_RT';
              'reference' 'labelled_ref';
             'xp_32' 'use_experimental_32' ;
             'inner_diam. [mm]' 'ID';
             'default_back_step'  'default_back_step';
             'default_integ_length' 'default_integ_length';
             'background_type' 'background_type';
             'RB_type_override' 'RB_type_override';
             'deflt_LB_drv_trshld' 'default_LB_drv_trshld';
             'deflt_RB_drv_trshld' 'default_RB_drv_trshld'};
        
savloarch_log(potential_parameters_to_save,fullfile(saving_path,'fitting_log.csv'),'archive');

       
end
fclose all;

%%%%%%%%%%%%%%%%%%%%%%%
%% INTERNAL ROUTINES %%
%%%%%%%%%%%%%%%%%%%%%%%

% 
% function save_conditions
% global experimental_conditions_log_path sample_name data_source time_unit M m_mg rho SSA pO2 MFC_T bed_length pulse_volume flow_rate followed_masses ref_mass  f34_RT f36_RT labelled_ref use_experimental_32 ID
% experimental_conditions_log={sample_name, data_source, time_unit, m_mg, M, rho, SSA, pO2, MFC_T, bed_length, pulse_volume, flow_rate, followed_masses, ref_mass, f34_RT, f36_RT, labelled_ref, use_experimental_32, ID};
% experimental_conditions_log_columns={'sample' 'data_source' 'time_unit' 'm' 'M'  'rho' 'SSA' 'pO2' 'MFC_T' 'bed_length' 'pulse_vol' 'flow_rate' 'followed_masses' 'ref_mass' 'f34' 'f36' 'reference' 'xp_32' 'inner_diam'};
% experimental_conditions_log=cell2table(experimental_conditions_log);
% experimental_conditions_log.Properties.VariableNames = experimental_conditions_log_columns;
% writetable(experimental_conditions_log,experimental_conditions_log_path,'Delimiter',',');
% end
