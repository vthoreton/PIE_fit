function [new_file] = reconstruct_data_file(path)
fid= fopen(sprintf('%s',path),'r'); %open the file
k=1;  
while ~feof(fid) %while the end of the file is not reached,
l=fgetl(fid);%read each line after an other
n_col_i=sum(l==sprintf('\t')); % count delimiters (tabs) and calculate the number of columns
structure(k,1)=n_col_i;%record the number of columns for each line
 k=k+1;  %increment k for next line
end
n_col=max(structure);% number of columns
fid= fopen(sprintf('%s',path),'r');%open the file from the beginning again
file=fread(fid,'*char')'; %read the content to a working file (char array)
fclose(fid);
file = strrep(file,sprintf('\t'),';DELIMITER'); %changes the delimiters to comma and add some content to avoid loosing the delimiters of the structure 
file=textscan(file, '%s', 'Delimiter', sprintf('\r\n'));
file=file{1};
for i = 1:size(structure,1)
    if structure(i)<n_col
        for ii = 1:n_col-structure(i)
           file{i}=sprintf('%s%sDELIMITER',file{i},sprintf(';')); %removes the temporary content...
        end
    end
end
for i=1:size(structure,1)
new_file(i,:) = strsplit(file{i},';');
end
new_file=new_file(:,1:end-1);
new_file=strrep(new_file,'DELIMITER','');%removes the temporary content...
end