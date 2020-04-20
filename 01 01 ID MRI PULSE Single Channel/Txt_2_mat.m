close all
clear
clc
%% opening the text file and replacing "," with "." 
directory = 'H:\Montazeri\DATA_MRI';
matfiles = dir(fullfile(directory, '*.txt'));
nfiles = length(matfiles);
data  = cell(nfiles);
for kk = 1 : nfiles
Input = fopen( fullfile(directory, matfiles(kk).name) );
Input2 = fopen(matfiles(kk).name,'w+'); %open file for writing
i = 1;
tline = fgetl(Input);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(Input);
    A{i} = tline;
end

 modified_A=cell(1,numel(A));

for i=1: numel(A)
    modified_A{i}= strrep(A{i}, ',', '.');
end

fprintf(Input2,'%s\r\n', modified_A{:});

fclose(Input);
fclose(Input2);
%% Extracting column data

Input=fopen(matfiles(kk).name,'r');

i=1;
nextline=fgetl(Input);
while ~feof(Input)    
    if sscanf(nextline,'%s %s %s %s')
        column_cell{i} = sscanf(nextline,'%f %f %f %f');
        i=i+1;
    end 
    nextline=fgetl(Input);
end
data_matrix=cell2mat(column_cell)';

%% Saving the output

savefile=matfiles(kk).name(1:end-4);
eval(sprintf('%s=%s',savefile,'data_matrix'))
save(savefile,savefile);

fclose('all');
delete(matfiles(kk).name)

end












