function file_type =  SPM_smooth(filenames,sizefwhm)

% listofnames = [filenames '*' '.nii'];

% list = ls(listofnames);

list = cellstr(filenames);

for i = 1:size(list,1)
    
    spm_smooth(list{i,:},['s' num2str(sizefwhm) list{i,:}],sizefwhm);
    file_type{i} = ['s' num2str(sizefwhm) list{i,:}];
    
end