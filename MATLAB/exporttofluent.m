function [] = exporttofluent(inputData,varNames,fileName)
%EXPORTTOFLUENT Summary of this function goes here
%   - inputaData: matrix with a atructure: x y z Value1 Value2...
%   - varNames: name of the variabel in FLUENT as a cell structure to
%   overwrite (within FLUENT display/contour/contours of to view 
%   available variables. In order to display cells as one compartment,
%   the interpolation data should be an UDS
%   - fileName: Name of the file, which is a *.txt-file ('filename.txt').
%   After the file is created, it must be changed to a *.ip file by
%   altering the file format manually from *.txt to .*ip

%% Write header
fid = fopen(fileName,'w');
% first line: interpolation file version (1,2,3,4 or 5)
% second line: dimension (2 or 3)
% third line: number of points/cells
% fourth line: number field variables
fprintf(fid, '3\n3\n%i\n%i\n', size(inputData,1),size(inputData,2)-3);
fclose(fid);
% Loop to write the variable names
for i = 1:size(inputData,2)-3
    fid = fopen(fileName,'a');
    fprintf(fid,'%s\n',varNames{i});
    fclose(fid);    
end

%% write main
% Loop over all columns of inputData. Every column of inputData (x,y,z,...)
% is a single section
for i = 1:size(inputData,2)

    fid = fopen(fileName,'a');
    % Opening of the section
    fprintf(fid,'(');
    fclose(fid);
    writematrix(inputData(:,i),fileName,'WriteMode','append')
    fid = fopen(fileName,'a');
    % Closing of section (depends on first line of header file, see 
    % FLUENT Users Guide Format of the Interpolation File)
    fprintf(fid,')');
    fclose(fid);
end

end

