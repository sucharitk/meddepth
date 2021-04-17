function add_header_csv(filename, header_text)
%
% function add_header_csv(filensame, header_text)
%
% for csv files exported ?rom matlab to do analayses in R add the column
% names at the top
%

ft = fileread(filename);
ft = [header_text,newline,ft];
FID = fopen(filename, 'w');
if FID == -1, error('Cannot open file %s', filename); end
fwrite(FID, ft, 'char');
fclose(FID);

end