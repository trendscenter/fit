function ica_fuse_writeDupli(filename,Duplication)

fid = fopen(filename,'a');
dlm=',';
for i=1:50000
    dlm=[dlm,','];
end
for i=1:length(Duplication)

    y=Duplication{i};
    y=str2mat(y)';[r,c]=size(y);
    y=[y;dlm(1:c)];
    fprintf(fid,'%s \n',y);
    
    
end

fclose(fid)
disp('file saved')