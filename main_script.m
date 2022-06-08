clear;





% parfor i=1:36
%     pause(0.0001*i)
%     cmd=['set path=%path:C:\Program Files\MATLAB\R2020b\bin\win64;=% & bin\Release\frz_lattice_model.exe 2.6 0.1 16 ',num2str(i)];
%     system(cmd,'-echo');%outputs to current matlab directory
% end


%%
for i=1:1000
    j=1;
    file=fopen([num2str(i),'.txt']);
    tline=fgetl(file);
    while ~isempty(tline) && ischar(tline)
        z=textscan(tline,'%f','Delimiter','\t');%outputs cell array
        t(j)=z{1}(1);
        c{i,j}=z{1}(2:end);
        tline=fgetl(file);
        j=j+1;
    end

fclose(file);

end

save('results.mat','c','t');


%%



