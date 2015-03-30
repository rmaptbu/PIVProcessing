%% Write 10s

clear;clc;
% Agg
basepathHD='E:\280712\Pbs\';
basepath='C:\Users\J-\Desktop\280712\Pbs\';
Cases=[2 7 8 9 5 19 24 40]; % Agg
%% PBS
% basepathHD='E:\280712\Pbs\';
% basepath='C:\Users\J-\Desktop\280712\Pbs\';
% Cases=[3 4 7 8 9 17 21 31]; % Agg

basepathS=[basepathHD,'Pre\Str\'];    % strobe write path

writepathS=[basepath,'Pre\StrF\'];    % strobe write path

cd(basepathS)
cases=struct2cell(dir);
cases=cases(1,3:end);


for i=1:length(cases)%Cases   % for each case
    disp(['Case ', num2str(i)])
    basepathScase=[basepathS,char(cases(i)),'\'];
    cd(basepathScase)
    filesS=struct2cell(dir);
    filesS=filesS(1,3:end);   
    writepathScase=[writepathS,char(cases(i)),'\'];mkdir(writepathScase);

    for j=1:length(filesS)
        
        disp(['Write image : Case ', num2str(i),', Image ' num2str(j)])
        ImS=double(imread([basepathScase,char(filesS(j))]));
     
        ImS=uint16([ImS(:,481:end) zeros(950,32)]);

        imwrite(ImS,[writepathScase,char(filesS(j))])
        imwrite(ImL,[writepathLcase,char(filesL(j))])
    end
end