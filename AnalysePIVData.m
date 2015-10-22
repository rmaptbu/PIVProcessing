function [vel_x,true_v_x]=AnalysePIVData(filename)
%Outputs the x- projection of x-velocities of all files
%Format: vel_x(n_interrogation,type,frame,flowratio)
%n_interrogation=number of interrogation windows in the y-direction of PIV
%type=type of data of projection: [mean,max,min,diff(mean,true value)]
%frame=which frame of images, sorted by entropy
%flowratio=percentage of moving circles vs. stationary ones

%Read data from file
data=matfile(filename);
v=data.resultslist(7:8,:).';
names=data.filename(2:2:end,1);
entropies_str=regexp(names,'_s\d.\d+','match');
particle_numbers_str=regexp(names,'_N\d.\d+','match');
movement_ratios_str=regexp(names,'_r\d.\d+','match');
entropies=[];
particle_numbers=[];
movement_ratios=[];
for i=1:size(v,1)
    entropies_str(i)=(regexp(entropies_str{i},'\d.\d+','match'));
    particle_numbers_str(i)=(regexp(particle_numbers_str{i},'\d.\d+','match'));
    movement_ratios_str(i)=(regexp(movement_ratios_str{i},'\d.\d+','match'));
    
    entropies=cat(1,entropies,str2double(entropies_str{i}));
    particle_numbers=cat(1, particle_numbers,str2double( particle_numbers_str{i}));
    movement_ratios=cat(1,movement_ratios,str2double(movement_ratios_str{i}));
end

%Project velocity data along the x-axis (keeping mean,min,max)
vel_x=zeros(size(v{1},1),3,size(v,1));
vel_y=zeros(size(v{1},1),3,size(v,1));
for i=1:size(v,1)
    vel_x(:,1,i)=-mean(v{i,1},2);
    vel_x(:,2,i)=-min(v{i,1},[],2);
    vel_x(:,3,i)=-max(v{i,1},[],2);
    vel_y(:,1,i)=-mean(v{i,2},2);
    vel_y(:,2,i)=-min(v{i,1},[],2);
    vel_y(:,3,i)=-max(v{i,1},[],2);
end

%Sort velocities by entropy
[entropies,mapping]=sort(entropies);
vel_x_sort=vel_x;
vel_y_sort=vel_y;
particle_numbers_sort=particle_numbers;
movement_ratios_sort=movement_ratios;
for i=1:size(mapping,1)
    vel_x_sort(:,:,i)=vel_x(:,:,mapping(i));
    vel_y_sort(:,:,i)=vel_y(:,:,mapping(i));
    particle_numbers_sort(i)=particle_numbers(mapping(i));
    movement_ratios_sort(i)=movement_ratios(mapping(i));
end
vel_x=vel_x_sort;
vel_y=vel_y_sort;
particle_numbers=particle_numbers_sort;
movement_ratios=movement_ratios_sort;

%Create grant truth
true_v_x=zeros(size(v{1},1),1);
pos = 1:256;
true_shape=10*(pos/256); 
for i=1:15
    true_v_x(i)=mean(true_shape(16*i:16*(i+1)));
end

%split data by mixing ratios (flow/stationary)
%ignore the s_xx
vel_x_1=[];
s_95=[];
vel_x_2=[];
s_90=[];
vel_x_3=[];
s_80=[];
vel_x_4=[];
s_50=[];
vel_x_5=[];
for i=1:size(v,1)
    if movement_ratios(i)==0.99
       vel_x_1=cat(3,vel_x_1,vel_x(:,:,i));
       s_95=cat(1,s_95,entropies(i));
    elseif movement_ratios(i)==0.95
       vel_x_2=cat(3,vel_x_2,vel_x(:,:,i));
       s_90=cat(1,s_90,entropies(i));
    elseif movement_ratios(i)==0.90
       vel_x_3=cat(3,vel_x_3,vel_x(:,:,i));
       s_80=cat(1,s_80,entropy(i));
    elseif movement_ratios(i)==0.80
       vel_x_4=cat(3,vel_x_4,vel_x(:,:,i));
       s_50=cat(1,s_50,entropies(i));
    elseif movement_ratios(i)==0.60
       vel_x_5=cat(3,vel_x_5,vel_x(:,:,i));
       s_50=cat(1,s_50,entropies(i));
    end
end
vel_x=cat(4,vel_x_5,vel_x_4,vel_x_3,vel_x_2,vel_x_1);

%Subtract true value from measurement
for j=1:size(vel_x,3)
    for r=1:size(vel_x,4)
        vel_x(:,4,j,r)=vel_x(:,1,j,r)-true_v_x;
    end
end




