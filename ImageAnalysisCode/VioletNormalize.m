close all
clear all

time = [0,6,12,18,24];
strain1 = [192]; %Also, MR style
sample = {{[1,2],[1,2],[1,2,3],[1,2,3],[1,2,3]}};
Date1 = {'July_18_2020_1'};
red_channel = 1;
green_channel = 0;
blue_channel = 0;
violet_channel = 2;

discard_violet = 1;
violet_cut = 0;
violet_back = 2.2199e+03;

for s_l = 1:length(strain1)
    strain_i = strcat(['MR',num2str(strain1(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date1{s_l}]);
    
    violet = [];
    
    for j = 1:length(time)
        t = strcat(['t',num2str(time(j))]);
        close all
        
        timeDir = strcat([dateDir,'/',t]);
        d2 = dir([timeDir, '/*.mat']);
        %samples = length(d2);
        samples = sample{s_l}{j};
        file = strcat([timeDir,'/']);
        
        
        
        for k = samples
            
            s1 = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
            load(s1)
            clearvars part4_V_Normalized
            save(s1)
            
            
        end
    end
    
end

for s_l = 1:length(strain1)
    strain_i = strcat(['MR',num2str(strain1(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date1{s_l}]);
    
    violet = [];
    
    for j = 1:length(time)
        t = strcat(['t',num2str(time(j))]);
        close all
        
        timeDir = strcat([dateDir,'/',t]);
        d2 = dir([timeDir, '/*.mat']);
        %samples = length(d2);
        samples = sample{s_l}{j};
        file = strcat([timeDir,'/']);
        
        
        
        for k = samples
            
            s1 = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
            load(s1,'part4')

            if violet_channel==1
                
                samples = sample{s_l}{j};
                
                for l = 1:length(part4)
                    violet = [violet part4(l).Intensity_One];
                end
                violet(violet<0)=0;
                
            elseif violet_channel == 2
                
                for l = 1:length(part4)
                    violet = [violet part4(l).Intensity_Two];
                end
                violet(violet<0)=0;
                
            elseif violet_channel == 3
                samples = sample{s_l}{j};
               
                for l = 1:length(part4)
                    violet = [violet part4(l).Intensity_Three];
                end
                violet(violet<0)=0;
               
                
            elseif violet_channel == 4
                samples = sample{s_l}{j};
                
                for l = 1:length(part4)
                    violet = [violet part4(l).Intensity_Four];
                end
                violet(violet<0)=0;
                
            end
        end
    end
    
end

clearvars part4 part4_V_Normalized

for s_l = 1:length(strain1)
    strain_i = strcat(['MR',num2str(strain1(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date1{s_l}]);
    
    violet = [];
    
    for j = 1:length(time)
        t = strcat(['t',num2str(time(j))]);
        close all
        
        timeDir = strcat([dateDir,'/',t]);
        d2 = dir([timeDir, '/*.mat']);
        %samples = length(d2);
        samples = sample{s_l}{j};
        file = strcat([timeDir,'/']);
        
        
        for k = samples
            
            s1 = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
            varlist = who;
            varlist1 = varlist;
            varlist =strjoin(varlist','$|');
            load(s1,'-regexp', ['^(?!' varlist ')\w']);
            clearvars part4_V_Normalized
            violet_norm = [];
            
            part4V_length = length(part4);
            
            if discard_violet == 1
                violetvolume = zeros(part4V_length,1);
                for ni = 1:part4V_length
                    if violet_channel == 1
                        violetvolume(ni) = part4(ni).Intensity_One/(5*part4(ni).Volume);
                    elseif violet_channel == 2
                        violetvolume(ni) = part4(ni).Intensity_Two/(5*part4(ni).Volume);
                    elseif violet_channel == 3
                        violetvolume(ni) = part4(ni).Intensity_Three/(5*part4(ni).Volume);
                    elseif violet_channel == 4
                        violetvolume(ni) = part4(ni).Intensity_Four/(5*part4(ni).Volume);
                    end
                end
                
                [violetvolume,vindex] = sort(violetvolume);
                
                vindex(violetvolume<violet_back) = [];
                violetvolume(violetvolume<violet_back) = [];
                
                
                
                
                
                for n = 1:length(vindex)
                    
                    
                    if violet_channel == 1
                        violet = [violet part4(vindex(n)).Intensity_One];
                    elseif violet_channel == 2
                        violet = [violet part4(vindex(n)).Intensity_Two];
                    elseif violet_channel == 3
                        violet = [violet part4(vindex(n)).Intensity_Three];
                    elseif violet_channel == 4
                        violet = [violet part4(vindex(n)).Intensity_Four];
                    end
                end
                
%                 [violetvolume,vindex] = sort(violetvolume);
%                 v_cut = int32(violet_cut*part4V_length);
%                 violetvolume(1:v_cut) = [];
%                 vindex(1:v_cut) = [];
%                 
%                 
%                 
%                 
%                 for n = 1:length(vindex)
%                     
%                     
%                     if violet_channel == 1
%                         violet = [violet part4(vindex(n)).Intensity_One];
%                     elseif violet_channel == 2
%                         violet = [violet part4(vindex(n)).Intensity_Two];
%                     elseif violet_channel == 3
%                         violet = [violet part4(vindex(n)).Intensity_Three];
%                     elseif violet_channel == 4
%                         violet = [violet part4(vindex(n)).Intensity_Four];
%                     end
%                 end
                violet(violet<0)=0;
            end
        end
    end
end

clearvars part4 part4_V_Normalized
med_violet = median(violet);
 
for s_l = 1:length(strain1)
    strain_i = strcat(['MR',num2str(strain1(s_l))]);
    dateDir = strcat(['/Users/reyer/Data/SingleCellEpi/',strain_i,'/',Date1{s_l}]);
    
    for j = 1:length(time)
        t = strcat(['t',num2str(time(j))]);
        close all
        
        timeDir = strcat([dateDir,'/',t]);
        d2 = dir([timeDir, '/*.mat']);
        %samples = length(d2);
        samples = sample{s_l}{j};
        file = strcat([timeDir,'/']);
        
        
        for k = samples
            clearvars part4_V_Normalized
            field1 = 'Cell';
            field2 = 'Volume';
            field3 = 'Center';
            field4 = 'Intensity_One';
            field5 = 'Intensity_Two';
            field6 = 'Intensity_Three';
            field7 = 'Intensity_Four';
            field8 = 'Original_One';
            field9 = 'Original_Two';
            field10 = 'Original_Three';
            field11 = 'Original_Four';
            
            varlist1 = who;
            
            part4_V_Normalized = struct(field1,[],field2,[],field3,[],field4,[],field5,[],field6,[],field7,[],field8,[],field9,[],field10,[],field11,[]);
            
            s1 = strcat([timeDir,'/sample_00',num2str(k),'.mat']);
            varlist = who;
            %varlist1 = varlist;
            varlist =strjoin(varlist','$|');
            load(s1,'-regexp', ['^(?!' varlist ')\w']);
            violet_norm = [];
            
            
            part4V_length = length(part4);
            part4_V_Normalized = struct(field1,[],field2,[],field3,[],field4,[],field5,[],field6,[],field7,[],field8,[],field9,[],field10,[],field11,[]);
            
            if discard_violet == 1
                violetvolume = zeros(part4V_length,1);
                for ni = 1:part4V_length
                    if violet_channel == 1
                        violetvolume(ni) = part4(ni).Intensity_One/(5*part4(ni).Volume);
                    elseif violet_channel == 2
                        violetvolume(ni) = part4(ni).Intensity_Two/(5*part4(ni).Volume);
                    elseif violet_channel == 3
                        violetvolume(ni) = part4(ni).Intensity_Three/(5*part4(ni).Volume);
                    elseif violet_channel == 4
                        violetvolume(ni) = part4(ni).Intensity_Four/(5*part4(ni).Volume);
                    end
                end
                
                [violetvolume,vindex] = sort(violetvolume);
                
                vindex(violetvolume<violet_back) = [];
                violetvolume(violetvolume<violet_back) = [];
                
%                 [violetvolume,vindex] = sort(violetvolume);
%                 v_cut = int32(violet_cut*part4V_length);
%                 violetvolume(1:v_cut) = [];
%                 vindex(1:v_cut) = [];
                
                
                
                
                for n = 1:length(vindex)
                    part4_V_Normalized(n).Cell = vindex(n);
                    part4_V_Normalized(n).Volume = part4(vindex(n)).Volume;
                    part4_V_Normalized(n).Center = part4(vindex(n)).Center;
                    
                    if violet_channel == 1
                        violet_norm(n) = med_violet/part4(vindex(n)).Intensity_One;
                        part4_V_Normalized(n).Intensity_One = violet_norm(n)*part4(vindex(n)).Intensity_One;
                        part4_V_Normalized(n).Original_One = part4(vindex(n)).Intensity_One;
                    elseif violet_channel == 2
                        violet_norm(n) = med_violet/part4(vindex(n)).Intensity_Two;
                        part4_V_Normalized(n).Intensity_Two = violet_norm(n)*part4(vindex(n)).Intensity_Two;
                        part4_V_Normalized(n).Original_Two = part4(vindex(n)).Intensity_Two;
                    elseif violet_channel == 3
                        violet_norm(n) = med_violet/part4(vindex(n)).Intensity_Three;
                        part4_V_Normalized(n).Intensity_Three = violet_norm(n)*part4(vindex(n)).Intensity_Three;
                        part4_V_Normalized(n).Original_Three = part4(vindex(n)).Intensity_Three;
                    elseif violet_channel == 4
                        violet_norm(n) = med_violet/part4(vindex(n)).Intensity_Four;
                        part4_V_Normalized(n).Intensity_Four = violet_norm(n)*part4(vindex(n)).Intensity_Four;
                        part4_V_Normalized(n).Original_Four = part4(vindex(n)).Intensity_Four;
                    end
                    
                    if blue_channel == 1
                        part4_V_Normalized(n).Intensity_One = part4(vindex(n)).Intensity_One;
                        part4_V_Normalized(n).Original_One = part4(vindex(n)).Intensity_One;
                    elseif blue_channel == 2
                        part4_V_Normalized(n).Intensity_Two = part4(vindex(n)).Intensity_Two;
                        part4_V_Normalized(n).Original_Two = part4(vindex(n)).Intensity_Two;
                    elseif blue_channel == 3
                        part4_V_Normalized(n).Intensity_Three = part4(vindex(n)).Intensity_Three;
                        part4_V_Normalized(n).Original_Three = part4(vindex(n)).Intensity_Three;
                    elseif blue_channel == 4
                        part4_V_Normalized(n).Intensity_Four = part4(vindex(n)).Intensity_Four;
                        part4_V_Normalized(n).Original_Four = part4(vindex(n)).Intensity_Four;
                    end
                    
                    if green_channel == 1
                        part4_V_Normalized(n).Intensity_One = part4(vindex(n)).Intensity_One*violet_norm(n);
                        part4_V_Normalized(n).Original_One = part4(vindex(n)).Intensity_One;
                    elseif green_channel == 2
                        part4_V_Normalized(n).Intensity_Two = part4(vindex(n)).Intensity_Two*violet_norm(n);
                        part4_V_Normalized(n).Original_Two = part4(vindex(n)).Intensity_Two;
                    elseif green_channel == 3
                        part4_V_Normalized(n).Intensity_Three = part4(vindex(n)).Intensity_Three*violet_norm(n);
                        part4_V_Normalized(n).Original_Three = part4(vindex(n)).Intensity_Three;
                    elseif green_channel == 4
                        part4_V_Normalized(n).Intensity_Four = part4(vindex(n)).Intensity_Four*violet_norm(n);
                        part4_V_Normalized(n).Original_Four = part4(vindex(n)).Intensity_Four;
                    end
                    
                    if red_channel == 1
                        part4_V_Normalized(n).Intensity_One = part4(vindex(n)).Intensity_One*violet_norm(n);
                        part4_V_Normalized(n).Original_One = part4(vindex(n)).Intensity_One;
                    elseif red_channel == 2
                        part4_V_Normalized(n).Intensity_Two = part4(vindex(n)).Intensity_Two*violet_norm(n);
                        part4_V_Normalized(n).Original_Two = part4(vindex(n)).Intensity_Two;
                    elseif red_channel == 3
                        part4_V_Normalized(n).Intensity_Three = part4(vindex(n)).Intensity_Three*violet_norm(n);
                        part4_V_Normalized(n).Original_Three = part4(vindex(n)).Intensity_Three;
                    elseif red_channel == 4
                        part4_V_Normalized(n).Intensity_Four = part4(vindex(n)).Intensity_Four*violet_norm(n);
                        part4_V_Normalized(n).Original_Four = part4(vindex(n)).Intensity_Four;
                    end
                    
                end
                
            else
                for n = 1:part4V_length
                    part4_V_Normalized(n).Cell = n;
                    part4_V_Normalized(n).Volume = part4(n).Volume;
                    part4_V_Normalized(n).Center = part4(n).Center;
                    
                    if violet_channel == 1
                        violet_norm(n) = med_violet/part4(n).Intensity_One;
                        part4_V_Normalized(n).Intensity_One = violet_norm(n)*part4(n).Intensity_One;
                    elseif violet_channel == 2
                        violet_norm(n) = med_violet/part4(n).Intensity_Two;
                        part4_V_Normalized(n).Intensity_Two = violet_norm(n)*part4(n).Intensity_Two;
                    elseif violet_channel == 3
                        violet_norm(n) = med_violet/part4(n).Intensity_Three;
                        part4_V_Normalized(n).Intensity_Three = violet_norm(n)*part4(n).Intensity_Three;
                    elseif violet_channel == 4
                        violet_norm(n) = med_violet/part4(n).Intensity_Four;
                        part4_V_Normalized(n).Intensity_Four = violet_norm(n)*part4(n).Intensity_Four;
                    end
                    
                    if blue_channel == 1
                        part4_V_Normalized(n).Intensity_One = part4(n).Intensity_One;
                    elseif blue_channel == 2
                        part4_V_Normalized(n).Intensity_Two = part4(n).Intensity_Two;
                    elseif blue_channel == 3
                        part4_V_Normalized(n).Intensity_Three = part4(n).Intensity_Three;
                    elseif blue_channel == 4
                        part4_V_Normalized(n).Intensity_Four = part4(n).Intensity_Four;
                    end
                    
                    if green_channel == 1
                        part4_V_Normalized(n).Intensity_One = part4(n).Intensity_One*violet_norm(n);
                    elseif green_channel == 2
                        part4_V_Normalized(n).Intensity_Two = part4(n).Intensity_Two*violet_norm(n);
                    elseif green_channel == 3
                        part4_V_Normalized(n).Intensity_Three = part4(n).Intensity_Three*violet_norm(n);
                    elseif green_channel == 4
                        part4_V_Normalized(n).Intensity_Four = part4(n).Intensity_Four*violet_norm(n);
                    end
                    
                    if red_channel == 1
                        part4_V_Normalized(n).Intensity_One = part4(n).Intensity_One*violet_norm(n);
                    elseif red_channel == 2
                        part4_V_Normalized(n).Intensity_Two = part4(n).Intensity_Two*violet_norm(n);
                    elseif red_channel == 3
                        part4_V_Normalized(n).Intensity_Three = part4(n).Intensity_Three*violet_norm(n);
                    elseif red_channel == 4
                        part4_V_Normalized(n).Intensity_Four = part4(n).Intensity_Four*violet_norm(n);
                    end
                    
                    
                    
                    
                    
                end
            end
            
            
            %save(filename, '-regexp', '^(?!(varlist|varlist1)$).')
            save(s1,'part4_V_Normalized','-append')
            %save(s, '-regexp', '^(?!(varlist1)$).')
            clearvars('-except',varlist1{:})
        end
        
        
    end
end