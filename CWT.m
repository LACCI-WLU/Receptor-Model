clear all
clc

% This program contains the core function of CWT methond.
% If the Mapping Toolbox was not installed, the M_Map can be a replacement. Below are links: 

% https://www.eoas.ubc.ca/~rich/mapug.html#p2
% http://www.rsmas.miami.edu/personal/miskandarani/MatlabScripts/M_map/private/mapug.html


% Copyright (C) Aug 2017 - present Xianda Gong, Leibniz Institute for Tropospheric Research (TROPOS)


%% read in the trajectories:
sourcedirectory=['C:\LAB\Data\MarParCloud\BT\BT_1\Sample1\'];
raw=dir([sourcedirectory,'*.txt']);
raw=sortObj(raw); % sort raw data according to its name


for i=1:length(raw)
    [versn, name, ext] = fileparts(raw(i).name);
    fid=fopen([sourcedirectory,name,'.txt']);
    data1=textscan(fid,'%s', 'delimiter', '\n');        %reads each line to find "1 PRESSURE" 
    data2=data1{1};
    for j=1:20   
        data3=cell2mat(data2(j,:));
        k=strcmp('7 PRESSURE THETA    AIR_TEMP RAINFALL MIXDEPTH RELHUMID SUN_FLUX',data3);
        if k == 1, break, end
    end
    fclose(fid);
    fid=fopen([sourcedirectory,name,'.txt']);
    data4=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter', '\t','HeaderLines',j);      %reads the data
    fclose(fid);    
    data5=cell2mat(data4);
    if size(data5,1)==723
        data5_h1=data5(1:3:end,:); %select useful data
        data5_h2=data5(2:3:end,:);
        data5_h3=data5(3:3:end,:);    
        CVAOtime(i,:)=datenum(num2str(data5(1,3:6)),'yy mm dd HH');    
        CVAOlat_h1(:,i)=data5_h1(:,10);
        CVAOlong_h1(:,i)=data5_h1(:,11);
        CVAOheight_h1(:,i)=data5_h1(:,12);
        CVAOprecip_h1(:,i)=data5_h1(:,16);
        CVAOlat_h2(:,i)=data5_h2(:,10);
        CVAOlong_h2(:,i)=data5_h2(:,11);
        CVAOprecip_h2(:,i)=data5_h2(:,16);
        CVAOheight_h2(:,i)=data5_h2(:,12);
        CVAOlat_h3(:,i)=data5_h3(:,10);
        CVAOlong_h3(:,i)=data5_h3(:,11); 
        CVAOprecip_h3(:,i)=data5_h3(:,16);
        CVAOheight_h3(:,i)=data5_h3(:,12);
    end
end
clear ans data1 data2 data2 data3 data4 data5 data5_h1 data5_h2 data5_h3 i j k name raw sourcedirectory versn ext fid;



%% load the PM10, dN20avg and CCN data, or what ever data you want; neeed time and conc
% The CWT need to combine the backtrajection with a certain parameter.
% This parameter can be total particle number concentration, PM10, PM1, CCN number concentration and so on.

load('C:\LAB\Data\MarParCloud\CVAO_aerosol\PNSD\Step1_3_Combine_analysis.mat','tiavg','dN500lavg','dN500l','Ti');
load('C:\LAB\Data\MarParCloud\CVAO_aerosol\CCN\Step2_3_Further_analysis_APS.mat','CCN_conc_015','CCN_conc_02','CCN_conc_03','CCN_conc_05','CCN_conc_07','time_CCN_015','time_CCN_02','time_CCN_03','time_CCN_05','time_CCN_07');

%% match the trajectories with the concentration time wise:
CWT_value=CCN_conc_03;
CWT_time=time_CCN_03;
time_all=[CVAOtime;MVtime];
long_all=[CVAOlong_h2(1:145,:) MVlong_h2(1:145,:)];
lat_all=[CVAOlat_h2(1:145,:) MVlat_h2(1:145,:)];
precip_all=[CVAOprecip_h2(1:145,:) MVprecip_h2(1:145,:)];
height_all=[CVAOheight_h2(1:145,:) MVheight_h2(1:145,:)];
CWT_match=NaN(size(time_all));

CWT_value=CCN_conc_03;
CWT_time=time_CCN_03;
time_all=[CVAOtime];
long_all=[CVAOlong_h2(1:145,:)];
lat_all=[CVAOlat_h2(1:145,:)];
precip_all=[CVAOprecip_h2(1:145,:)];
height_all=[CVAOheight_h2(1:145,:)];
CWT_match=NaN(size(time_all));


for i=1:length(time_all)
    [rm(i,1) cm(i,1)]=min(abs(CWT_time-time_all(i)));
    if rm(i)<1/24 %& cm(i)>3 & cm(i)<length(PSCF_value)-3
        CWT_match(i,1)=CWT_value(cm(i));
    end
end

[row, col] = find(isnan(CWT_match));
CWT_match(row)=[];
time_all(row)=[];
long_all(:,row)=[];
lat_all(:,row)=[];
height_all(:,row)=[];
precip_all(:,row)=[];


%% Precipitation and height filter:
precip_gesamt=sum(precip_all);
ind=zeros(size(precip_all,1),size(precip_all,2));
for i=1:size(precip_all,2)
    for j=1:size(precip_all,1)
        if precip_all(j,i)>0.5 && precip_all(j)>10% || height_all(j,i)<10 || height_all(j,i)>2000       % here you can adjust the threshold values for precipitation and height
            ind(j:end,i)=99999;
            lat_all(j:end,i)=NaN;
            long_all(j:end,i)=NaN;
        end
    end
end


% calculation
long_start=-80;
long_end=40;
lat_start=-10;
lat_end=80;
p=5;   
q=5;

Vm=[];
id1=0;
id2=0;

CWTB=zeros((long_end-long_start)/p,(lat_end-lat_start)/q);

for i=long_start:p:long_end-p
    id1=id1+1;
    for j=lat_start:q:lat_end-q
        id2=id2+1;
        for k=1:length(CWT_match) % present the BT
            m=0;
             for w=1:size(long_all) % present each point
                 if long_all(w,k)>=i & long_all(w,k)<i+p & lat_all(w,k)>=j & lat_all(w,k)<j+q
                     m=m+1;
                 end
             end
             Vm=[Vm m];
        end
        if sum(Vm)>=10
            CWTB(id1,id2)=Vm*CWT_match/sum(Vm);
        else
            CWTB(id1,id2)=NaN;
        end
        
        Vm=[];
    end
    fprintf('%s %f %s\n','Completed', (id1)/(long_end-long_start)*p*100,'%');
    id2=0;
end



                 
long_grid=NaN(size(CWTB,1),size(CWTB,2));
m=0;
for i=long_start:p:long_end-p
    m=m+1;
    long_grid(m,:)=i;
end
lat_grid=NaN(size(CWTB,1),size(CWTB,2));
m=0;
for i=lat_start:q:lat_end-q
    m=m+1;
    lat_grid(:,m)=i;
end
clear i j w z m n id;
% plot the result
% combine

% figure('Color','white');
% worldmap([-10 80],[-90 40]);
% load coast;  
% pcolorm(lat_grid,long_grid,CWTB);
% plotm(lat, long,'black','LineWidth',2);
% hold on
% c=colorbar('eastoutside');
% ylabel(c,'CWT');
% text(0.25,2.5,'(b)','color','black','FontSize',14);

figure
m_proj('miller','lat',[-10 80],'long',[-90 40]);
m_coast('linewidth',2,'color','k');
m_coast('patch',[.6 .6 .6]);
m_grid('linewi',2,'tickdir','in','XaxisLocation','bottom','backcolor',[1 1 1]);
hold on
m_pcolor(long_grid,lat_grid,CWTB)%;shading flat%;colormap(map);%;shading interp;                                  
shading interp
m_gshhs_i('color','k');
hold on
m_coast('linewidth',2,'color','k');
c=colorbar;
ylabel(c,'CCN')
hold on
m_line(23.35,-71.95,'marker','square','markersize',6,'color','r');
title('CCN03')
break