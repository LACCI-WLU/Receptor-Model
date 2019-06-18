clear all
close all


% This program contains the core function of PSCF methond.
% If the Mapping Toolbox was not installed, the M_Map can be a replacement. Below are links: 

% https://www.eoas.ubc.ca/~rich/mapug.html#p2
% http://www.rsmas.miami.edu/personal/miskandarani/MatlabScripts/M_map/private/mapug.html


% Copyright (C) Aug 2017 - present Xianda Gong, Leibniz Institute for Tropospheric Research (TROPOS)

%% read in the trajectories:
sourcedirectory=['C:\BT\'];
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
    data5_h1=data5(1:3:end,:);
    data5_h2=data5(2:3:end,:);
    data5_h3=data5(3:3:end,:);
    
    time(i,:)=datenum(num2str(data5(1,3:6)),'yy mm dd HH');  
    
    lat_h1(:,i)=data5_h1(:,10);
    long_h1(:,i)=data5_h1(:,11);
    precip_h1(:,i)=data5_h1(:,16);
    height_h1(:,i)=data5_h1(:,12);
    
    lat_h2(:,i)=data5_h2(:,10);
    long_h2(:,i)=data5_h2(:,11);
    precip_h2(:,i)=data5_h2(:,16);
    height_h2(:,i)=data5_h2(:,12);
    
    lat_h3(:,i)=data5_h3(:,10);
    long_h3(:,i)=data5_h3(:,11);
    precip_h3(:,i)=data5_h3(:,16);
    height_h3(:,i)=data5_h3(:,12);
    fprintf('%s %f %s\n','Completed', i/(length(raw))*100,'%');
    
end
% final data:
time_all=[time;time;time];
lat_all=[lat_h1 lat_h2 lat_h3];
long_all=[long_h1 long_h2 long_h3];
precip_all=[precip_h1 precip_h2 precip_h3];
height_all=[height_h1 height_h2 height_h3];   



%% Trim backtrajectory:
ln=6;    %New backtrajectory in days
lat_all(ln*24+1:end,:)=[];
long_all(ln*24+1:end,:)=[];



%% Plot the trajectories

figure
m_proj('stereographic','latitude',90,'radius',80,'rotangle',0);
%m_proj('miller','lat',[0 90],'long',[-60 60]);
m_coast('linewidth',2,'color','k');
m_grid('linewi',2,'tickdir','in','XaxisLocation','top');
m_coast('patch',[.85 .77 .65]);
m_grid('linewi',2,'backcolor',[0.9 0.99 1]);
hold on
hc=m_plot(long_h1,lat_h1,'b','linewidth',1);

clear ans data1 data2 data3 data5_h1 data5_h2 data5_h3 ext fid i j k name raw sourcedirectory versn;

%% load the PM10, dN20avg and CCN data
% The PSCF need to combine the backtrajection with a certain parameter.
% This parameter can be total particle number concentration, PM10, PM1, CCN number concentration and so on.
load('C:\LAB\Data\Cyprus_Campaign\Classification\Classification_version_04.mat','CCN_conc_02','CCN_conc_03','CCN_conc_05','CCN_conc_07','time_CCN_02','time_CCN_03','time_CCN_05','time_CCN_07','Time1','dN20avg','dN70avg','dN80avg','dN200avg','PM10');

%% get the background 
BStart=[datenum('04-Apr-2017 21:00:00');datenum('05-Apr-2017 06:00:00');datenum('05-Apr-2017 20:00:00');datenum('06-Apr-2017 06:00:00');datenum('09-Apr-2017 22:00:00');datenum('10-Apr-2017 20:00:00');datenum('11-Apr-2017 21:00:00');datenum('12-Apr-2017 05:00:00');datenum('12-Apr-2017 21:00:00');datenum('14-Apr-2017 23:00:00');datenum('15-Apr-2017 02:00:00');datenum('15-Apr-2017 04:00:00');datenum('15-Apr-2017 22:00:00');datenum('16-Apr-2017 05:00:00');datenum('16-Apr-2017 21:00:00');datenum('17-Apr-2017 05:00:00');datenum('17-Apr-2017 21:00:00');datenum('18-Apr-2017 21:00:00');datenum('19-Apr-2017 05:00:00');datenum('19-Apr-2017 21:00:00');datenum('20-Apr-2017 05:00:00');datenum('20-Apr-2017 22:00:00');datenum('21-Apr-2017 20:00:00');datenum('22-Apr-2017 04:00:00');datenum('24-Apr-2017 22:00:00');datenum('25-Apr-2017 05:00:00');datenum('25-Apr-2017 20:00:00');datenum('26-Apr-2017 05:00:00');datenum('26-Apr-2017 21:00:00');datenum('27-Apr-2017 05:00:00');datenum('27-Apr-2017 21:00:00')];
BEnd=[datenum('05-Apr-2017 03:00:00');datenum('05-Apr-2017 17:00:00');datenum('06-Apr-2017 04:00:00');datenum('09-Apr-2017 17:00:00');datenum('10-Apr-2017 16:00:00');datenum('11-Apr-2017 17:00:00');datenum('12-Apr-2017 03:00:00');datenum('12-Apr-2017 13:00:00');datenum('14-Apr-2017 22:00:00');datenum('15-Apr-2017 01:00:00');datenum('15-Apr-2017 03:00:00');datenum('15-Apr-2017 18:00:00');datenum('16-Apr-2017 04:00:00');datenum('16-Apr-2017 18:00:00');datenum('17-Apr-2017 03:00:00');datenum('17-Apr-2017 17:00:00');datenum('18-Apr-2017 20:00:00');datenum('19-Apr-2017 03:00:00');datenum('19-Apr-2017 18:00:00');datenum('20-Apr-2017 04:00:00');datenum('20-Apr-2017 18:00:00');datenum('21-Apr-2017 18:00:00');datenum('22-Apr-2017 03:00:00');datenum('24-Apr-2017 18:00:00');datenum('25-Apr-2017 03:00:00');datenum('25-Apr-2017 17:00:00');datenum('26-Apr-2017 03:00:00');datenum('26-Apr-2017 16:00:00');datenum('27-Apr-2017 03:00:00');datenum('27-Apr-2017 17:00:00');datenum('28-Apr-2017 01:00:00')];
id=[];
for i=1:length(BStart)
    id=[id;find(Time1>datenum('04-Apr-2017 17:00:00') & Time1<datenum('23-Apr-2017 00:00:00') & Time1>=BStart(i) & Time1<BEnd(i))];
end


%% match the trajectories with the concentration time wise:
PSCF_value=PM10(id);
PSCF_time=Time1(id);
PSCF_match=NaN(size(time_all));

for i=1:length(time_all)
    [rm(i,1) cm(i,1)]=min(abs(PSCF_time-time_all(i)));
    if rm(i)<1/24 %& cm(i)>3 & cm(i)<length(PSCF_value)-3
        PSCF_match(i,1)=PSCF_value(cm(i));
    end
end
PSCF_p=prctile(PSCF_match,75);
[row, col] = find(isnan(PSCF_match));
PSCF_match(row)=[];
time_all(row)=[];
long_all(:,row)=[];
lat_all(:,row)=[];


%% Precipitation and height filter:
precip_gesamt=sum(precip_all);
ind=zeros(size(precip_all,1),size(precip_all,2));
for i=1:size(precip_all,2)
    for j=1:size(precip_all,1)
        if precip_all(j,i)>0.1 && precip_gesamt(j)>5 || height_all(j,i)<10 || height_all(j,i)>1000       % here you can adjust the threshold values for precipitation and height
            ind(j:end,i)=99999;
            lat_all(j:end,i)=NaN;
            long_all(j:end,i)=NaN;
        end
    end
end



%% Create analysis grid:
% This part help you difine the grid. p and q are unit of grid. 5 is recommended to use. 
long_start=-180;    
long_end=180;
lat_start=0;
lat_end=90;
p=5;   
q=5;
m=0;    
n=0;

PSCF_m=NaN((long_end-long_start)/p,(lat_end-lat_start)/q);
PSCF_n=NaN((long_end-long_start)/p,(lat_end-lat_start)/q);
PSCF=NaN((long_end-long_start)/p,(lat_end-lat_start)/q);
long_grid=NaN((long_end-long_start)/p,(lat_end-lat_start)/q);
lat_grid=NaN((long_end-long_start)/p,(lat_end-lat_start)/q);
precip_grid=NaN((long_end-long_start)/p+1,(lat_end-lat_start)/q+1);
precip_p75=NaN((long_end-long_start)/p+1,(lat_end-lat_start)/q+1);


for i=long_start:p:long_end-p
    m=m+1;
    for j=lat_start:q:lat_end-q
        n=n+1;
        long_gitter=[i;i;i+p;i+p];
        lat_gitter=[j;j+q;j+q;j];
        [am bm]=find(long_all>i & long_all<i+p & lat_all>j & lat_all<j+q); 
        long_grid(m,n)=i+p/2;
        lat_grid(m,n)=j+q/2;
        PSCF_m(m,n)=length(bm);
        [cm dm]=find(PSCF_match(bm)>PSCF_p);
        PSCF_n(m,n)=length(dm);
        PSCF(m,n)=PSCF_n(m,n)/PSCF_m(m,n);
    end
    n=0;
end



%% Weighting function method one

log_PSCF_m=log10(PSCF_m+1);
log_PSCF_m_max=log10(max(PSCF_m(:)));
ratio=log_PSCF_m./log_PSCF_m_max;

%discrete:
w=NaN(size(ratio,1),size(ratio,2));
for i=1:size(ratio,1)
    for j=1:size(ratio,2)
        if ratio(i,j)>0.85
            w(i,j)=1;
        end
        if ratio(i,j)>0.6 & ratio(i,j)<0.85
            w(i,j)=0.725;
            PSCF(i,j)=PSCF(i,j)*0.725;
            %precip_mean(i,j)=precip_mean(i,j)*0.725;
        end
        if ratio(i,j)>0.35 & ratio(i,j)<0.6
            w(i,j)=0.475;
            PSCF(i,j)=PSCF(i,j)*0.475;
            %precip_mean(i,j)=precip_mean(i,j)*0.475;
        end
        if ratio(i,j)<0.35 
            w(i,j)=0.175;
            PSCF(i,j)=PSCF(i,j)*0.175;
            %precip_mean(i,j)=precip_mean(i,j)*0.175;
        end
    end
end

% %% Weighting function method two
% for i=1:size(PSCF,1)
%     for j=1:size(PSCF,2)
%         if PSCF_m(i,j)<=10
%             PSCF(i,j)=NaN;
%         end
%         if PSCF_m(i,j)<30 & PSCF_m(i,j)>10
%                     PSCF(i,j)=PSCF(i,j)*0.7;
%         end
%         if PSCF_m(i,j)<=10 & PSCF_m(i,j)>5
%                     PSCF(i,j)=PSCF(i,j)*0.4;
%         end
%         if PSCF_m(i,j)<=5
%                     PSCF(i,j)=PSCF(i,j)*0.2;
%         end
%     end
% end


%% Edit the grid so that the plots show everything well:
for i=1:size(lat_grid,1)
    lat_grid1(i)=90;
    PSCF1(i)=mean(PSCF(:,end));
end
lat_grid=[lat_grid,lat_grid1'];
long_grid=[long_grid(:,1),long_grid];
% PSCF=[PSCF(:,1),PSCF];
PSCF=[PSCF,PSCF1'];
PSCF_m=[PSCF_m(:,1),PSCF_m];

for i=1:size(long_grid,2)
    long_grid1(i)=-180;
    long_grid2(i)=180;
end
long_grid=[long_grid1; long_grid; long_grid2];
lat_grid=[lat_grid(1,:); lat_grid; lat_grid(1,:)];
PSCF=[(PSCF(1,:)+PSCF(end,:))/2; PSCF; (PSCF(1,:)+PSCF(end,:))/2];
PSCF_m=[PSCF_m(1,:); PSCF_m; PSCF_m(1,:)];



%% plotten 

% figure
% m_proj('lambert','lon',[-40 -180],'lat',[45 90]);
% m_coast('linewidth',2,'color','k');
% m_coast('patch',[.6 .6 .6]);
% m_grid('linewi',2,'tickdir','in','XaxisLocation','bottom','backcolor',[1 1 1]);
% hold on
% m_pcolor(long_grid,lat_grid,PSCF);shading interp;                                 
% m_gshhs_i('color','k');;
% hold on
% %m_line(-148.72,70.3,'marker','square','markersize',6,'color','r')
% %m_coast('linewidth',2,'color','k');
% c=colorbar
% ylabel(c,'PSCF') 
% %m_grid('linewi',2,'tickdir','in','XaxisLocation','bottom','backcolor',[1 1 1]);

figure
%m_proj('lambert','lon',[180 -180],'lat',[-90 -20]);
m_proj('stereographic','latitude',90,'radius',80,'rotangle',0);
%m_proj('stereographic','longitude',-133,'latitude',90,'radius',50)%,'rotangle',45);   
%m_proj('albers equal-area','lat',[-90 -30],'long',[-180 180],'rect','on');
%m_proj('miller','lat',[-90 0]); 
m_coast('linewidth',2,'color','k');
m_coast('patch',[.6 .6 .6]);
m_grid('linewi',2,'tickdir','in','XaxisLocation','bottom','backcolor',[1 1 1]);
hold on
m_pcolor(long_grid,lat_grid,PSCF)%;shading flat%;colormap(map);%;shading interp;                            
m_gshhs_i('color','k');;
hold on
m_coast('linewidth',2,'color','k');
c=colorbar
ylabel(c,'PSCF')
hold on
m_line(23.35,-71.95,'marker','square','markersize',6,'color','r')
break


