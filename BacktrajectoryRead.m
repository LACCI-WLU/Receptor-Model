clear all
clc

% This program aim to read raw HYSPLIT back trajectory info.
% % Copyright (C) Aug 2017 - present Xianda Gong, Leibniz Institute for Tropospheric Research (TROPOS)

%% read in the trajectories:
sourcedirectory=['C:\LAB\Data\Cyprus_Campaign\Backtrajectory\Sample1\'];
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
    if size(data5,1)==482
        data5_h1=data5(1:2:end,:); %select useful data
        data5_h2=data5(2:2:end,:);    
        time(i,:)=datenum(num2str(data5(1,3:6)),'yy mm dd HH');    
        lat_h1(:,i)=data5_h1(:,10);
        long_h1(:,i)=data5_h1(:,11);
        height_h1(:,i)=data5_h1(:,12);
        precip_h1(:,i)=data5_h1(:,16);
        lat_h2(:,i)=data5_h2(:,10);
        long_h2(:,i)=data5_h2(:,11);
        precip_h2(:,i)=data5_h2(:,16);
        height_h2(:,i)=data5_h2(:,12);
    end
end
clear ans data1 data2 data2 data3 data4 data5_h1 data5_h2 i j k name raw sourcedirectory versn ext fid;


%% plot backtrajectory

worldmap([-10 90],[-120 60]);
load coast;  
plotm(lat, long);
land = shaperead('landareas', 'UseGeoCoords', true);  
geoshow(land, 'FaceColor', [0.2 0.8 0.2]);  
lakes = shaperead('worldlakes', 'UseGeoCoords', true);  
geoshow(lakes, 'FaceColor', 'blue');  
rivers = shaperead('worldrivers', 'UseGeoCoords', true);  
geoshow(rivers, 'Color', 'blue');
cities = shaperead('worldcities', 'UseGeoCoords', true);
geoshow(cities,'Marker','.', 'Color','red');

StartT=datenum('01-Apr-2017 00:00:00');
EndT=datenum('01-May-2017 00:00:00');
id=find(time>=StartT & time<=EndT);

for i=1:length(id)
    hold off
    worldmap([-10 90],[-120 60]);
    load coast;  
    plotm(lat, long);
%     plotm(lat_h1(1:241,id(i)),long_h1(1:241,id(i)),'black','LineWidth',1);
%     hold on
    plotm(lat_h2(1:241,id(i)),long_h2(1:241,id(i)),'red','LineWidth',1);
    legend(datestr(time(id(i))));
    pause;
end

clear i id;

