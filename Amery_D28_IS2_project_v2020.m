% Code to analyze ICESat-2 data (ATL06 and ATL03 data products)
% Author: Catherine C. Walker
% Accompaniment to Walker/Becker/Fricker paper submitted to GRL,
% October/revised December 2020
% Uses ICESat-2 ATL06 and ATL03 granules, downloaded from NSIDC, to
% determine ice surface topography and change over the Amery Ice Shelf
% ice calving front before and after the calving of D-28 in September 2019


%%
foldername = '/../ATL06/'; %directory where ATL06 files are stored
addpath(foldername) % add to path
files=dir('*.h5'); % determine number of files to loop through

%%

% TO EXPORT TO GIS SOFTWARE FOR MAPPING/FIGURES
% Loop to read in each granule to get lat/lon/elevation to write to shapefile 

for j = 1:size(files)
    clear D6 % clear current data
    

    try
% Get the filename
ATL06filename=files(j).name;
% Get rid of the extension
filename_1=ATL06filename(1:end-3);
% Read in the data (using subroutine read_ATL06_alt)
D6=read_ATL06_alt(ATL06filename);

for beam=1:3
    clear S_ATL06 % clear shapefile of previous data 

    for pair=1:2

        lats_to_write = D6(beam).latitude(:,pair);%
        lons_to_write = D6(beam).longitude(:,pair);%
        heights_to_write = D6(beam).h_li(:,pair);%

        lats_to_write_a = lats_to_write(1:10:end);
        lons_to_write_a = lons_to_write(1:10:end);
        heights_to_write_a = heights_to_write(1:10:end);

        for i = 1:length(lats_to_write_a)
            S_ATL06(i).Geometry = 'point';
            S_ATL06(i).Lat = lats_to_write_a(i);
            S_ATL06(i).Lon = lons_to_write_a(i);
            S_ATL06(i).Elev = heights_to_write_a(i);
        end
        filename_1 = [ATL06filename(1:end-3) '_b' num2str(beam) '_p' num2str(pair)];
        shapewrite(S_ATL06, filename_1);
    end
end

    catch

    disp(j) % these didn't work
    end
end


%%
% READ IN THE DATA and place into structure "Amery"
Amery=[];
for i = 1:length(files)
    
ATL06_file1=[foldername files(i).name];
callName = files(i).name(11:end-10); % get the granule name

try
D6=read_ATL06_alt(ATL06_file1);


for b = 1:3
beamNum = ['B' num2str(b)];

for p = 1:2
    % place data into structure by beam (1-3) and pair (1-2) number
    Amery.(callName).(beamNum).x_atc(:,p) = D6(b).x_atc(:,p);
    Amery.(callName).(beamNum).h_li(:,p) = D6(b).h_li(:,p);
    Amery.(callName).(beamNum).latitude(:,p) = D6(b).latitude(:,p);
    Amery.(callName).(beamNum).longitude(:,p) = D6(b).longitude(:,p);
end

end
catch
    disp(i)
end
clear D6
end



%%
% De-trend background 
for i = 1:length(files)
    
ATL06_file1=[foldername files(i).name];
callName = files(i).name(11:end-10); % get the granule name

try

for b = 1:3
beamNum = ['B' num2str(b)];

for p = 1:2
    
    Amery.(callName).(beamNum).x_atc(:,p) = D6(b).x_atc(:,p);
    Amery.(callName).(beamNum).h_li(:,p) = D6(b).h_li(:,p);

    
    h_li = Amery.(callName).(beamNum).h_li(:,p);
    h_li(h_li>100)=NaN;h_li(h_li<-30)=NaN;
    x_atc = Amery.(callName).(beamNum).x_atc(:,p);
    x_atc(h_li>100)=NaN;x_atc(h_li<-30)=NaN;

    
    Good = isnan(x_atc) + isnan(h_li);% 
    f2 = fit(h_li(Good == 0),x_atc(Good == 0),'smoothingspline','SmoothingParam',0.1);
    Amery.(callName).(beamNum).h_li_fit(:,p) = x_atc(:,p)-f2(h_li(:,p));
    
    clear h_li
    clear x_atc
    clear latitude
    clear longitude
    
    clear fit_h_li
    clear fit_x_atc
    clear fit_latitude
    clear fit_longitude

end

end
catch
    disp(i) % these didn't work
end
clear D6
end



%%
% PLOTS

%example, all beams along RGT 0104 on April 05, 2019 (cycle 3):
granuleCall='ATL06_20190405024508_01040312';
figure; 
scatter3(Amery.(granuleCall).B1.latitude_x,-Amery.(granuleCall).B1.longitude_x,Amery.(granuleCall).B1.h_li_x,20,Amery.(granuleCall).B1.h_li_x,'filled');hold on;
hold on;
scatter3(Amery.(granuleCall).B2.latitude_x,-Amery.(granuleCall).B2.longitude_x,Amery.(granuleCall).B2.h_li_x,20,Amery.(granuleCall).B2.h_li_x,'filled');hold on;
hold on;
scatter3(Amery.(granuleCall).B3.latitude_x,-Amery.(granuleCall).B3.longitude_x,Amery.(granuleCall).B3.h_li_x,20,Amery.(granuleCall).B3.h_li_x,'filled');hold on;

% example, plot all the lines over Amery
fn=fieldnames(Amery);
figure, 
for j = 1:length(fn)
    
    dateCall = char(fn(j));
    
    for b=1:3
        beamNum = ['B' num2str(b)];
        try
        for p = 1:2
            scatter3(Amery.(dateCall).(beamNum).longitudeN(:,p), Amery.(dateCall).(beamNum).latitudeN(:,p), Amery.(dateCall).(beamNum).h_liN(:,p),5,Amery.(dateCall).(beamNum).h_liN(:,p));hold on;
        end
        catch
            disp([dateCall '_' beamNum])
        end
        
    end
end


%% 
% ATL06 reader originally written by Ben Smith 2018
% Github source: https://github.com/SmithB/ICESat2/blob/6ca3e5bd7a23a6fc942f838fab1ee964d05daf17/matlab/readATL06.m

function D3=read_ATL06(filename, pair, fields)

if ~exist('pair', 'var')
    pair=[1,2,3];
end

if ~exist('fields','var')
    for k=1:length(pair)
        try
            temp=h5info(filename, sprintf('/gt%dl/land_ice_segments', pair(k)));
            fields={temp.Datasets.Name};
        end
        if exist('fields','var')
            break
        end
    end
    fields{end+1}='/ground_track/x_atc';
    fields{end+1}='/ground_track/y_atc';   
end

fields=fields(~ismember(fields,'segment_id'));


beams={'l','r'};
for kP=1:length(pair)
    
    for kB=1:2
        ID{kB}=h5read(filename,sprintf('/gt%d%s/land_ice_segments/%s', pair(kP), beams{kB}, 'segment_id')); 
    end
    ID_both=unique(cat(1, ID{:}));
    
    for kB=1:2
        [~, in_ind{kB}, out_ind{kB}]=intersect(ID{kB}, ID_both);
    end
    
    for kF=1:length(fields)
        temp=cell(1,2);
        for kB=1:2
            temp{kB}=h5read(filename,sprintf('/gt%d%s/land_ice_segments/%s', pair(kP), beams{kB}, fields{kF})); 
        end
        D3_name=strrep(fields{kF},'/ground_track/','');
        D3(kP).(D3_name)=NaN(numel(ID_both), 2);
        
        for kB=1:2
            D3(kP).(D3_name)(out_ind{kB}, kB)=temp{kB}(in_ind{kB});
        end
    end
end


fields=fieldnames(D3);
for kP=1:3
    bad=D3(kP).h_li>1e30;
    for kf=1:length(fields);
        if isa(D3(kP).(fields{kf}),'single') ||  isa(D3(kP).(fields{kf}),'double')
            D3(kP).(fields{kf})(bad)=NaN;        
        end
    end
end
end



%% Front advance calculations

% distances measured along Vw and Ve in GIS (Fig. 1a)
west_front_L = [2977 2250 2542 1742 1379 2395 2467];
east_front_L = [2414 2517 2157 1497 1027 2106 1646];

% dates used 
dates_spots=[datenum(2008,1,2) datenum(2010,3,28) datenum(2012,2,7) datenum(2013,11,7) datenum(2014,12,19) datenum(2016,3,20) datenum(2018,2,8) datenum(2019,9,12)];

% figure out date change
for i=1:length(dates_spots)-1
dates_change(i) = dates_spots(i+1)-dates_spots(i);
end

% determine m/day change
for i=1:length(west_front_L)
west_front_rate(i) = west_front_L(i)./dates_change(i);
east_front_rate(i) = east_front_L(i)./dates_change(i);
end



