%Created by: Rachel Holser (rholser@ucsc.edu) and Arina Favilla
%(afavilla@ucsc.edu)
%Last Updated: 18-Sep-2022

%Completes calculations needed and generates Figures 4 and 5 of Holser et
%al. 2022

%Requires:
%           DayOfTripColormap.mat
%           findseq.m
%           (https://www.mathworks.com/matlabcentral/fileexchange/28113-findseq)

%% Load data from .mat files for PB2017 deployments

load('2017001_6018_RRHmod_TV3.mat')
A=load('2017002_U 20_RRHmod_TV4_alpha_2.mat');
B=load('2017003_5712_RRHmod_TV4_alpha_2.mat');
C=load('2017004_B794_RRHmod_TV4_alpha_2.mat');
D=load('2017005_6871_RRHmod_TV4_alpha_2.mat');
%E=load('2017006_5950_RRHmod_TV3.mat'); - no dive data
F=load('2017007_6767_RRHmod_TV4_alpha_2.mat'); 
G=load('2017008_7238_RRHmod_TV4_alpha_2.mat');
H=load('2017009_6260_RRHmod_TV4_alpha_2.mat');

%% Figure 5 - Kernel densities of dive statistics, calculated day and night

%Define grid over which to estimate density for Duration vs Depth
gridx1=0:10:1000; %depths
gridx2=0:2:40; %duration (minutes)
[X,Y] = meshgrid(gridx1, gridx2);
x1 = X(:);
x2 = Y(:);
xi = [x1 x2];

%Create data array for ksdensity calculation, using NIGHT dives only (based
%on solar elevation).  Data in array needs to match the xi dimensions
data=[A.DiveStat.Maxdepth(A.DiveStat.SolarEl<0) A.DiveStat.Dduration(A.DiveStat.SolarEl<0)/60];
data=[data;B.DiveStat.Maxdepth(B.DiveStat.SolarEl<0) B.DiveStat.Dduration(B.DiveStat.SolarEl<0)/60;...
    C.DiveStat.Maxdepth(C.DiveStat.SolarEl<0) C.DiveStat.Dduration(C.DiveStat.SolarEl<0)/60;...
    F.DiveStat.Maxdepth(F.DiveStat.SolarEl<0) F.DiveStat.Dduration(F.DiveStat.SolarEl<0)/60;...
    G.DiveStat.Maxdepth(G.DiveStat.SolarEl<0) G.DiveStat.Dduration(G.DiveStat.SolarEl<0)/60];

%Calculate kernel density of NIGHT duration against depth for all other 2017 animals
[f,xi2] = ksdensity(data,xi,'Function','pdf');

%Reshape ksdensity output for contour
count=1;
for i=1:size(Y,2)
    fn(:,i)=f(count:count+20);
    count=count+21;
end

%Create data array for ksdensity calculation, using DAYTIME dives only 
%(based on solar elevation).  Data in array needs to match the xi dimensions
data2=[A.DiveStat.Maxdepth(A.DiveStat.SolarEl>=0) A.DiveStat.Dduration(A.DiveStat.SolarEl>=0)/60];
data2=[data2;B.DiveStat.Maxdepth(B.DiveStat.SolarEl>=0) B.DiveStat.Dduration(B.DiveStat.SolarEl>=0)/60;...
    C.DiveStat.Maxdepth(C.DiveStat.SolarEl>=0) C.DiveStat.Dduration(C.DiveStat.SolarEl>=0)/60;...
    F.DiveStat.Maxdepth(F.DiveStat.SolarEl>=0) F.DiveStat.Dduration(F.DiveStat.SolarEl>=0)/60;...
    G.DiveStat.Maxdepth(G.DiveStat.SolarEl>=0) G.DiveStat.Dduration(G.DiveStat.SolarEl>=0)/60];

%Calculate kernel density of DAY duration against depth for all other 2017 animals
[f2,xi22] = ksdensity(data2,xi);

%Reshape ksdensity output for contour
count=1;
for i=1:size(Y,2)
    fd(:,i)=f2(count:count+20);
    count=count+21;
end

%Define grid over which to estimate kernel density for PDI~Duration
gridx1=0:0.1:100; %PDI
gridx2=0:2:40; %duration (minutes)
[X2,Y2] = meshgrid(gridx1, gridx2);
x1 = X2(:);
x2 = Y2(:);
xi = [x1 x2];

%pull out DAY PDI and duration data (columns need to match order in xi above)
dataDay=[A.DiveStat.PDI(A.DiveStat.SolarEl>=0)/60 A.DiveStat.Dduration(A.DiveStat.SolarEl>=0)/60];
dataDay=[dataDay;B.DiveStat.PDI(B.DiveStat.SolarEl>=0)/60 B.DiveStat.Dduration(B.DiveStat.SolarEl>=0)/60;...
    C.DiveStat.PDI(C.DiveStat.SolarEl>=0)/60 C.DiveStat.Dduration(C.DiveStat.SolarEl>=0)/60;...
    D.DiveStat.PDI(D.DiveStat.SolarEl>=0)/60 D.DiveStat.Dduration(D.DiveStat.SolarEl>=0)/60;...
    F.DiveStat.PDI(F.DiveStat.SolarEl>=0)/60 F.DiveStat.Dduration(F.DiveStat.SolarEl>=0)/60;...
    G.DiveStat.PDI(G.DiveStat.SolarEl>=0)/60 G.DiveStat.Dduration(G.DiveStat.SolarEl>=0)/60];

[fD,xi2] = ksdensity(dataDay,xi,'Function','pdf');

%reshape output array for use in contour plot
count=1;
for i=1:size(Y2,2)
    fDb(:,i)=fD(count:count+20);
    count=count+21;
end

%pull out NIGHT PDI and duration data (columns need to match order in xi above)
dataNight=[A.DiveStat.PDI(A.DiveStat.SolarEl<0)/60 A.DiveStat.Dduration(A.DiveStat.SolarEl<0)/60];
dataNight=[dataNight;B.DiveStat.PDI(B.DiveStat.SolarEl<0)/60 B.DiveStat.Dduration(B.DiveStat.SolarEl<0)/60;...
    C.DiveStat.PDI(C.DiveStat.SolarEl<0)/60 C.DiveStat.Dduration(C.DiveStat.SolarEl<0)/60;...
    D.DiveStat.PDI(D.DiveStat.SolarEl<0)/60 D.DiveStat.Dduration(D.DiveStat.SolarEl<0)/60;...
    F.DiveStat.PDI(F.DiveStat.SolarEl<0)/60 F.DiveStat.Dduration(F.DiveStat.SolarEl<0)/60;...
    G.DiveStat.PDI(G.DiveStat.SolarEl<0)/60 G.DiveStat.Dduration(G.DiveStat.SolarEl<0)/60];

[fN,xi2] = ksdensity(dataNight,xi,'Function','pdf');

%reshape output array for use in contour plot
count=1;
for i=1:size(Y2,2)
    fNb(:,i)=fN(count:count+20);
    count=count+21;
end

%plot kernel density contour lines with 6018 data points scattered on top
%(of night time dives for both). Use 9 contour lines to divide the kernel 
%density into 10 equal "slices." Use log scale on y-axis of PDI plots.

load('DayOfTripColormap.mat'); 
DiveStat=sortrows(DiveStat,8,'descend');

figure(1)
t=tiledlayout(2,2);
t.TileSpacing='compact';

nexttile
ax=gca;
hold on
[con,conax]=contour(Y,X,fd,9);
scatter(DiveStat.Dduration(DiveStat.SolarEl>=0)/60,DiveStat.Maxdepth(DiveStat.SolarEl>=0),20,(DiveStat.JulDate(DiveStat.SolarEl>=0)-min(DiveStat.JulDate)),'filled')
conax.LineWidth=2;
conax.LineColor='k';
ax.Colormap=DayOfTripColormap;
ax.YLabel.String='Dive Depth (m)';
ax.Color='w';
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=20;
ax.XLim=[0,35];
ax.YLim=[0,900];
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.Color='w';
ax.CLim=[0,100];

%Plot daytime Duration~PDI
nexttile
ax=gca;
hold on
[con,conax]=contour(Y2,X2,fDb,9);
%50% is 0.061414, 90% is 0.1228, 75% is 0.03070
scatter(DiveStat.Dduration(DiveStat.SolarEl>=0)/60,DiveStat.PDI(DiveStat.SolarEl>=0)/60,20,(DiveStat.JulDate(DiveStat.SolarEl>=0)-min(DiveStat.JulDate)),'filled')
conax.LineWidth=2;
conax.LineColor='k';
ax.Colormap=DayOfTripColormap;
c=colorbar;
c.Color='k';
c.FontSize=20;
c.Label.String='Day of Trip';
ax.Color='w';
ax.XColor='k';
ax.YColor='k';
ax.FontSize=20;
ax.XLim=[0,35];
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='PDI (min)';
ax.Title.Color='w';
ax.YLim=[0.5,400];
ax.YScale='log';
ax.CLim=[0,100];

%Plot night Depth~Duration
nexttile
ax=gca;
hold on
[con,conax]=contour(Y,X,fn,9);
scatter(DiveStat.Dduration(DiveStat.SolarEl<0)/60,DiveStat.Maxdepth(DiveStat.SolarEl<0),20,(DiveStat.JulDate(DiveStat.SolarEl<0)-min(DiveStat.JulDate)),'filled')
conax.LineWidth=2;
conax.LineColor='k';
ax.Colormap=DayOfTripColormap;
ax.YLabel.String='Dive Depth (m)';
ax.XLabel.String='Dive Duration (min)';
ax.Color='w';
ax.XColor='k';
ax.YColor='k';
ax.YDir='reverse';
ax.FontSize=20;
ax.XLim=[0,35];
ax.YLim=[0,900];
ax.LabelFontSizeMultiplier = 1.2;
ax.Title.Color='w';
ax.CLim=[0,100];

%Plot night Duration~PDI
nexttile
ax=gca;
hold on
[con,conax]=contour(Y2,X2,fNb,9);
scatter(DiveStat.Dduration(DiveStat.SolarEl<0)/60,DiveStat.PDI(DiveStat.SolarEl<0)/60,20,(DiveStat.JulDate(DiveStat.SolarEl<0)-min(DiveStat.JulDate)),'filled')
conax.LineWidth=2;
conax.LineColor='k';
ax.Colormap=DayOfTripColormap;
c1=colorbar;
c1.Color='k';
c1.FontSize=20;
c1.Label.String='Day of Trip';
ax.XLabel.String='Dive Duration (min)';
ax.Color='w';
ax.XColor='k';
ax.YColor='k';
ax.FontSize=20;
ax.XLim=[0,35];
ax.LabelFontSizeMultiplier = 1.2;
ax.YLabel.String='PDI (min)';
ax.Title.Color='w';
ax.YScale='log';
ax.YLim=[0.5,400];
ax.CLim=[0,100];


%% Load temperature data needed for Figure 4

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [26, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["CorrectedDepth", "time", "depth", "etemp", "alight"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Temp data
ts = readtable('2017001_1590113_full_iknos_raw_data.csv', opts);

ts.DateUTC = datetime(ts.time, 'ConvertFrom','datenum'); % keeps fractions of seconds

% remove fractions of seconds
[ts.Year, ts.Month, ts.Day, ts.Hour, ts.Min, ts.Sec]...
    =datevec(ts.time);
ts.Sec=round(ts.Sec);
ts.DateUTC = datetime(ts.Year, ts.Month, ts.Day, ts.Hour, ts.Min, ts.Sec);

ts = ts(:,1:6);

% check for NaN in timeseries and remove
NaNdepth = isnan(ts.CorrectedDepth);
if sum(NaNdepth)>0
    ts = ts(~NaNdepth,:);
end
clear NaNdepth

% Mat files
load('2017001_6018_RRHmod_TV3.mat');
DiveStat.DayOfTrip = ceil(DiveStat.JulDate)-floor(MetaData.DepartDate);

Track_Best.DateUTC = datetime(datestr(Track_Best.JulDate)); % fix Date format

% Crop ts to match DiveStat
DiveStat.DateUTC = datetime(DiveStat.Year, DiveStat.Month, DiveStat.Day, DiveStat.Hour,DiveStat.Min, DiveStat.Sec);
DSstart = find(ts.DateUTC>=DiveStat.DateUTC(1),1,'first');
DSend = find(ts.DateUTC>=(DiveStat.DateUTC(end)+seconds(DiveStat.Dduration(end))),1,'first');

ts=ts(DSstart:DSend,:);
clear DSstart DSend

%% Fix water temperature data
% (1) Apply correction suggested by Simmons et al. 2009 for MK9 fast-response thermistor

ts.ExtTempShift=ts.etemp;

% Subtract 0.05 degC and apply 1-s time lag (Simmons et al. 2009)
ts.ExtTempShift(2:end)=ts.etemp(1:end-1)-0.05;
ts.WaterTemp=ts.ExtTempShift;

%% Surface Temp Stats
% (2) Determine average water temperature during surface interval

% surface interval defined as ?
water_surface=cell(length(DiveStat.DiveNumber),1);
for w=1:length(DiveStat.DiveNumber)
    water_surface{w}=find(ts.DateUTC>=(DiveStat.DateUTC(w)+seconds(DiveStat.Dduration(w))) & ts.DateUTC<=(DiveStat.DateUTC(w)+seconds(DiveStat.Dduration(w))+seconds(DiveStat.PDI(w))));
end

% calculate stats of water surface layer - takes awhile to run
SurfaceTempStats = table('Size', [length(DiveStat.DiveNumber) 8], 'VariableTypes', ...
    {'double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'DiveNum','AscPts','DescPts','Avg','Std','Min','Median','Max'});

SurfaceTempStats.DiveNum = DiveStat.DiveNumber;
for i=1:length(DiveStat.DiveNumber)
    i
    wt = ts.CorrectedDepth(water_surface{i})>0;
    
    %     wt = Dive.Dry(water_surface{i})==0;
    asc_end = find(wt==0,1,'first')-1;
    desc_st = find(wt==0,1,'last')+1;
    if ~isempty(desc_st) && desc_st > length(wt)
        wt_idx = water_surface{i}([1:asc_end length(wt)]);
    else
        wt_idx = water_surface{i}([1:asc_end desc_st:length(wt)]);
    end
    if isempty(wt_idx)
        wt_idx = water_surface{i};
    else
        water_surface{i,2} = wt_idx;
        SurfaceTempStats.AscPts(i) = asc_end;
        SurfaceTempStats.DescPts(i) = length(wt)-desc_st+1;
        SurfaceTempStats.Avg(i) = mean(ts.WaterTemp(wt_idx));
        SurfaceTempStats.Std(i) = std(ts.WaterTemp(wt_idx));
        SurfaceTempStats.Min(i) = min(ts.WaterTemp(wt_idx));
        SurfaceTempStats.Median(i) = median(ts.WaterTemp(wt_idx));
        SurfaceTempStats.Max(i) = max(ts.WaterTemp(wt_idx));
        ts.WaterTemp(water_surface{i,1}) = SurfaceTempStats.Avg(i);
    end
end

%% Bottom temp stats
% surface interval defined as ?
water_bottom=cell(length(DiveStat.DiveNumber),1);
for w=1:length(DiveStat.DiveNumber)
    water_bottom{w}=find(ts.DateUTC>=(DiveStat.DateUTC(w)+seconds(DiveStat.DescTime(w)))...
        & ts.DateUTC<=(DiveStat.DateUTC(w)+seconds(DiveStat.DescTime(w))+seconds(DiveStat.Botttime(w))));
end

BottomTempStats = table('Size', [length(DiveStat.DiveNumber) 8], 'VariableTypes', ...
    {'double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'DiveNum','AscPts','DescPts','Avg','Std','Min','Median','Max'});

BottomTempStats.DiveNum = DiveStat.DiveNumber;
for i=1:length(DiveStat.DiveNumber)
    i
    wt_idx = water_bottom{i};
    water_bottom{i,2} = wt_idx;
    BottomTempStats.Avg(i) = mean(ts.WaterTemp(wt_idx));
    BottomTempStats.Std(i) = std(ts.WaterTemp(wt_idx));
    BottomTempStats.Min(i) = min(ts.WaterTemp(wt_idx));
    BottomTempStats.Median(i) = median(ts.WaterTemp(wt_idx));
    BottomTempStats.Max(i) = max(ts.WaterTemp(wt_idx));
    ts.WaterTemp(water_bottom{i,1}) = BottomTempStats.Avg(i);
end


%% Figure 4 - Fit curve to dive duration of normal animals

days=3; %number of days to cut off of end of record to eliminate "arrival" dives from curve

DurationData=[(A.DiveStat.JulDate(A.DiveStat.JulDate<(A.MetaData.ArriveDate-days))...
    -A.MetaData.DepartDate) A.DiveStat.Dduration(A.DiveStat.JulDate<(A.MetaData.ArriveDate-days))/60;...
    (B.DiveStat.JulDate(B.DiveStat.JulDate<(B.MetaData.ArriveDate-days))...
    -B.MetaData.DepartDate) B.DiveStat.Dduration(B.DiveStat.JulDate<(B.MetaData.ArriveDate-days))/60;...
    (C.DiveStat.JulDate(C.DiveStat.JulDate<(C.MetaData.ArriveDate-days))...
    -C.MetaData.DepartDate) C.DiveStat.Dduration(C.DiveStat.JulDate<(C.MetaData.ArriveDate-days))/60;...
    (D.DiveStat.JulDate(D.DiveStat.JulDate<(D.MetaData.ArriveDate-days))...
    -D.MetaData.DepartDate) D.DiveStat.Dduration(D.DiveStat.JulDate<(D.MetaData.ArriveDate-days))/60;...
    (F.DiveStat.JulDate(F.DiveStat.JulDate<(F.MetaData.ArriveDate-days))...
    -F.MetaData.DepartDate) F.DiveStat.Dduration(F.DiveStat.JulDate<(F.MetaData.ArriveDate-days))/60;...
    (G.DiveStat.JulDate(G.DiveStat.JulDate<(G.MetaData.ArriveDate-days))...
    -G.MetaData.DepartDate) G.DiveStat.Dduration(G.DiveStat.JulDate<(G.MetaData.ArriveDate-days))/60;...
    (H.DiveStat.JulDate(H.DiveStat.JulDate<(H.MetaData.ArriveDate-days))...
    -H.MetaData.DepartDate) H.DiveStat.Dduration(H.DiveStat.JulDate<(H.MetaData.ArriveDate-days))/60];

DurationData=sortrows(DurationData,1);

%fit curve to normal seal dive duration
fitobj=fit(DurationData(:,1),DurationData(:,2),'poly3');
fitci=confint(fitobj);
p22 = predint(fitobj,DurationData(:,1),0.95,'functional','off','observation','on');

%fit curve to 6018 dive duration
DiveStat.DayOfTrip = ceil(DiveStat.JulDate)-floor(MetaData.DepartDate);
DiveStat.DayOfTrip2=DiveStat.JulDate-floor(MetaData.DepartDate);
fitobj2=fit(DiveStat.DayOfTrip2,DiveStat.Dduration/60,'poly3');
p22b = predint(fitobj2,DiveStat.DayOfTrip2,0.95,'functional','off','observation','on');

%Plot PDI and Duration in two panel figure
figure(2)
t=tiledlayout(2,1);
t.TileSpacing='compact';

%Plot PDI ~ Day of Trip with temperature
nexttile
ax=gca;
s1=scatter(DiveStat.DayOfTrip2(DiveStat.SolarEl>0),DiveStat.PDI(DiveStat.SolarEl>0)/60,...
    50,SurfaceTempStats.Avg(DiveStat.SolarEl>0),'Linewidth',1.5);
hold on
s2=scatter(DiveStat.DayOfTrip2(DiveStat.SolarEl<0),DiveStat.PDI(DiveStat.SolarEl<0)/60,...
    35,SurfaceTempStats.Avg(DiveStat.SolarEl<0),'filled');
colormap(turbo)
c=colorbar;
p5=yline(2.1,'k--','LineWidth',2.5);
p6=yline(2.7,'k-','LineWidth',2.5);
ax.CLim=[9,16];
ax.XLim=[0,91];
ax.YScale='log';
ax.YLim=[1,300];
%ax.XLabel.String='Day of Trip';
ax.YLabel.String='PDI (min)';
c.Label.String = 'Surface Water Temp (\circC)';
ax.FontSize=24;
l=legend([s1 s2 p5 p6],'Day','Night','Normal','6018');
l.FontSize=16;
l.Location='Northwest';
l.NumColumns=2;

%Plot Duration ~ Day of Trip with solar elevation
nexttile
ax=gca;
scatter(DiveStat.DayOfTrip2,DiveStat.Dduration/60,20,DiveStat.SolarEl,'filled')
hold on
ax.YLabel.String='Dive Duration (min)';
ax.XLabel.String='Day of Trip';
ax.Colormap=parula;
p3=plot(DiveStat.DayOfTrip2,p22b,'k-','LineWidth',2);
p4=plot(DurationData(:,1),p22,'k--','LineWidth',2);
ax.XColor='k';
ax.YColor='k';
ax.YDir='normal';
c=colorbar;
c.Label.String='Solar Elevation';
ax.FontSize=24;
ax.CLim=[-70,80];
ax.YLim=[0,35];
ax.XLim=[0,91];
ax.LabelFontSizeMultiplier = 1.2;
