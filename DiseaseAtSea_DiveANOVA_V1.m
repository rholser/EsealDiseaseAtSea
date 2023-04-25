%Rachel Holser (rholser@ucsc.edu)
%Last Updated: 19-Sep-2022

%Compile data for ANOVAs to compare dive metrics of seal 6018 to other individuals

%Uses TV3.mat file for 6018 due to missing chunk of foieGras track 
%(affects day/night seaparation) 

%Excludes first n dives (100) from analysis due to transit off the shelf

%% Compile All Dives
clear
cd('d:\Dropbox\MATLAB\Sick Seal 6018\TV4_alpha ANNU\')
files = dir('*.mat');

pb=[];
pb_track=[];
n=100; %number of dives to remove from start of record
pb_count=0;
pbt_count=0;

for j=1:size(files,1)
%Compile tracking data, sorted into post-molt and post-breeing
    A=load(files(j).name);
    if A.MetaData.Group.Season==1
        if ~isempty(A.DiveStat)==1
            pb_count=pb_count+1;
            %DiveStats
            DayOfTrip=ceil(A.DiveStat.JulDate)-floor(A.MetaData.DepartDate);
            DayOfYear=floor(A.DiveStat.JulDate)-datenum(num2str(A.DiveStat.Year),...
                'yyyy')+1;
            temp=[DayOfTrip DayOfYear A.DiveStat.JulDate...
                A.DiveStat.Year...
                A.DiveStat.Hour...
                A.DiveStat.Maxdepth...
                A.DiveStat.Dduration...
                A.DiveStat.Botttime...
                A.DiveStat.DWigglesBott...
                A.DiveStat.SolarEl...
                A.DiveStat.Efficiency...
                A.DiveStat.DescRate...
                A.DiveStat.DescTime...
                A.DiveStat.AscRate...
                A.DiveStat.AscTime...
                A.DiveStat.TotVertDistBot...
                A.DiveStat.BottRange...
                A.DiveStat.PDI];
            temp(:,19)=A.TOPPID;
            m=size(temp,1);
            %remove first n dives
            pb=[pb;temp(n:m-n,:)];
            clear DayOfTrip DayOfYear temp TOPPID
        end
        if ~isempty(A.Track_Best)==1
            pbt_count=pbt_count+1;
            %Track
            DayOfTrip=ceil(A.Track_Best.JulDate)-floor(A.MetaData.DepartDate);
            DayOfYear=floor(A.Track_Best.JulDate)-datenum(datestr...
                (A.Track_Best.JulDate,'yyyy'),'yyyy')+1;
            temp=[DayOfTrip DayOfYear table2array(A.Track_Best(:,1:3)) table2array(A.Track_Best(:,6:8))];            
            temp(:,9)=A.TOPPID;
            pb_track=[pb_track;temp];
            clear DayOfTrip DayOfYear temp
        end
    end
    clear A
end

%Convert to table and create rounded depth and duration
pb=array2table(pb,'VariableNames',{'DayOfTrip','DayOfYear','JulDate','Year',...
    'Hour','Maxdepth','Duration','BottTime','DWigglesBott','SolarEl',...
    'Efficiency','DescRate','DescTime','AscRate','AscTime','TVD','BottRange','PDI','TOPPID'});
pb.MaxdepthR=roundn(pb.Maxdepth,1);
pb.DurationM=pb.Duration/60;
pb.BottTimeM=pb.BottTime/60;
pb.PDIM=pb.PDI/60;

%Remove data outliers (very short or extremely deep dives)
pb(pb.Duration<120,:)=[];
pb(pb.Maxdepth>1900,:)=[];
pb(pb.PDI>43200,:)=[];
pb(pb.TOPPID==2015002,:)=[];

%Convert to table and create date vector
pb_track=array2table(pb_track,'VariableNames',{'DayOfTrip','DayOfYear','JulDate','Lat',...
    'Long','Long360','TransitRate','Dist','TOPPID'});
[pb_track.Year,pb_track.Month,pb_track.Day,pb_track.Hour,pb_track.Min,...
    pb_track.Sec]=datevec(pb_track.JulDate);

%% Separate 2017 and pre-2021 and Day/Night
PB2017=pb(pb.Year==2017,:);
AllPB=pb(pb.Year<2021,:);

AllPB_Day=AllPB(AllPB.SolarEl>=0,:);
AllPB_Night=AllPB(AllPB.SolarEl<0,:);

PB2017_Day=PB2017(PB2017.SolarEl>=0,:);
PB2017_Night=PB2017(PB2017.SolarEl<0,:);

clear pm_count pb_count model model2 j m n count

%% Run ANOVA for Durations

[~,~,stats]=anova1(AllPB.DurationM,AllPB.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Day.DurationM,AllPB_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Night.DurationM,AllPB_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s


[~,~,stats]=anova1(PB2017.DurationM,PB2017.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Day.DurationM,PB2017_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Night.DurationM,PB2017_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

%% Run ANOVA for MaxDepth
[~,~,stats]=anova1(AllPB.MaxdepthR,AllPB.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Day.MaxdepthR,AllPB_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Night.MaxdepthR,AllPB_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017.MaxdepthR,PB2017.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Day.MaxdepthR,PB2017_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Night.MaxdepthR,PB2017_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

%% Run ANOVA for Bottom Time
[~,~,stats]=anova1(AllPB.BottTimeM,AllPB.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Day.BottTimeM,AllPB_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Night.BottTimeM,AllPB_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017.BottTimeM,PB2017.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Day.BottTimeM,PB2017_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Night.BottTimeM,PB2017_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

%% Run ANOVA for AscRate

[~,~,stats]=anova1(AllPB.AscRate,AllPB.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Day.AscRate,AllPB_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Night.AscRate,AllPB_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017.AscRate(PB2017.TOPPID<2017009),PB2017.TOPPID(PB2017.TOPPID<2017009));
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Day.AscRate,PB2017_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Night.AscRate,PB2017_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

%% Run ANOVA for DescRate

[~,~,stats]=anova1(AllPB.DescRate,AllPB.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Day.DescRate,AllPB_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Night.DescRate,AllPB_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017.DescRate(PB2017.TOPPID<2017009),PB2017.TOPPID(PB2017.TOPPID<2017009));
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Day.DescRate,PB2017_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Night.DescRate,PB2017_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

%% Run ANOVA for PDI

[~,~,stats]=anova1(AllPB.PDIM,AllPB.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Day.PDIM,AllPB_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(AllPB_Night.PDIM,AllPB_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017.PDIM,PB2017.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Day.PDIM,PB2017_Day.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

[~,~,stats]=anova1(PB2017_Night.PDIM,PB2017_Night.TOPPID);
multcompare(stats);
mean(stats.means)
stats.s

