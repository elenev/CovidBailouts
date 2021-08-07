% File makes credit spread graph
% svn; 4/29/2020
% ve; 8/3/2020 (updated with new data)

close all;
clear;
clc;

forslides=false;

% FRED identifiers
ids = {'BAMLH0A0HYM2','BAMLC0A4CBBB','BAMLC0A1CAAA','DGS1','DGS5','DGS10','AAA10Y'};
varnames = {'HYspread','BBBspread','AAAspread','CMT1','CMT5','CMT10','convyield'};
labels = {'High Yield','BBB','AAA'};

% Date range
startdate = datetime(2020,1,1);
enddate = datetime();

%% Load data from FRED
fname = 'fred_spreads.mat';
reload = false;
if reload
	data = table(zeros(0,1));
	data.Properties.VariableNames{1}='matlabDate';
	c = fred();
	for i=1:length(ids)
		obj = fetch(c,ids{i},startdate,enddate);
		tmptable = array2table(obj.Data);
		tmptable.Properties.VariableNames = {'matlabDate',varnames{i}};
		data = outerjoin(data,tmptable,'Keys',{'matlabDate'},'MergeKeys',true);
	end

	% Compute datetime
	data.Date = datetime(data.matlabDate,'convertFrom','datenum');
	save(fname,'data')
else
	load(fname,'data');
end
% % Load credit spread data
% mdata      = xlsread('ICE_BofA_OAS.xlsx','FRED Graph','a14:d6177','basic'); 
% 
% Tstart = 6080;
% Tend = 6164;
% ttime=mdata(Tstart:Tend,1)+693960;
% HYspread=mdata(Tstart:Tend,2);
% BBBspread=mdata(Tstart:Tend,3);
% AAAspread=mdata(Tstart:Tend,4);

%% Plot
% For axis trimming
minmax = @(vector) [min(vector), max(vector)];

% Plot credit spreads from Jan 1 2020 until now
figure('units','normalized','outerposition',[0 0 1 1]);
if forslides
    fontfig = 16; % use 24 for slides and 14 for paper
    fileout = 'creditspreads_slides';
else
    fontfig = 11; % use 24 for slides and 14 for paper
    fileout = 'creditspreads';
end
subplot(1,3,1)
plot(data.Date,data.AAAspread,'LineWidth',1.4)
% datetick('x','mm-dd')
% xlim([2019.9,2020.4]);
ylim([0,12]);
xlim(minmax(data.Date));
ylabel('% per year')
title('AAA Spread')
set(gca,'FontSize',fontfig);
subplot(1,3,2)
plot(data.Date,data.BBBspread,'LineWidth',1.4)
% xlim([2019.9,2020.4]);
ylim([0,12]);
xlim(minmax(data.Date));
ylabel('% per year')
title('BBB Spread')
set(gca,'FontSize',fontfig);
subplot(1,3,3)
plot(data.Date,data.HYspread,'LineWidth',1.4)
% xlim([2019.9,2020.4]);
ylim([0,12]);
xlim(minmax(data.Date));
ylabel('% per year')
title('High Yield Spread')
set(gca,'FontSize',fontfig);

savefig(fileout)
printPanel(gcf,fileout,[1,3]);
print('creditspreads','-dpng')
movefile([fileout,'*'], '..\Figures');



%% Treasury and CDS graphs
spread_cols = {'Spread_1Y','Spread_5Y','Spread_10Y'};
% Load CDS sovereign spreads for the US
cdsdata = readtable('CDS_sovUS_2.xlsx');
cdsdata = cdsdata(cdsdata.InstrumentCurrency=="USD" & cdsdata.DocumentClause=="CR14",:);

% Re-scale in percent
cdsdata{:,spread_cols} = cdsdata{:,spread_cols} * 100;

data = innerjoin(data,cdsdata,'LeftKeys',{'Date'},'RightKeys',{'DataContributionDate'}, ...
	'RightVariables',{'Spread_1Y','Spread_5Y','Spread_10Y'});

% Plot Treasury yields from Jan 1 2020 until April 27 2020
figure('units','normalized','outerposition',[0 0 1 1]);
if forslides
    fileout = 'TreasCDSspreads_slides';
else
    fileout = 'TreasCDSspreads';
end
subplot(1,3,1)
plot(data.Date,data{:,{'CMT1','CMT5','CMT10'}},'LineWidth',1.4)
legend('1-yr','5-yr','10-yr')
% datetick('x','mm-dd')
% xlim([2019.9,2020.4]);
ylim([0,2]);
xlim(minmax(data.Date));
ylabel('% per year');
title('Treasury Yields');
set(gca,'FontSize',fontfig);

subplot(1,3,2)
plot(data.Date,data{:,spread_cols},'LineWidth',1.4)
legend('1-yr','5-yr','10-yr')
ylim([0,0.5]);
xlim(minmax(data.Date));
ylabel('% per year');
title('CDS Spreads');
set(gca,'FontSize',fontfig);

subplot(1,3,3)
plot(data.Date,data.convyield,'LineWidth',1.4)
ylim([0,5]);
xlim(minmax(data.Date));
ylabel('% per year');
title('Convenience Yield');
set(gca,'FontSize',fontfig);

savefig(fileout)

%% Print
printPanel(gcf,fileout,[1,3]);
movefile([fileout,'*'], '../Results');



function printPanel( f, printfile, format )
	orient(f,'landscape')
	set(f,'PaperUnits','inches');
	set(f,'PaperSize',[4*format(2) 4*format(1)]);
	set(f,'PaperPosition',[0 0 4*format(2) 4*format(1)]);
	set(f,'PaperPositionMode','manual');
	%set(gcf,'PaperSize',[10*format(1), 3*format(2)]);
	if ~isempty(printfile)
		print('-dpdf',[printfile,'.pdf']);
		print('-depsc',[printfile,'.eps']);
	end
end
