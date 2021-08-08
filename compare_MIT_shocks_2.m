clear;
close all;

% Define what data to get
respath='./';
outpath='./Results/';

startecon = 'xi88midxigrid';
otherecon = 'xi88pandemic';

% Define what data to plot
from_irf =		4; % Must be a scalar
from_mit =		[1:4,6,5];
%from_mit_alt =	[6,7:9,5];
from_mit_alt = [6,7:10,5];
coloridx =      [1:4,5,6];
coloridx_alt = [5,7:10,6];

% Define labels
axlabel = {'Fin. Rec.','One-time Pandemic','New Normal'};
pollabel = {'No Covid Policy','PPP','MSLP','CCF','Combo','CBL'};
keylabel = {'donothing','ppp','mslp','ccf','combo','cbl'};

pollabel_alt = {'Combo','Combo-sized PPP','wo/ Forgiveness','wo/ Guarantees','Combo + Extra PPP','CBL'};
keylabel_alt = {'combo','pppcombosize','pppcombosizenof','pppcombosizenog','combo2xppp','cbl'};
%pollabel_alt = {'Combo','Combo + Extra PPP','Combo + Transfers'};
%keylabel_alt = {'combo','combo2xppp','combotrasnfers'};

% Plot the first shock (the one from the IRFs)
leave_off_noshock = true;

N_shocks = 2; %startecon, otherecon -- hard coded
N_policies = length(from_mit);
N_policies_alt = length(from_mit_alt);

% Define file names
start_resfile = ['res_',startecon];
resfile = ['res_',otherecon];

outfile=['GRbar_',resfile];

% Load files
irfs = load([respath, 'GR_',    start_resfile, '.mat']);
mit1 = load([respath, 'MIT_',   start_resfile, '.mat']);
mit2 = load([respath, 'MITPT_', resfile,       '.mat']);

mit1.simseries_mean = mit1.simseries_mean(1:11);

N_irf = numel( irfs.simseries_mean );
N_mit = numel( mit1.simseries_mean );

plot_shocks = [from_irf, ...
			   N_irf + from_mit, ...
			   N_irf + N_mit + from_mit ];
		   
plot_shocks_alt = [from_irf, ...
			   N_irf + from_mit_alt, ...
			   N_irf + N_mit + from_mit_alt ];

% Extract and arrange imported series
simseries_mean_start = irfs.simseries_mean;
simseries_mean = [simseries_mean_start; mit1.simseries_mean; mit2.simseries_mean];
qC = [nan(size(simseries_mean_start)); mit1.qC; mit2.qC];
indexmap = irfs.indexmap;


% struct as input for plotting functions
simdata.series=simseries_mean;
simdata.names=indexmap;
simdata.N_shocks=N_shocks;
simdata.N_policies=length(pollabel);


% Define transforms of series
Y0 = simseries_mean{from_irf+1}(1,indexmap.get('Y'));

      levels.fcn  = @(xvec,basevec) 100 * xvec;
    levelsY0.fcn  = @(xvec,basevec) 100 * xvec / Y0;
  arithdelta.fcn  = @(xvec,basevec) 100 * (xvec-xvec(1,:) );
arithdeltaY0.fcn  = @(xvec,basevec) 100 * (xvec-xvec(1,:))/ Y0;
   geomdelta.fcn  = @(xvec,basevec) 100 * (xvec./xvec(1,:) - 1);
 
      levels.lab  = '';
    levelsY0.lab  = 'Pct of t=0 GDP';
  arithdelta.lab  = 'Pp Change from t=0';
arithdeltaY0.lab  = 'Change from t=0 As Pct of t=0 GDP';
   geomdelta.lab  = 'Pct Change from t=0';

cumarithdelta.fcn   = @(xvec,basevec) 100 * cumsum( xvec - basevec, 'reverse' );
cumgeomdelta.fcn   = @(xvec,basevec) 100 * cumsum( xvec - basevec, 'reverse' ) ./ xvec(1,:);
 
cumarithdelta.lab   = 'Cum Pp Change from t=0';
 cumgeomdelta.lab   = 'Cum Pct Change from t=0';


%% Compute welfare
baseline_idx = [(N_irf+1)*ones(1,N_mit), (N_irf+N_mit+1)*ones(1,N_mit)];

qCB = arrayfun(@(i)qC{i}.qCB, baseline_idx);
qCS = arrayfun(@(i)qC{i}.qCS, baseline_idx);

VB_policy = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VB')), N_irf + (1:2*N_mit) );
VS_policy = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VS')), N_irf + (1:2*N_mit) );

VB_base = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VB')), baseline_idx);
VS_base = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VS')), baseline_idx);

cvwelfare_all = (VB_policy ./ VB_base - 1) .* qCB + ...
			 (VS_policy ./ VS_base - 1) .* qCS;

cvwelfare = cvwelfare_all(plot_shocks(2:end) - N_irf);

%% Compute welfare (alternative policies relative to combo)
basealt_idx = [(N_irf+from_mit_alt(1))*ones(1,N_mit), (N_irf+N_mit+from_mit_alt(1))*ones(1,N_mit)];
qCB_alt = arrayfun(@(i)qC{i}.qCB, basealt_idx);
qCS_alt = arrayfun(@(i)qC{i}.qCS, basealt_idx);
VB_basealt = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VB')), basealt_idx);
VS_basealt = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VS')), basealt_idx);

cvwelfare_allalt = (VB_policy ./ VB_basealt - 1) .* qCB_alt + ...
			 (VS_policy ./ VS_basealt - 1) .* qCS_alt;

cvwelfare_alt = cvwelfare_allalt(plot_shocks_alt(2:end) - N_irf);

%% Define plots for paper
%varsel{1} = {'Y','C','X','Gmlbt'};
%varsel{1} = {'Y','C','X','tax_byY'};
varsel{1} = {'Y','C','X'};
%varsel{1} = {};
varsel{2} = {'C','X','DWL_byY'};
varsel{3} = {'bailout','programspend','tax','rD','Gdef_prim','Gmlbt'}; %
varsel{4} = {'Drate','Lrate','fracS','rB','Lspr','LsprT'};
varsel{5} = {'brupt','Ibkasset','WI_byY'};
varsel{6} = {'cB','cS','C','VB','VS','cvwelfare'};
varsel{7} = {'fiscspend','Gdbtsvc','Gdef_prim','Gmlbt'}; %

%titles{1} = {'GDP','Consumption','Investment','Govt Debt'};
titles{1} = {'GDP','Consumption','Investment'};
titles{2} = {'Consumption','Investment','DWL / GDP'};
titles{3} = {'Bailouts','Cost of Program','Tax Revenue', ... %
				'Safe Rate','Primary Deficit','Govt Debt'};
titles{4} = {'Default Rate','Loss Rate','Frac Savers','Loan Rate','Loan Spread','Dur-Adj Loan Spread'};
titles{5} = {'Intermediary Failures','Intermediary Assets','Intermediary Net Worth / GDP'};
titles{6} = {'Consumption, B','Consumption, S','Aggr. Consumption', ...
				'Value Fcn, B', 'Value Fcn, S', 'CEV Welfare Rel. to No Covid Policy'};
titles{7} = {'Fiscal Spend','Debt Service','Primary Deficit','Mkt Val Debt'};
 
labelYaxis = {false, true, true, false, true, true, true};

transforms{1} = { cumgeomdelta, cumgeomdelta, cumgeomdelta, [] };
transforms{2} = { geomdelta, geomdelta, arithdelta };
transforms{3} = {  arithdeltaY0, levelsY0, arithdeltaY0, ...
					levels, levelsY0, arithdeltaY0 };
transforms{4} = { levels, levels, levels, levels, levels, levels };
transforms{5} = { levels, geomdelta, arithdelta };
transforms{6} = { geomdelta, geomdelta, geomdelta, geomdelta, geomdelta, levelsY0 };
transforms{7} = { arithdeltaY0, arithdeltaY0, levelsY0, arithdeltaY0 };

%irftransforms{1} = {geomdelta, geomdelta, geomdelta, arithdeltaY0};

experiments = {'newnormal','newnormal','newnormal','newnormal','newnormal','newnormal','newnormal'};
includePandemicTransition = {false, true, true, true, true, true, true};


%% Plot MIT bars
close all;

plotspec.varsel=varsel;
plotspec.titles=titles;
plotspec.transforms=transforms;
plotspec.axlabel=axlabel;
plotspec.pollabel=pollabel;
plotspec.keylabel=keylabel;
plotspec.labelYaxis=labelYaxis;
plotspec.experiments=experiments;
plotspec.includePandemicTransition=includePandemicTransition;
plotspec.fsize=12;
plotspec.coloridx=coloridx;

[~,cellDict_newnormal]=plotMITShocks(simdata,cvwelfare,plot_shocks,plotspec,[outpath,outfile,'_MITbar'],{startecon,otherecon});

plotspec.experiments = {'onetime','onetime','onetime','onetime','onetime','onetime','onetime'};
[~,cellDict_onetime]=plotMITShocks(simdata,cvwelfare,plot_shocks,plotspec,[outpath,outfile,'_MITbarOnetime'],{startecon,otherecon});

%% Plot MIT bars with alternative programs
simdata_alt = simdata;
simdata_alt.N_policies = N_policies_alt;
plotspec_alt.varsel= { {'Drate','X','brupt','programspend','Gmlbt','cB','cS','cvwelfare' } };
plotspec_alt.titles={ {'Default Rate','Investment','Intermediary Failures','Cost of Corp Program', ...
					'Govt Debt','Consumption B','Consumption S','CEV Welfare Rel. to Combo' } };
plotspec_alt.transforms={ { levels, geomdelta, levels, arithdeltaY0, arithdeltaY0, geomdelta, geomdelta, levelsY0} };
plotspec_alt.axlabel=axlabel;
plotspec_alt.pollabel=pollabel_alt;
plotspec_alt.keylabel=keylabel_alt;
plotspec_alt.labelYaxis={ true };
plotspec_alt.experiments={ 'newnormal' };
plotspec_alt.fsize=12;
plotspec_alt.coloridx=coloridx_alt;

[~,cellDictAlt]=plotMITShocks(simdata_alt,cvwelfare_alt,plot_shocks_alt,plotspec_alt,[outpath,outfile,'_MITbarAlt']);


%% Store MIT results in simres
cellDict = [cellDict_newnormal; cellDict_onetime, cellDictAlt];
load([outpath,'/simres.mat'],'mapDict');

% merge and save
newKeys = cellDict(:,1);
if any( isKey( mapDict, newKeys ) )
	warning('simres Key names already exist. Overwriting.');
end
addlmapDict = containers.Map( newKeys,cellDict(:,2) );
mapDict = [mapDict; addlmapDict];
save([outpath,'/simres.mat'],'mapDict');

%% Plot IRFs
% Only plot panel 1
irfvarsel = {'Y','C','X','Gmlbt'};
irftitles = {'GDP','Consumption','Investment','Gov. Debt'};
irftransforms = {geomdelta, geomdelta, geomdelta, geomdelta};
% Only plot some policies
irfplotidx = [1,N_policies-1,N_policies];
irfplotidx_list = {irfplotidx,N_policies+irfplotidx};

plotspecIRF.varsel=irfvarsel;
plotspecIRF.titles=irftitles;
plotspecIRF.transforms=irftransforms;
plotspecIRF.pollabel=pollabel(irfplotidx);
plotspecIRF.irfplotidx=irfplotidx_list{1};
plotspecIRF.coloridx=coloridx(irfplotidx);
plotspecIRF.fsize=10;

plotMITIRF(simdata,plot_shocks,plotspecIRF,[outpath,outfile,'_MITIRF1']);

plotspecIRF.irfplotidx=irfplotidx_list{2};

plotMITIRF(simdata,plot_shocks,plotspecIRF,[outpath,outfile,'_MITIRF2']);

%% Plot MIT bars for slides

varsel=[];
varsel{1} = {'Y','C','X'};
%varsel{1} = {};
varsel{2} = {'C','X','DWL_byY','Gmlbt'};
varsel{3} = {'programspend_byY','bailout_byY','rD','Gdef_prim_byY'}; %
varsel{4} = {'Drate','Lrate','LsprT'};
varsel{5} = {'brupt','Ibkasset','WI_byY'};
varsel{6} = {'cB','cS','cvwelfare'};

titles=[];
titles{1} = {'GDP','Consumption','Investment'};
%titles{1} = {'GDP','Consumption','Investment','Tax Revenue / GDP'};
titles{2} = {'Consumption','Investment','DWL / GDP','Govt Debt'};
titles{3} = {'Cost of Program / GDP','Bailouts / GDP','Safe Rate','Primary Deficit/GDP'};
titles{4} = {'Default Rate','Loss Rate','Loan Spread'};
titles{5} = {'Intermediary Failures','Intermediary Assets','Intermediary Net Worth / GDP'};
titles{6} = {'Consumption, B','Consumption, S','CEV Welfare Rel. to Do Nothing'};

           
transforms=[];
transforms{1} = { cumgeomdelta, cumgeomdelta, cumgeomdelta, [] };
transforms{2} = { geomdelta, geomdelta, arithdelta, geomdelta };
transforms{3} = { levels, levels, levels, levels};
transforms{4} = { levels, levels, levels };
transforms{5} = { levels, geomdelta, arithdelta };
transforms{6} = { geomdelta, geomdelta, levelsY0 };
       
labelYaxis = {false, true, false, false, true, true};

close all;
plotspec.varsel=varsel;
plotspec.titles=titles;
plotspec.experiments = {'newnormal', 'newnormal', 'newnormal', 'newnormal', 'newnormal', 'newnormal'};
plotspec.transforms=transforms;
plotspec.axlabel=axlabel;
plotspec.pollabel=pollabel;
plotspec.labelYaxis=labelYaxis;
plotspec.fsize=14;
plotspec.coloridx=coloridx;

plotMITShocks(simdata,cvwelfare,plot_shocks,plotspec,[outpath,outfile,'_MITbarslides']);

% Do again, but for the "Onetime"
plotspec.experiments = {'onetime', 'onetime', 'onetime', 'onetime', 'onetime', 'onetime'};
plotMITShocks(simdata,cvwelfare,plot_shocks,plotspec,[outpath,outfile,'_MITbarslidesOnetime']);

%% Plot MIT bars with alternative programs for slides
varsel=[]; transforms=[]; titles=[];

varsel{1} = {'Drate','brupt','programspend','Gmlbt'};
varsel{2} = {'cB','cS','X','cvwelfare' };

transforms{1} = { levels, levels, arithdeltaY0, arithdeltaY0 };
transforms{2} = { geomdelta, geomdelta, geomdelta, levelsY0  };

titles{1} = {'Default Rate','Intermediary Failures','Cost of Corp Program', 'Govt Debt'} ;
titles{2} = {'Consumption B','Consumption S','Investment','CEV Welfare Rel. to Combo' };

plotspec_alt.varsel= varsel;
plotspec_alt.titles=titles;
plotspec_alt.transforms=transforms;
plotspec_alt.axlabel=axlabel;
plotspec_alt.pollabel=pollabel_alt;
plotspec_alt.labelYaxis={ true, true };
plotspec_alt.experiments={ 'newnormal', 'newnormal' };
plotspec_alt.fsize=14;
plotspec_alt.coloridx=coloridx_alt;

plotMITShocks(simdata_alt,cvwelfare_alt,plot_shocks_alt,plotspec_alt,[outpath,outfile,'_MITbarAlt']);



%% Plot IRFs for slides

irfvarsel = {'Y','C','X','Gmlbt','DWL_byY','Gdef_prim_byY'};
irftitles = {'GDP','Consumption','Investment','Gov. Debt','DWL / GDP','Primary Deficit/GDP'};
irftransforms = {geomdelta, geomdelta, geomdelta, geomdelta, arithdelta, arithdelta};

irfplotidx = [1,N_policies-1,N_policies];
irfplotidx_list = {N_policies+irfplotidx(1),N_policies+irfplotidx(1:2),N_policies+irfplotidx};

plotspecIRF.varsel=irfvarsel;
plotspecIRF.titles=irftitles;
plotspecIRF.transforms=irftransforms;
plotspecIRF.fsize=14;

for n=1:length(irfplotidx_list)
    plotspecIRF.irfplotidx=irfplotidx_list{n};
	plotspecIRF.pollabel=pollabel(irfplotidx(1:n));
    plotspecIRF.coloridx=coloridx(irfplotidx(1:n));
    plotMITIRF(simdata,plot_shocks,plotspecIRF,[outpath,outfile,'_MITIRFslide',num2str(n)]);
end

irfvarsel = {'cB','cS','fracS','WSm'};
irftitles = {'Cons B','Cons S','Saver corp. share (EOP)','Saver wealth (BOP)'};
irftransforms = {geomdelta, geomdelta, arithdelta, geomdelta};

plotspecIRF.varsel=irfvarsel;
plotspecIRF.titles=irftitles;
plotspecIRF.transforms=irftransforms;
plotspecIRF.irfplotidx=irfplotidx_list{3};
plotspecIRF.coloridx=coloridx(irfplotidx(1:n));
plotMITIRF(simdata,plot_shocks,plotspecIRF,[outpath,outfile,'_MITIRFslidecons']);

close all;
