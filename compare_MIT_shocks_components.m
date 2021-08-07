clear;
close all;

% Define what data to get
respath='./';
outpath='./Results/';

startecon = 'xi88midxigrid';
otherecon = 'xi88pandemic';

% Define what data to plot
from_irf =		1; % Must be a scalar
from_mit =		[12:13,1];
from_mitpt =	1;
coloridx =      [1:4,5,6];
coloridx_alt = [5,7:10,6];

% Define labels
pollabel = {'Markov','+ Idiosync. Vol','+ Idiosync. Mean','+ New Normal'};
keylabel = {'markov','vol','mean','newnormal'};

% Define file names
start_resfile = ['res_',startecon];
resfile = ['res_',otherecon];

outfile=['GRbar_',resfile];

% Load files
irfs  = load([respath, 'GR_',    start_resfile, '.mat']);
mit   = load([respath, 'MIT_',   start_resfile, '.mat']);
mitpt = load([respath, 'MITPT_', resfile,       '.mat']);

N_irf = numel( irfs.simseries_mean );
N_mit = numel( mit.simseries_mean );

% Extract and arrange imported series
simseries_mean_start = irfs.simseries_mean;
simseries_mean = [simseries_mean_start(from_irf); mit.simseries_mean(from_mit); mitpt.simseries_mean(from_mitpt)];
plot_shocks = 1:length(simseries_mean);

qC = [nan(size(simseries_mean_start)); mit.qC; mitpt.qC];
indexmap = irfs.indexmap;

N_policies = length(plot_shocks)-1;

% struct as input for plotting functions
simdata.series=simseries_mean;
simdata.names=indexmap;
simdata.N_shocks=1;
simdata.N_policies=N_policies;


% Define transforms of series
Y0 = simseries_mean{from_irf}(1,indexmap.get('Y'));

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
baseline_idx = ones(1,N_policies);

%qCB = arrayfun(@(i)qC{i}.qCB, baseline_idx);
%qCS = arrayfun(@(i)qC{i}.qCS, baseline_idx);

VB_policy = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VB')), plot_shocks(2:end) );
VS_policy = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VS')), plot_shocks(2:end) );

VB_base = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VB')), baseline_idx);
VS_base = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('VS')), baseline_idx);

qCB = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('qCB')), baseline_idx);
qCS = arrayfun(@(i)simseries_mean{i}(2,indexmap.get('qCS')), baseline_idx);

cvwelfare = (VB_policy ./ VB_base - 1) .* qCB + ...
			 (VS_policy ./ VS_base - 1) .* qCS;

%% Define plots for paper
%varsel{1} = {'Y','C','X','Gmlbt'};
%varsel{1} = {'Y','C','X','tax_byY'};
varsel{1} = {'Y','C','X'};
%varsel{1} = {};
varsel{2} = {'Y','C','X','Drate','WI_byY','Gmlbt'};

%varsel{7} = {'tax_byY','Gdbtsvc','Gdef_prim','Gmlbt'}; %

%titles{1} = {'GDP','Consumption','Investment','Govt Debt'};
titles{1} = {'GDP','Consumption','Investment'};
titles{2} = {'GDP','Consumption','Investment','Default Rate','Intermediary Net Worth / GDP','Govt Debt'};


 
labelYaxis = {true, true, true, false, true, true};

transforms{1} = { cumgeomdelta, cumgeomdelta, cumgeomdelta, [] };
transforms{2} = { geomdelta, geomdelta, geomdelta, levels, arithdelta, arithdeltaY0};


%irftransforms{1} = {geomdelta, geomdelta, geomdelta, arithdeltaY0};

experiments = {'onetime','onetime'};
includePandemicTransition = {false, false};


%% Plot MIT bars
close all;

plotspec.varsel=varsel(2:end);
plotspec.titles=titles(2:end);
plotspec.transforms=transforms(2:end);
plotspec.axlabel={'None','Shock Components'};
plotspec.pollabel=pollabel;
plotspec.keylabel=keylabel;
plotspec.labelYaxis=labelYaxis;
plotspec.experiments=experiments;
plotspec.includePandemicTransition=includePandemicTransition;
plotspec.fsize=12;
plotspec.coloridx=coloridx;

plotMITShocks(simdata,cvwelfare,plot_shocks,plotspec,[outpath,outfile,'_MITbarComponents'],{startecon,otherecon});