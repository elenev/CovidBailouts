if usejava('desktop')
   clear; 
else
    ver;
end
close all;

%% Define exercise

% file with model
respath='./';
outpath='./Results/';

%resfile='res_20200403_nuS02xi88_s130';
%start_resfile = resfile;

%resfile = 'res_20210210_nuS02xi88pandemic_s130';
%start_resfile = 'res_20200403_nuS02xi88_s130';

%start_resfile = 'res_20210210_xi88midxigrid_s130';
start_resfile = 'res_20210210_xi88safe_s130';

%resfile = 'res_20210210_xi88pandemic_s130';
resfile = start_resfile;

% number of periods
N_runs=10000;
NT_sim=25;

full_paths = true;

% output table file
if full_paths
	prefix = 'MIT';
else
	prefix = 'MITshort';
end

if ~strcmp(resfile,start_resfile)
	prefix = [prefix,'PT'];
end
outfile=[prefix,'_',resfile];

import_initial_guess_from_PT = false;

%% Load files
varlist={'simseries','statevec','indexmap','varnames','params'};

% Load to res file
load([respath,resfile,'.mat']);

% Load from sim_res file
load(['sim_',start_resfile],varlist{:});
varnames_store = varnames;

% Load from res file
if ~strcmp(resfile,start_resfile)
	tmp = load(start_resfile,'mobj');
	startmobj=tmp.mobj;
	
	% Temporarily load sim_res of new 
	tmp = load(['sim_',resfile],varlist{:});
	to_simseries = tmp.simseries;
	to_statevec = tmp.statevec;
else
	startmobj=mobj;
	to_simseries = simseries;
	to_statevec = statevec;
end

if import_initial_guess_from_PT
	PT = load(['PT_',resfile],'simseries_mean');
end

%% Set starting point
gv = @(x)indexmap.get(x);
start_ini=5;
% orig_statevec=statevec(2:end);

% N_vars=length(startvals);
% simulate for one period
NT_sim_ini=1;
%warning('NT_sim_ini = 0');
Nvarsim=startmobj.NSTEX + startmobj.NSTEN + startmobj.NSOL + startmobj.NV + ...
			startmobj.NADD + startmobj.NCOND + startmobj.NZNS;
		
% Find states such that statevec = start_ini
% Note, statevec(2) corresponds to simseries(1)
idx_ini = find(statevec(2:end)==start_ini);

% Define initial values
startvals=mean(simseries(idx_ini,:));

% Simulate forward to get actual, not averaged, values for variables other
% than state variables
if NT_sim_ini>0
	firstrow=[start_ini, startvals(1:Nvarsim)];
%	startpt_vec=startvals([1,1+[gv('KB'),gv('LB'),gv('WI'),gv('BG')]]);
	startpt_vec=[start_ini,startvals([gv('KB'),gv('LB'),gv('WI'),gv('BG')])];
	transprob=cumsum(startmobj.Exogenv.mtrans(start_ini,:));
	shock_prob=transprob(start_ini);
	if start_ini>1
		shock_prob_minus=transprob(start_ini-1);
	else
		shock_prob_minus=0;
	end
	rvar_next=(shock_prob+shock_prob_minus)/2;
	[simseries_ini,varnames_ini,~,~,~]=startmobj.simulate(NT_sim_ini,0,startpt_vec,0,rvar_next*ones(NT_sim_ini,1));

	[simseries_ini, varnames_ini] = startmobj.computeSimulationMoments(simseries_ini,varnames_ini,[],firstrow);
	startvals=simseries_ini(end,:);
end
N_vars=length(startvals);

% initial point in terms of quantities
% Because we don't have start-of-period AI or bI, use the previous period's
% end-of-period variables net of growth rate
%enstates = simseries(:,[gv('KB'),gv('LB'),gv('WI'),gv('BG')]);
quantenstates = startvals([gv('KB_g'), gv('AB_g'), gv('AI_g'), gv('bI'), gv('BG_g')]) / startmobj.Params.mu_G;
quantpoint = [start_ini, quantenstates];

% Prep initial guess for MIT state
guess_statevec=to_statevec(2:end);
guess_enstates = to_simseries(:, [ gv('KB'), gv('LB'), gv('WI'), gv('BG') ] );

% Compute average wagebill
wagebill=simseries(:,indexmap.get('wbill'));

%% Define shocks
%Shock to bank wealth
%shockcell{1}.exst = 4;
%shockcell{1}.params = {'WIDWL',-0.8381 * startvals(:,gv('WI'))};

% Define pure shocks
om = 0.25;
mu_om = 1-0.036; %0.9550;
%rho = 0.02;
if params.rho==0
	rho = 0.04;
else
	rho = params.rho;
end
sigma_eps=.05;
lscale = 0; %0.05;


if strcmp(resfile,start_resfile)
	% Use this shock when res_ file does NOT have a pandemic state
	pureshockcell{1}.exst = 6;
	pureshockcell{1}.params = {'rho',rho;'sigma_eps',sigma_eps;};
	pureshockcell{1}.pts = {'Om',om^2;'mu_om',mu_om,;'Lscale',mobj.Params.Lscale(1)*(1-lscale);};

	pureshockcell{2}.exst = 6;
	pureshockcell{2}.params = {'rho',rho;'sigma_eps',sigma_eps;};
	pureshockcell{2}.pts  ={};
	
	pureshockcell{3}.exst = 6;
	pureshockcell{3}.params = {'rho',rho;'sigma_eps',sigma_eps;};
	pureshockcell{3}.pts  ={'Om',om^2};
	
	pureshockcell{4}.exst = 6;
	pureshockcell{4}.params = {'rho',rho;'sigma_eps',sigma_eps;};
	pureshockcell{4}.pts  ={'Om',om^2;'mu_om',mu_om};
	
	pureshockcell{5}.exst = 6;
	pureshockcell{5}.params = {'rho',rho;'sigma_eps',sigma_eps;};
	pureshockcell{5}.pts  ={'Om',om^2;'mu_om',mu_om;'Lscale',mobj.Params.Lscale(1)*(1-lscale);};

else
	% Use this shock when res_ file does have a pandemic state
	pureshockcell{1}.exst = 11;
	pureshockcell{1}.params = {'rho',rho;'sigma_eps',sigma_eps;}; %
	pureshockcell{1}.pts = {};
end


% pureshockcell{2}.exst = 11;
% pureshockcell{2}.params = {'rho',rho};
% pureshockcell{2}.pts = {};

% 
% pureshockcell{2}.exst = exst;
% pureshockcell{2}.params = {'rho',rho};
% pureshockcell{2}.pts = {'Om',om^2;'mu_om',mu_om};

% pureshockcell{3}.exst = exst;
% pureshockcell{3}.params = {'rho',rho};
% pureshockcell{3}.pts = {'Om',om^2;'mu_om',mu_om; ...
% 	'Lscale',mobj.Params.Lscale(1)*(1-lscale);};

% Define policies
% ombarT=0.93;
% disp(['Targeted firms: ', num2str(1-ProductionModel.Mpayoff_gamma(ombarT,mu_om,om^2))]);

ppp_size = 0.031;
mslp_size = 0.028;

ccf_size_gdp = 0.0598; % qB*AG_g in the combo economy

ppp_guarantee = 1;
mslp_guarantee = 0.95;

ppp_rate = 0.01;
mslp_rate = 0.03;

ppp_forgive = -1;
mslp_forgive = 0;

wgtavg = @(ppp,mslp) ( ppp_size * ppp + mslp_size * mslp ) / ( ppp_size + mslp_size);
wgtavgN = @(ppp,mslp,n) ( n*ppp_size * ppp + mslp_size * mslp ) / ( n*ppp_size + mslp_size);

policies{1}.params = {}; % No policy
policies{2}.params = {'ABridge',ppp_size;'rBridge',ppp_rate;'fBridge',ppp_forgive; ...
					  'Ibr',ppp_guarantee;'ombar',Inf}; % Untargeted PPP
policies{3}.params = {'ABridge',mslp_size;'rBridge',mslp_rate;'fBridge',mslp_forgive; ...
					  'Ibr',mslp_guarantee;'ombar',Inf}; % Untargeted MSLP
%policies{4}.params = {'ABridge',.06;'rBridge',.02;'fBridge',-.5;'Ibr',0.975;'ombar',Inf}; % Untargeted both
%policies{5}.params = {'ABridge',.03;'rBridge',.01;'fBridge',-1;'Ibr',1;'ombar',.93}; % Targeted PPP
%policies{6}.params = {'ABridge',.03;'rBridge',.03;'fBridge',0;'Ibr',0.95;'ombar',.93}; % Targeted MSLP
%policies{7}.params = {'ABridge',.06;'rBridge',.02;'fBridge',-.5;'Ibr',0.975;'ombar',.93}; % Targeted both
policies{4}.params = {'fracG',0.089}; % Bond buying
policies{5}.params = {'bridge',.075;'rBridge',wgtavg(ppp_rate,mslp_rate); ...
				      'fBridge',wgtavg(ppp_forgive,mslp_forgive); ...
					  'Ibr',wgtavg(ppp_guarantee,mslp_guarantee);'ombar',0.67}; % Conditional both
policies{6}.params = {'ABridge',ppp_size+mslp_size; ...
					  'rBridge',wgtavg(ppp_rate,mslp_rate); ...
					  'fBridge',wgtavg(ppp_forgive,mslp_forgive); ...
					  'Ibr',wgtavg(ppp_guarantee,mslp_guarantee); ...
					  'ombar',Inf;'fracG',0.089}; % Combo: Untargeted + bond buying

policies{7}.params = {'ABridge',ppp_size+mslp_size+ccf_size_gdp; ...
					  'rBridge',ppp_rate; ...
					  'fBridge',ppp_forgive; ...
					  'Ibr',ppp_guarantee; ...
					  'ombar',Inf}; % Combo-sized PPP
policies{8}.params = {'ABridge',ppp_size+mslp_size+ccf_size_gdp; ...
					  'rBridge',ppp_rate; ...
					  'fBridge',0; ...
					  'Ibr',ppp_guarantee; ...
					  'ombar',Inf}; % Combo-sized PPP without forgiveness			  
policies{9}.params = {'ABridge',ppp_size+mslp_size+ccf_size_gdp; ...
					  'rBridge',ppp_rate; ...
					  'fBridge',ppp_forgive; ...
					  'Ibr',0; ...
					  'ombar',Inf}; % Combo-sized PPP without guarantees
				  
policies{10}.params = {'ABridge',2*ppp_size+mslp_size; ...
					  'rBridge',wgtavgN(ppp_rate,mslp_rate,2); ...
					  'fBridge',wgtavgN(ppp_forgive,mslp_forgive,2); ...
					  'Ibr',wgtavgN(ppp_guarantee,mslp_guarantee,2); ...
					  'ombar',Inf;'fracG',0.089}; % Combo w 2x PPP: Untargeted + bond buying		 
policies{11}.params = {'ABridge',ppp_size+mslp_size; ...
					  'rBridge',wgtavg(ppp_rate,mslp_rate); ...
					  'fBridge',wgtavg(ppp_forgive,mslp_forgive); ...
					  'Ibr',wgtavg(ppp_guarantee,mslp_guarantee); ...
					  'ombar',Inf;'fracG',0.089; ...
					  'T',mobj.Params.T+0.0276}; % Combo: Untargeted + bond buying				  
				  
% Combine
N_s = numel(pureshockcell);
N_p = numel(policies);
N_shock = N_p + N_s - 1;
shockcell = cell( N_shock, 1 );
for ss = 1:N_s
	if ss==1
		for pp = 1:N_p
			idx = N_p*(ss-1) + pp;
			shockcell{ idx }.exst = pureshockcell{ss}.exst;
			shockcell{ idx }.params = [ pureshockcell{ss}.params; ...
										policies{pp}.params ];
			shockcell{ idx }.pts = pureshockcell{ss}.pts;
		end
	else
		idx = N_p + ss - 1;
		shockcell{ idx }.exst = pureshockcell{ss}.exst;
		shockcell{ idx }.params = pureshockcell{ss}.params;
		shockcell{ idx }.pts = pureshockcell{ss}.pts;
	end
end
				   
%% Initialize variables to stores impulse responses
simseries_median = cell(N_shock,1);
simseries_mean = cell(N_shock,1);
simseries_std = cell(N_shock,1);
qC = cell(N_shock,1);

%% Create parallel pool
if usejava('desktop') && full_paths
    if ~exist('no_par_processes','var')
        no_par_processes=min( 16, feature('numcores') );
	end
	cp=gcp('nocreate');
    if ~isempty(cp) 
        if  cp.NumWorkers~=no_par_processes
            delete(cp);
            if no_par_processes>0
                parpool(no_par_processes);
            end
        end
    else
        if no_par_processes>0
            parpool(no_par_processes);
        end
    end
elseif full_paths
    scheduler = parcluster('local');
    parpool(scheduler, str2double(getenv('NTHREADS')));
end

N_runs = N_runs * full_paths + 1 * (~full_paths);

%% Launch
for s=1:N_shock
	disp(['Shock ',num2str(s),' of ',num2str(N_shock)]);
    tens_simseries = zeros(NT_sim+1,N_vars,N_runs);

    % compute entry of random number matrix that sets first state
    % deterministically to start_shock
	start_shock = shockcell{s}.exst;

    % Create shock matrix
    rng(1);
    %shmatfull = rand(NT_sim*N_runs,1);
    shmatfull = lhsdesign(N_runs,NT_sim);

	exnpt=mobj.Exogenv.exnpt;
    SDFmat=zeros(NT_sim,exnpt);
	
	% Compute first period economy
%	idx_switch = 1 + find( orig_statevec(2:end)==start_shock & orig_statevec(1:end-1)==start_ini);
	idx_switch = 1 + find( guess_statevec(2:end)==start_shock & guess_statevec(1:end-1)==start_ini);
	if import_initial_guess_from_PT
		guessenstates = PT.simseries_mean{1}(2,[gv('KB'),gv('LB'),gv('WI'),gv('BG')]);
	else
		guessenstates = mean(guess_enstates(idx_switch,:)); 
	end
	
	[transmat,simvec,mitparams,qC{s}] = computeFirstPeriod(mobj,quantpoint,guessenstates,shockcell{s},0);
	
	paramNames = fieldnames(mitparams);
	periodsAfterShock = full_paths * (NT_sim-1);
	for p = paramNames'
		p = p{:};
		origval = mobj.Params.(p);
		mitval = mitparams.(p);
		if ~isa( origval, 'function_handle')
			origval = origval(:)';
			mitval = mitval(:)';
			mitparams.(p) = [origval; mitval; repmat( origval, periodsAfterShock, 1 )];
		end
	end
	
	% Select Markov transition probabilities
	transprob=cumsum(mobj.Exogenv.mtrans(start_shock,:));
	
	if full_paths
		fprintf([repmat('.',1,100) '\n\n']);
		parfor n=1:N_runs       
		%for n=1:N_runs
			%--------------------------------------------------------------------------
			% start simulation
			%--------------------------------------------------------------------------
			%fprintf('Run %d - Start \n',n);
			% simulate
			shmat = shmatfull(n,:)';


			exnext=find(transprob-shmat(1)>0,1,'first');
			startpt_vec = transmat(exnext,:);

			[simseries,varnames,~,~,~]=mobj.simulate(NT_sim-1,0,[exnext,startpt_vec],0,shmat(2:end));
			simseries = [simvec;simseries];

			simseries_orig=simseries;
			varnames_orig=varnames;
			statevec = simseries(:,1);
			%fprintf('Run %d - After simulation \n',n);
			startvals_orig = [start_ini,startvals(1:size(simseries,2)-1)];

			[simseries, varnames] = mobj.computeSimulationMoments(simseries,varnames,[],startvals_orig, mitparams);

			nvars = length(varnames);
			%fprintf('Run %d - After computation \n',n);
	%         disp(size(startvals))
	%         disp(size(simseries))
			tens_simseries(:,:,n) = [startvals; simseries];
			if mod(n,N_runs/100)==0
				%disp([num2str(n),'/',num2str(N_runs),': ',num2str(round(1000*n/N_runs)/10),'% complete']);
				fprintf('\b|\n');
			end

		end
		fprintf('\n');
	else
		% Just report period 1 results (nothing stochastic here)
		startvals_orig = [start_ini,startvals(1:length(simvec)-1)];
		simseries = simvec;
		varnames = ['exst',varnames_store(1:length(simvec)-1)];
		
		warning('off','MATLAB:singularMatrix');
		warning('off','MATLAB:nearlySingularMatrix');
		[simseries, varnames] = mobj.computeSimulationMoments( ...
			simseries,varnames,[],startvals_orig, mitparams);
		warning('on','MATLAB:singularMatrix');
		warning('on','MATLAB:nearlySingularMatrix');
		
		nvars = length(varnames);
		tens_simseries(1:2,:,:) = [startvals; simseries];

	end
    varnames = varnames_store;
    nvars = length(varnames);

    % make HashMap with mapping of names to indices
    indexmap=java.util.HashMap;
    for i=1:nvars
        indexmap.put(varnames{i},i);
    end
%     varst=zeros(length(startpt_vec)-1,1);
%     for i=1:length(startpt_vec)-1
%         varst(i)=indexmap.get(mobj.En_names{i});
%     end

    %save(outfile,'tens_simseries','indexmap');

    simseries_median{s} = median(tens_simseries,3);
    simseries_mean{s} = mean(tens_simseries,3);
    simseries_std{s} = std(tens_simseries,[],3);

end

save(outfile,'simseries_mean','simseries_median','simseries_std','indexmap','NT_sim','N_shock','qC');

% if usejava('desktop')
%    plot_trans; 
% end


function [transmat,simnext,params,qC] = computeFirstPeriod(mobj,quantpoint,guessenstates,shockstruct,n)
	exst=shockstruct.exst;
	params = mobj.Params;
	paramDeviations = shockstruct.params;
	for ii=1:size(paramDeviations,1)
		params.(paramDeviations{ii,1}) = paramDeviations{ii,2};
	end
	
	exogenv = mobj.Exogenv;
	exogDeviations = shockstruct.pts;
	for ii=1:size(exogDeviations,1)
		colIdx = find( strcmp( exogDeviations{ii,1}, exogenv.exnames ) );
		exogenv.pts_perm(exst,colIdx) = exogDeviations{ii,2};
		exogenv.pts_all(exst,colIdx) = exogDeviations{ii,2};
	end

	% adjust initial guesses
	adjust = cell( mobj.NSOL, 1);
	adjust(:) = { @(x) x };
	% adjust initial guess for consumption
	adjust([6,11]) = { @(x) x - log(1 + params.demandShock * exp(x) ) };
	
	exnpt = mobj.Exogenv.exnpt;
	guesspoint = [exst, guessenstates];
	quantpoint(1) = exst;
	solguess=mobj.evaluatePol(guesspoint);
	for ii=1:mobj.NSOL
		solguess(ii) = adjust{ii}(solguess(ii));
	end
	transguess=mobj.evaluateTrans(guesspoint);
	valguess=mobj.evaluateVal(guesspoint);
	stateguess=guessenstates(2:3)';
	guess=[solguess;transguess([1,(exnpt+1):(4*exnpt)]);valguess(6);stateguess];
	objfun = @(x)computeMITShockState(quantpoint,x,mobj,params,exogenv);

	options=optimset('Display','off','TolX',1e-15,'TolFun',1e-12,...
			'MaxIter',100,'MaxFunEvals',100^2,'FinDiffType','central');
	[sol,~,exfl] = fsolve(objfun,guess,options);
	if exfl<1
		error('No Eqm Found in Run %d',n);
	end
	[~,simnext, transmat, qC]=objfun(sol);
	simnext = [simnext,nan(1,mobj.NCOND)];
end