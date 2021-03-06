% Store all stats in CSV file

% Clear
clc; clear;

% Economies
bench_econ = 'xi88midxigrid';
other_econs = {'xi88pandemic'};

		   
econs = [{bench_econ},other_econs];

simres_filename = 'simres';
simres_path = ['Results\',simres_filename];
fid=0;%fid = fopen([simres_path,'.csv'],'w');
%fprintf(fid,'%s,%s\n','key','value');
sampleIdentifiers = {'u','e','r','c','l'};
irfSampleIdentifiers = {'u','r','c'};

% Set column headers
colnames={'mean','std','corrG','corrG_1','corrOm','corrOm_1', ...
	'corrY','corrY_1','corrY_gr','corrY_gr_1','AC', ...
	'min','p01','p05','p10','p25','p50','p75','p90','p95','p99','max'};

unitless = {'corrG','corrG_1','corrOm','corrOm_1', ...
	'corrY','corrY_1','corrY_gr','corrY_gr_1','AC'};
isUnitless = ismember(colnames,unitless);

bench_simres = load(['sim_res_',bench_econ,'.mat']);

deffmt = '%s,%0.3f\n';

%% Set parameter formats
paramNames = fieldnames(bench_simres.params);
paramFormats = cell( size( paramNames) );
paramFactors = cell( size( paramNames ) );
paramFormats(:) = {'%0.2f'};
paramFactors(:) = {@(x)x};

match = @(target) cellfun( @(x) contains(x,target), paramNames );
isPct = match('PB');
isPct = isPct | strcmp('PS',paramNames);
isPct = isPct | strcmp('gammaB',paramNames);
isPct = isPct | strcmp('gammaS',paramNames);
isPct = isPct | strcmp('sig2_A',paramNames);
isPct = isPct | strcmp('sigma_eps',paramNames);
isPct = isPct | strcmp('deltaK',paramNames);
isPct = isPct | strcmp('gammaG',paramNames);
isPct = isPct | strcmp('T',paramNames);
isPct = isPct | strcmp('tau0',paramNames);
isPct = isPct | strcmp('tauPi',paramNames);
isPct = isPct | strcmp('tauD',paramNames);
isPct = isPct | strcmp('etaI',paramNames);
isPct = isPct | strcmp('zetaI',paramNames);

isSquare = match('sig2_A');
isSquare = isSquare | match('sig2_om');

paramFormats( isPct ) = {'%0.3f'};
paramFactors( isPct & ~isSquare) = {@(x)100*x};
paramFactors( isPct & isSquare) = {@(x)100*sqrt(x)};
paramFactors( ~isPct & isSquare) = {@(x)sqrt(x)};



cellDict = cell( 1e8, 2 );
ii = 1;
for econCounter = 1:numel(econs)
	econ = econs{econCounter};
	if econCounter > 1
		simres = load(['sim_res_',econ,'.mat']);
	else
		simres = bench_simres;
		
		% Parameters
		for paramCounter=1:numel(paramNames)
			param = paramNames{paramCounter};
			vecValue = simres.params.(param);
			for paramElementCounter = 1:length(vecValue)
				paramId = sprintf('param - %s - %d', param,paramElementCounter);
				value = paramFactors{paramCounter}( vecValue( paramElementCounter ) );
				cellDict(ii,:) = writePair(paramId, value, fid, ['%s,',paramFormats{paramCounter},'\n']) ;
				ii = ii + 1;
			end
		end
	end
	
	[~,varIdx,varIdx_bench] = intersect( simres.dispnames, ...
			bench_simres.dispnames, 'stable' );
	
	% Set factors and formats (see function at the bottom of the file)
	[factors, formats] = defineFormats( simres );
	
	% Simulation Results
	for sampleCounter = 1:numel(simres.tabout_exog)
		tab = simres.tabout_exog{sampleCounter};
		bench_tab = bench_simres.tabout_exog{sampleCounter};
		
		sampleId = sampleIdentifiers{sampleCounter};
		
		%[have_irf, irf_idx] = ismember(sampleId, irfSampleIdentifiers);
        have_irf=false;
		if have_irf
			irfseries = irfres.simseries_mean{irf_idx};
		end
		
		for intersectVarCounter = 1:length(varIdx)
			varCounter = varIdx( intersectVarCounter );
			benchVarCounter = varIdx_bench( intersectVarCounter );
			var = simres.dispnames{varCounter};
			%var = strrep( var, '_', '' );
			%var = strrep( var, '0', 'zero');
			%var = strrep( var, '1', 'one');
			%var = strrep( var, '2', 'sq');
			fmtStr = formats{varCounter};
			factor = factors{varCounter};
			for statCounter = 1:size(tab,2)
				%resId = sprintf('sim - %s - %s - %s - %s',econ,...
				%	sampleId, var, colnames{statCounter} );
				resId = ['sim - ',econ,' - ',sampleId,' - ',var,' - ',colnames{statCounter}];
				if isUnitless(statCounter)
					tmpFactor = 1;
				else
					tmpFactor = factor;
				end
				value = tmpFactor * tab(varCounter,statCounter);
				cellDict(ii,:) = writePair(resId, value, fid, ['%s,',fmtStr,'\n']) ;
				
				%compId = sprintf('comp - %s - %s - %s - %s',econ,...
				%	sampleId, var, colnames{statCounter} );
				compId = ['comp - ',econ,' - ',sampleId,' - ',var,' - ',colnames{statCounter}];
				%if factor == 100
				%	value = tab(varCounter,statCounter) - ...
				%		bench_tab(varCounter,statCounter);
				%else
					value = tab(varCounter,statCounter) ./ ...
						bench_tab(benchVarCounter,statCounter) - 1;
				%end
				value = 100 * value;
				cellDict(ii+1,:) = writePair(compId, value, fid, deffmt) ;
				ii = ii + 2;
			end
			
			% IRFs
			if have_irf && ~isempty( irfres.indexmap.get(var) )
				irf_values = irfseries( [1,2], irfres.indexmap.get(var) );
				%levelId = sprintf('irf - %s - %s - %s - %s',econ,...
				%	sampleId, var, 'mean' );
				levelId = ['irf - ',econ,' - ',sampleId,' - ',var,' - mean'];
				value = factor * irf_values(end);
				cellDict(ii,:) = writePair(levelId, value, fid, ['%s,',fmtStr,'\n']) ;
				
				%changeId = sprintf('irf - %s - %s - %s - %s',econ,...
				%	sampleId, var, 'change' );
				changeId = ['irf - ',econ,' - ',sampleId,' - ',var,' - change'];
				if factor == 100
					value = diff( irf_values );
				else
					value = irf_values(end) ./ irf_values(1) - 1;
					value = 100*value;
				end
				cellDict(ii+1,:) = writePair(changeId, value, fid, deffmt) ;
				ii = ii + 2;
			end
			
		end
    end
    
    isFLmodel=simres.params.CBS;
    
	% Manual
    if ~isFLmodel
        % IQR of fracS changes
        cellDict(ii,:) = writePair( ...
            ['sim - ',econ,' - u - d_fracS - IQR'], ...
            -100*diff( prctile( simres.simseries( :, simres.indexmap.get('d_fracS') ), [75,25] ) ), ...
            fid, '%s,%0.2f\n' );
        ii = ii+1;
    end

	% Fiscal policy sensitivities
	logY = simres.simseries( :, simres.indexmap.get('logY') );
	regY = @(X) corr(X,logY) * std(X) / std(logY);
	for var = {'lgammaG_byY','lT_byY','lLtax_byY'}
		var = var{:};
		cellDict(ii,:) = writePair( ...
			['sim - ',econ,' - u - ',var,' - regY'], ...
			regY( simres.simseries( :, simres.indexmap.get(var) ) ), ...
			fid, deffmt );
		ii = ii+1;
	end

	% IQR of productivity XS stds
	sdIQR = std(diff(norminv([0.75,0.25], ...
		simres.simseries(:,simres.indexmap.get('mu_om')), ...
		simres.simseries(:,simres.indexmap.get('Om')).^0.5),[],2));

	cellDict(ii,:) = writePair( ...
		['sim - ',econ,' - u - prodIQR - std'], ...
		100*sdIQR, ...
		fid, deffmt );
	ii = ii+1;
	
    if ~isFLmodel
        % Franchise value
        franch = mean( simres.simseries(:, simres.indexmap.get('VI') ) ) ./ ...
            mean( simres.simseries(:, simres.indexmap.get('WI') ) ) - 1;
        cellDict(ii,:) = writePair( ...
			['sim - ',econ,' - u - franch - nmean'], ...
            100*franch, ...
            fid, deffmt );
        ii = ii+1;
        
        % Accounting ROE
        accROE = mean( simres.simseries(:, simres.indexmap.get('Imasset') ) ) .* ...
            mean( simres.simseries(:, simres.indexmap.get('Lspr') ) ) ./ ...
            mean( simres.simseries(:, simres.indexmap.get('WI') ) );
        cellDict(ii,:) = writePair( ...
            ['sim - ',econ,' - u - accROE - nmean'], ...
            100*accROE, ...
            fid, deffmt );
        ii = ii+1;
        
        if econCounter==1
            bench_franch = franch;
            bench_accROE = accROE;
        end
        
        cellDict(ii,:) = writePair( ...
            ['comp - ',econ,' - u - franch - nmean'], ...
            100*(franch./bench_franch-1), ...
            fid, '%s,%0.2f\n' );
        ii = ii+1;
        cellDict(ii,:) = writePair( ...
            ['comp - ',econ,' - u - accROE - nmean'], ...
            100*(accROE./bench_accROE-1), ...
            fid, deffmt );
        ii = ii+1;
	end
% 	
% 	% MIT shocks
% 	try
% 		mitfilename = ['MIT_res_',econ,'.mat'];
% 		mitres = load(mitfilename);
% 		
% 		
% 	catch err
% 		if strcmp(err.identifier,'MATLAB:load:couldNotReadFile')
% 			warning(['File ',mitfilename,' doesn''t exist. Skipping...']);
% 		end
% 	end
% 	
% 	% MIT shocks combined with transitions from bench
% 	try
% 		mitPTfilename = ['MITPT_res_',econ,'.mat'];
% 		mitPTres = load(mitPTfilename);
% 		
% 		addlCellDict = processMITshock( mitPTres.simseries_mean, 
% 	catch err
% 		if strcmp(err.identifier,'MATLAB:load:couldNotReadFile')
% 			warning(['File ',mitPTfilename,' doesn''t exist. Skipping...']);
% 		end
% 	end
% 	
	
end

try
	% Process welfare
	welfare = readtable( ['Results\welfare_',bench_econ,'.xls'],'Sheet','NumClust 10');
	for econCounter=1:size(welfare,1)
		wId = sprintf('comp - %s - u - cvwelfare - mean', welfare.Row{econCounter} );
		value = welfare.Total(econCounter);
		cellDict(ii,:) = writePair(wId, 100*value, fid, deffmt );
		ii = ii+1;
	end
catch err
	if strcmp(err.identifier,'MATLAB:spreadsheet:book:fileOpen')
		warning('Welfare XLS file not found. Skipping CEV Welfare measures');
	end
end

% Process benchmark errors
errprct = readtable( ['Results\errstats_res_',bench_econ,'.xls'] );
labels = errprct.Properties.VariableNames;
for eqnCounter = 1:size(errprct,1)
	for prctCounter = 1:size(errprct,2)
		cellDict(ii,:) = writePair( ...
			sprintf('err - bench - u - eq%02d - %s', eqnCounter, labels{prctCounter} ), ...
			errprct{ eqnCounter, labels{prctCounter} }, ...
			fid, '%s,%0.4f\n' );
		ii = ii + 1;
	end
end

% Wrap up
cellDict = cellDict(1:ii-1,:);
mapDict = containers.Map( cellDict(:,1), cellDict(:,2) );
%fclose(fid);
save([simres_path,'.mat'],'mapDict');

% Create latex db
%path_to_datatooltk = '"C:\Program Files (x86)\datatooltk\bin\datatooltk"';
%command = sprintf('%s --csv %s.csv --csvheader --sort key --output %s.dbtex',...
%	path_to_datatooltk,simres_path,simres_path);
%system(command);

function keyvalcell = writePair( key, val, varargin)
	keyvalcell = { key, val };
	%if nargin > 2
	%	fprintf(fid, fmt, key, val);
	%end
end

% Set variable formats
function [factors, formats, isPct] = defineFormats( simres )
	formats = cell( size( simres.dispnames ) );
	factors = cell( size( simres.dispnames ) );
	formats(:) = {'%0.3f'};
	factors(:) = {1};

	match = @(target) cellfun( @(x) contains(x,target), simres.dispnames );
	isPct = match('_gr');
	isPct = isPct | match('log');
	isPct = isPct | match('lvg');
	isPct = isPct | match('rate');
	isPct = isPct | match('LGD');
	isPct = isPct | match('rB');
	isPct = isPct | match('rD');
	isPct = isPct | match('rT');
	isPct = isPct | match('ret');
	isPct = isPct | match('frac');
	isPct = isPct | match('Lspr');
	isPct = isPct | match('LsprT');
	isPct = isPct | match('Tspr');
	isPct = isPct | match('_byY');
	isPct = isPct | match('X_byKB');
	isPct = isPct | match('exp');
	isPct = isPct | match('bind');
	isPct = isPct | strcmp('brupt', simres.dispnames);
	isPct = isPct | strcmp('franch', simres.dispnames);
	isPct = isPct | strcmp('Imasset', simres.dispnames);
	isPct = isPct | strcmp('expRIE', simres.dispnames);
	isPct = isPct | strcmp('wacc', simres.dispnames);
	isPct = isPct | strcmp('WI', simres.dispnames);


	formats( isPct ) = {'%0.3f'};
	factors( isPct ) = {100};
end
