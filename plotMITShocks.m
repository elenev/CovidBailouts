function [handles,cellDictOut]=plotMITShocks(simdata,cvwelfare,plot_shocks,plotspec,varargin)

varsel=plotspec.varsel;
titles=plotspec.titles;
transforms=plotspec.transforms;
labelYaxis=plotspec.labelYaxis;
experiments=plotspec.experiments;
axlabel=plotspec.axlabel;
pollabel=plotspec.pollabel;
keylabel=plotspec.keylabel;
fsize=plotspec.fsize;
coloridx=plotspec.coloridx;

simseries_mean=simdata.series;
indexmap=simdata.names;
N_shocks=simdata.N_shocks;
N_policies=simdata.N_policies;

if nargin<4
    disp('Not enough inputs');
    return;
end

leave_off_noshock=true;
econs_output=[];
printName=[];
if nargin>4
    printName=varargin{1};
    if nargin>5
        econs_output=varargin{2};
        if nargin>6
            leave_off_noshock=varargin{3};
        end
    end
end

cellDictOut = [];
if ~isempty(econs_output)  
    cellDictOut = cell(500,2);
    %keylabel = {'donothing','ppp','mslp','ccf','combo','cbl'};
    headers = {'mit','mitpt'};
    econs = econs_output;
    dictDivider = {' - '};
    key_firstpart = strcat( headers(  kron(1:N_shocks ,ones(1,N_policies) ) ), ...
        dictDivider( ones(1,N_policies * N_shocks) ), ...
        econs(    kron(1:N_shocks ,ones(1,N_policies) ) ), ...
        dictDivider( ones(1,N_policies * N_shocks) ), ...
        keylabel( kron(ones(1,N_shocks) , 1:N_policies) ), ...
        dictDivider( ones(1,N_policies * N_shocks) ) );
end
    

divider = {': '};
tableCaptions = strcat(axlabel(kron(2:length(axlabel),ones(1,N_policies))), ...
	divider( ones(1,N_policies * (length(axlabel)-1 )) ), ...
	pollabel(kron(ones(1,length(axlabel)-1),1:N_policies)) );
tableCaptions = [axlabel(1), tableCaptions];

counter = 0;

N_panels = numel(varsel);
handles=cell(N_panels,1);

% Obtain default color order
cls=get(gca,'colororder');

for p=1:N_panels
    handles{p}=figure;
	N_plots = numel( varsel{p} );
	
	if N_plots==0
		continue;
	end
	
	if N_plots <= 3
		N_rows = 1;
		N_columns = N_plots;
	elseif N_plots == 4
		N_rows = 2;
		N_columns = 2;
	else
		N_rows = 2;
		N_columns = ceil( N_plots / N_rows );
	end
	
	tmpvarsel = varsel{p};
	tmptitles = titles{p};
	tmptransform = transforms{p};
	tmplabelYaxis = labelYaxis{p};
	tmpExperiments = experiments{p};
	
	graphValuesTable = zeros( N_plots, length(plot_shocks) );
	
	for i=1:N_plots
		currvarsel = tmpvarsel(i);
		
		if isempty( tmptransform{i} )
			continue;
		end
		
		if ~strcmp( currvarsel{:}, 'cvwelfare' )
			varidx = indexmap.get(tmpvarsel{i});

			ssvals = arrayfun( @(j) simseries_mean{1}(:,varidx), plot_shocks, ...
				'UniformOutput', false);
			vals = arrayfun( @(j) simseries_mean{j}(:,varidx), plot_shocks, ...
				'UniformOutput', false );

			ssvals = [ssvals{:}];
			vals = [vals{:}];
		else
			ssvals = nan(2, length(plot_shocks));
			vals = ssvals;
			vals(2,2:end) = cvwelfare;
		end
		
		allvals = tmptransform{i}.fcn(vals,ssvals);
		plotvals = allvals(2,:);
		graphValuesTable(i,:) = plotvals;
		
		plotvals_noshock = [plotvals(1), nan(1, N_policies-1)];
		plotvals_shocks = reshape( plotvals(2:end), N_policies, N_shocks )';

		if leave_off_noshock
			plotvals = plotvals_shocks;
			axlabel_plot = axlabel(2:end);
		else
			plotvals = [plotvals_noshock; plotvals_shocks];
			axlabel_plot = axlabel;
		end
		
		if strcmp(tmpExperiments,'onetime')
			plotvals = plotvals(1,:);
			axlabel_plot = axlabel_plot(1);
		elseif strcmp(tmpExperiments,'newnormal')
			plotvals = plotvals(2:end,:);
			axlabel_plot = axlabel_plot(2:end);
		end

		subplot(N_rows,N_columns,i);
			
		if size(plotvals,1)==1
			hold on;
			h{i} = cell(1,size(plotvals,2));
			for bb=1:size(plotvals,2)
				h{i}{bb}=bar(bb,plotvals(bb), ...
					'FaceColor', ...
					cls(mod(coloridx(bb)-1,7)+1,:));
			end
			hold off;
			set(gca,'xtick',1:bb);
			set(gca,'xticklabel',pollabel);
			set(gca,'xticklabelrotation',45);
		else
			h{i}=bar(plotvals, ...
					'FaceColor','flat');
				
				for jj=1:size(plotvals,2)
					for kk=1:size(plotvals,1)
						h{i}(jj).CData(kk,:) = cls(mod(coloridx(jj)-1,7)+1,:);
					end
				end
					
			set(gca,'xticklabel',axlabel_plot);
		end


		title( tmptitles{i} );
        set(gca,'FontSize',fsize);

		if tmplabelYaxis
			ylabel( tmptransform{i}.lab );
		end
		
		% Add to cellDict
        if ~isempty(econs_output)
            key_lastpart = strcat(currvarsel,dictDivider,{'mean'});
            combinedKeys = strcat( key_firstpart, ...
                key_lastpart( ones(1,N_policies * N_shocks) ) );
            outvals = plotvals';
            outvals = outvals(:);
            if strcmp(tmpExperiments,'both')
                cellDictOut(counter + (1:N_policies*N_shocks), 1) = combinedKeys';
                cellDictOut(counter + (1:N_policies*N_shocks), 2) = num2cell(outvals);
                counter = counter + N_policies*N_shocks;
			elseif strcmp(tmpExperiments,'onetime')
                cellDictOut(counter + (1:N_policies), 1) = combinedKeys(1:N_policies)';
                cellDictOut(counter + (1:N_policies), 2) = num2cell(outvals);
                counter = counter + N_policies;
			elseif strcmp(tmpExperiments,'newnormal')
				cellDictOut(counter + N_policies + (1:N_policies), 1) = combinedKeys(N_policies+1:2*N_policies)';
                cellDictOut(counter + N_policies + (1:N_policies), 2) = num2cell(outvals);
                counter = counter + N_policies;
            end
		end
        
    end
	
    if ~isempty(printName)
        printPanel( gcf, [printName,num2str(p)], [N_rows,N_columns] );
	end
	
	graphValuesTable = array2table(graphValuesTable);
	graphValuesTable.Properties.RowNames = tmptitles;
	
	writetable(cell2table(tableCaptions),[printName,'.xlsx'], ...
		'Sheet',['Panel ',num2str(p)],...
		'WriteRowNames',false,...
		'WriteVariableNames',false,...
		'Range','B1:Z1');
	pause(0.25);
	
	writetable(graphValuesTable,[printName,'.xlsx'], ...
		'Sheet',['Panel ',num2str(p)],...
		'WriteRowNames',true,...
		'WriteVariableNames',false,...
		'Range','A2:Z100');
end

if ~isempty(econs_output)
    cellDictOut=cellDictOut(1:counter,:);
	cellDictOut = cellDictOut(~cellfun(@isempty,cellDictOut(:,1)),:);
end



end

% Print
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