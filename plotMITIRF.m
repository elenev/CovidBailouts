function handle=plotMITIRF(simdata,plot_shocks,plotspec,varargin)

varsel=plotspec.varsel;
titles=plotspec.titles;
irftransforms=plotspec.transforms;
pollabel=plotspec.pollabel;
irfplotidx=plotspec.irfplotidx;
coloridx=plotspec.coloridx;
fsize=plotspec.fsize;

simseries_mean=simdata.series;
indexmap=simdata.names;

N_per = size(simseries_mean{1},1);
N_plots = length(varsel);

printName=[];
if nargin>3
    printName=varargin{1};    
end

handle=figure;
if N_plots <= 3
    N_rows = 1;
    N_columns = N_plots;
elseif N_plots == 4
    N_rows = 2;
    N_columns = 2;
else
    N_columns = 3;
    N_rows = ceil( N_plots / N_columns );
end

% Obtain default color order
cls=get(gca,'colororder');

for i=1:N_plots
    varidx = indexmap.get(varsel{i});
    
    ssvals = arrayfun( @(j) simseries_mean{1}(:,varidx), plot_shocks, ...
        'UniformOutput', false); % not needed
    vals = arrayfun( @(j) simseries_mean{j}(:,varidx), plot_shocks, ...
        'UniformOutput', false );
    
    ssvals = [ssvals{:}];
    vals = [vals{:}];
    
    allvals = irftransforms{i}.fcn(vals,ssvals);
    subplot(N_rows,N_columns,i);
    hold on;
    for cci=1:length(irfplotidx)
        cc=irfplotidx(cci);
        ccc=coloridx(cci);
        plot(0:N_per-1,allvals(:,1+cc),'o-','LineWidth',1.5, ...
            'Color',cls(mod(ccc-1,7)+1,:));
    end
    hold off;
    title( titles{i} );
    xlim([0,N_per-1]);
    ylabel( irftransforms{i}.lab );
    
    if i==N_plots
        legend(pollabel,'location','northeast');
    end
    set(gca,'FontSize',fsize);
end


    if ~isempty(printName)
        printPanel( gcf, printName, [N_rows,N_columns] );
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
