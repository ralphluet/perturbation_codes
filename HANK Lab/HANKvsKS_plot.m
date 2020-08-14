aggrshock='uncertainty'
%% Load KS IRFs
load(['../Krusell Smith Lab/IRF_SS_BASELINE_HANC_UNC.mat'])

%% Load HANK IRFs

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Y", "C", "I"];
opts.VariableTypes = ["double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("export_df.csv", opts);

% Convert to output type
IRF_Y_HANK = tbl.Y;
IRF_C_HANK = tbl.C;
IRF_I_HANK = tbl.I;


% Clear temporary variables
clear opts


%%
mpar.maxlag=40

figurename=['IRF_Y_HANKvsKS_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 800 800])
plot(1:mpar.maxlag-1,IRF_Y,'b-','LineWidth',3.5)
hold on
plot(1:mpar.maxlag-1,IRF_Y_HANK(1:mpar.maxlag-1),'r--','LineWidth',3.5)
legend('Krusell-Smith economy','HANK with illiquid capital','Location','SouthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
ylim([-0.4 0.101])
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black','HandleVisibility','off')
set(gca, 'FontName','arial','FontSize',40);
set(gca,'Position',[0.2 0.15 0.75 0.8])
printpdf(gcf,['../latex/' figurename])

figurename=['IRF_C_HANKvsKS_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 800 800])
plot(1:mpar.maxlag-1,IRF_C,'b-','LineWidth',3.5)
hold on
plot(1:mpar.maxlag-1,IRF_C_HANK(1:mpar.maxlag-1),'r--','LineWidth',3.5)
legend('Krusell-Smith economy','HANK with illiquid capital','Location','SouthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black','HandleVisibility','off')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])

figurename=['IRF_I_HANKvsKS_' aggrshock];
figure('Name',figurename,'NumberTitle','off','Position',[1 1 800 800])
plot(1:mpar.maxlag-1,IRF_I,'b-','LineWidth',3.5)
hold on
plot(1:mpar.maxlag-1,IRF_I_HANK(1:mpar.maxlag-1),'r--','LineWidth',3.5)
legend('Krusell-Smith economy','HANK with illiquid capital','Location','SouthEast')
ylabel('Percent','Interpreter','none','FontName','arial','FontSize',40)
xlabel('Quarter','Interpreter','none','FontName','arial','FontSize',40)
hold on
plot(0:mpar.maxlag-1,zeros(1,mpar.maxlag),'--black','HandleVisibility','off')
set(gca, 'FontName','arial','FontSize',40); set(gca,'Position',[0.2 0.15 0.75 0.8]) 
printpdf(gcf,['../latex/' figurename])
