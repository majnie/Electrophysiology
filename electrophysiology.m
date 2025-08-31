function [Table,fig1,fig2] = electrophysiology(directory)

% The purpose of this function is to batch process data collected from
% electrophysiology experiments according to the time at which drugs were
% applied to the hippocampal slices in each experiment respectively.

% INPUT = an excel file cotaining drug application times per experiment
% OUTPUT = a table containing the filenames of all the data sets in the directory, 
% and the peak amplitude of the first fEPSP recorded per application of each respective drug.

% Assumptions of the this function include: the first
% stimulation was applied at the same time for all electrophysiology
% experiments (at 10ms); assume that the time interval between data collection is the same for all the
% experiments (9.76500000000000e-05); assume that the number of data points collected per
% recording is the same (2048).

% First call upon the data set by its directory path
directory = dir('notebook.xlsx');
file = fullfile(directory.folder, directory.name);
notebook = readtable(file);

% Extract all the wcp files from the directory 
fileList = dir('*.wcp');

% Loop through each wcp file and import using the import_wcp function
for i = 1:length(fileList)
    filename = fullfile(fileList(i).name); % Within the loop, filename is assigned the name of the current file being processed 
    % in each iteration to permit saving all the files into the ephys_data cell. 
    ephys_data(i) = import_wcp(filename);
end

% Ensure that the order of the files in the excel file and ephys_data array
% correspond to one another.
notebook_sorted = sortrows(notebook,1);
ephys_sorted = orderfields(ephys_data);

% Add a column to the notebook pertaining to the last recording of each
% experiment for ease of handling later.
for t = 1:length(ephys_sorted)
    End_times(t) = length(ephys_sorted(t).rec_index);
end

notebook_sorted.End = End_times';

% Extract one of the data files
sample_wcp = ephys_data(1); % isolate the details pertaining to the example experiment
sample_data = ephys_data(1).S{1}; 
name = sample_wcp.file_name; % find the times of drug application pertaining to the sample data by file name 
find_times = strcmp(notebook_sorted.Filename,name);
sample_times = notebook_sorted(find_times,:); % extract the corresponding drug application times

% Create a time axis for the sample data
data_size = size(sample_data,1);
time = (0:sample_wcp.t_interval:(data_size-1)*sample_wcp.t_interval)*1000;

% Extract the time interval in which the first fEPSP was produced.
index = time>11 & time<59;
% here is is assumed that the 1st stimulus occurs at 10ms and the second at
% 60ms. It is also assumed that an artifact directly follows stimulus
% application.
window = find(index);

% Remove experimental artifacts by making the baseline closer to zero and smooth the data to remove noise
num_conditions = width(sample_times{1,2:end});

smoothfactor = 11;

for s = 1:width(sample_data)
    baseline_period = 1:sum(time<11);
    baseline_average = mean(sample_data(baseline_period,s));
    new_sample_data = sample_data - baseline_average;
    column_data = new_sample_data(:,s);
    smooth_data(:,s) = smooth(column_data,smoothfactor); % Create a new array containing only the smoothed data from the example data set.
end 

% Create a figure displaying the example traces for each of the conditions
% from the smoothed data.
figure;
fig1 = gcf;

variablenames = {'Baseline 1st fEPSP peak', 'CADO 1st fEPSP peak', 'CADO+DPCPX 1st fEPSP peak', 'CADO+DPCPX+NBQX 1st fEPSP peak' };
% Important for the example traces to be clearly identified from one
% another.

for a = 1:num_conditions
    first_fEPSP{a} = smooth_data(window,sample_times{1,a+1}-1); % have to omit the first column 
    % as it contains the filenames and have to select the recording immediately before the next 
    % drug application to obtained the data pertaining to each of the drug equilibriums. 
    
    subplot(2,4,a);
    plot(time(index),first_fEPSP{a}); % only want to plot the data pertaining to the first fEPSP 
    % produced in response to each of the 4 treatments. 

    xlim([min(time(index)) max(time(index))]);
    ylim([-0.6 0.3]); % remove the artifacts 
    xlabel('time (ms)');
    ylabel(sample_wcp.channel_info{1}.unit);

    title(variablenames(a));

    hold on
end 


% Create a subplot in this multipanel figure displaying a time course plot 
subplot(2,4,5:8);
hold on

for o = 1:size(sample_data,2)
    j = smooth_data(window,o); % only want to plot time vs peak amplitude of the 1st fEPSP per recording
    n(o)= min(j);
end 

scatter(1:size(n,2),n(1,:));
xlabel('time (min)');
ylabel(sample_wcp.channel_info{1}.unit);
title('Time course plot');

% Extract the peak amplitude of the first fEPSP produced for each conditon
% per each experiment.
for m = 1:size(ephys_sorted,2) % all files in the directory must be processed

    for g = 1:width(ephys_sorted(m).S{1})
        b_period = 1:sum(time<11);
        b_average = mean(ephys_sorted(m).S{1}(b_period,g));
        data = ephys_sorted(m).S{1} - b_average;
        column = data(:,g);
        smoothed(:,g) = smooth(column,smoothfactor); % ensure noise is removed from all the data sets that are processed 
    end 

    for q = 1:num_conditions
        fEPSP = smoothed(window,notebook_sorted{1,q+1}-1); % only extract the infromation pertaining to the first fEPSP 
        % and to the points at which the drug action is at equilibrium.
        peak_fEPSP(m,q) = min(fEPSP); % find the peak of the fEPSP using the minimum function. 
    end

filenames{m} = ephys_sorted(m).file_name;

end 

% Create a table which will serve as the main output of the function
T = array2table(peak_fEPSP,"VariableNames",variablenames);
Table = addvars(T, filenames', 'Before', 1,'NewVariableNames','Filenames'); 

% Find the average 1st fEPSP peak per condition 
for z = 1:num_conditions
    average(1,z) = mean(peak_fEPSP(:,z)); % Use a for loop to find the averages of all the peaks stored in the table.
    STD(1,z) = std(peak_fEPSP(:,z)); % Standard devation can be used as error bars.
end 

% Plot bar graph to demonstrate the average peak fEPSP per contition
figure;
fig2 = gcf;
b = bar(average); 
b.FaceColor = [0.3010 0.7450 0.9330];
hold on
xticklabels({'Baseline', 'CADO', 'CADO+DPCPX', 'CADO+DPCPX+NBQX'});
xlabel('Condition','FontSize',14);
ylabel('fEPSP Amplitude (mV)','FontSize',14);
title('Average 1st fEPSP peak per condition','FontSize',17);
errorbar(average, STD, 'k.', 'LineWidth', 1); % add standrad deviation as error bar

% Carry out statistical analysis to identify significant differences
for w = 1:num_conditions
    normality(w) = lillietest(peak_fEPSP(:,w)); % use lillietest to identify normaility of the data
end

if sum(normality) > 0 % this means that if any data sets are not normally distributed, a non-parametric test will be used to analyse the data 
    [p, tbl, stats] = kruskalwallis(peak_fEPSP,variablenames);
elseif sum(normality) == 0 % this means that all the data sets are noramlly disributed and a parametric test can be used 
    [p, tbl, stats] = anova1(peak_fEPSP,variablenames);
end 

% A post-hoc multiple comparisons test will be carried out if the results
% present a signficant difference greater than 0.05
c = multcompare(stats,"Alpha",0.05,"CriticalValueType","bonferroni","Display","off");

end