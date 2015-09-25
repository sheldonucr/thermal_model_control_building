%% Intitalization
% Close all existing figures, plots, etc.

clear all;
close all;
clc;

% Variables
csv_location = 'data/useful data-with HVAC all the time with no other effects.csv';
startday = 26;
numdays = 160; % Number of days to test.

%% Load Data
fulldata = csvread(csv_location,1,1);
test_set = fulldata(startday:startday+numdays*24-1,1:end); % first 'numdays' amount of the fulldata set.
valid_set = fulldata(end-numdays*24+1:end, 1:end); %last numdays amount of the fulldata set.

[size, datapoints] = size(test_set); % Assuming test_set and valid_set are same size.

%% Generating scales
time = zeros(size,1);
for i = 0:size-1
    time(i+1,1) = mod(i,24)+1;
end


csvwrite_file = 'sweepresults3-useful data-with HVAC.csv';
%headers = {'order' 'test_rmse_1' 'test_rmse_2' 'test_rmse_3' 'test_rmse_4' 'test_rmse_5' 'valid_rmse_1' 'valid_rmse_2' 'valid_rmse_3' 'valid_rmse_4' 'valid_rmse_5'};
%csvwrite(csvwrite_file,headers,0,0);
%dlmwrite(csvwrite_file,headers);
% orders = 5:5:250;
orders = 1:1:100;

for i = 1:length(orders)
    order = orders(i)
    
    % Parse Testing data
    test_input = [time, test_set(1:end,1), test_set(1:end,12:26)];
    test_output = test_set(1:end,2:6);

    % Testing Phase
    SSM = my_subspace(test_input, test_output, time, order);
    testing_res = sim(SSM,test_input);
    
    % Parse Validation data
    valid_input = [time, valid_set(1:end,1), valid_set(1:end,12:26)];
    valid_output = valid_set(1:end,2:6);

    % Validation Phase
    valid_res = sim(SSM,valid_input);
    
    testing_absolute_rmse = sqrt(sum(abs(testing_res - test_output).^2)/size);
    valid_absolute_rmse = sqrt(sum(abs(valid_res - valid_output).^2)/size);
    
    %Save to CSV
    sweep_res = [order testing_absolute_rmse valid_absolute_rmse]
    %csvwrite(csvwrite_file,sweep_res,i,0);
    dlmwrite(csvwrite_file,sweep_res,'-append');
end
