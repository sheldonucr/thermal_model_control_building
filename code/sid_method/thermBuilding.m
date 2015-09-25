%% Intitalization
% Close all existing figures, plots, etc.

clear all;
close all;
clc;

% Variables
% csv_location = 'data/only consider ambient temperature (5zone model).csv';
csv_location = 'data/consider all the factors (5 zone model).csv';
% csv_location = 'data/Consider all facotrs(2 zone model).csv';
startday = 1; % set to 26 for first 2 csv locations.
numdays = 60; % Number of days to test.
order = 40;
rooms = 5;


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

titles = ['SPACE 1: Testing Phase - Simulation Results',
          'SPACE 2: Testing Phase - Simulation Results',
          'SPACE 3: Testing Phase - Simulation Results',
          'SPACE 4: Testing Phase - Simulation Results',
          'SPACE 5: Testing Phase - Simulation Results' ];

titles2 = ['SPACE 1: Validation Phase - Simulation Results',
          'SPACE 2: Validation Phase - Simulation Results',
          'SPACE 3: Validation Phase - Simulation Results',
          'SPACE 4: Validation Phase - Simulation Results',
          'SPACE 5: Validation Phase - Simulation Results' ];
      
      
%% Parse Testing data
% test_input = [time, test_set(1:end,1), test_set(1:end,7:26)]; % 5zone,only considering ambient temperature
test_input = [time, test_set(1:end,1), test_set(1:end,7:56)]; % 5zone,conder all the factors
% test_input = [time, test_set(1:end,1), test_set(1:end,4:13)];% 2zone, Consider all facotrs
% test_input = [time, test_set(1:end,1)];
test_output = test_set(1:end,2:6);

%% Calculate SSM 
SSM = my_subspace(test_input, test_output, time, order);

%% Testing Phase
testing_res = sim(SSM,test_input);
figure(1);
%sim(SSM,test_input);
for i=1:rooms
   subplot(rooms,1,i);
   plot(1:size,testing_res(1:end,i), 'red');
   hold on;
   plot(1:size,test_output(1:end,i),':blue');
   %plot(1:672,out_temp+(25-median(out_temp)),'green');
   legend('Estimated','EnergyPlus');
   xlabel('Time (Hours)');
   ylabel('Temperature (C)');
   title(titles(i,1:end));
   hold off;
end

%% Parse Validation data
% valid_input = [time, valid_set(1:end,1), valid_set(1:end,7:26)];% 5zone,only considering ambient temperature
valid_input = [time, valid_set(1:end,1), valid_set(1:end,7:56)];% 5zone,consider all the factors
% valid_input = [time, valid_set(1:end,1), valid_set(1:end,4:13)];% 2zone, Consider all facotrs
% valid_input = [time, valid_set(1:end,1)];
valid_output = valid_set(1:end,2:6);

%% Validation Phase
valid_res = sim(SSM,valid_input);
figure(2);
%sim(SSM,valid_input);
for i=1:rooms
   subplot(rooms,1,i);
   plot(1:size,valid_res(1:end,i),'red');
   hold on;
   plot(1:size,valid_output(1:end,i),':blue');
   %plot(1:672,out_temp,'green');
   legend('Estimated','EnergyPlus');
   xlabel('Time (Hours)');
   ylabel('Temperature (C)');
   title(titles2(i,1:end));
   hold off;
end

%% Calculating Errors

%Relative
% zza = (max(test_output) - min(test_output));
% zzb = abs(testing_res - test_output);
% for i = 1:length(zzb)
%    zzc(i,1:rooms) = zzb(i,1:end) ./ zza;
% end
zzc = abs(testing_res - test_output)./test_output;
testing_relative_rmse = sqrt(sum((zzc).^2)/size)
testing_max_relative_error = max(zzc)

testing_absolute_rmse = sqrt(sum(abs(testing_res - test_output).^2)/size)
testing_max_absolute_error = max(abs(testing_res - test_output))

% zza = (max(valid_output) - min(valid_output));
% zzb = abs(valid_res - valid_output);
% for i = 1:length(zzb)
%    zzc(i,1:rooms) = zzb(i,1:end) ./ zza;
% end
zzc = abs(valid_res-valid_output)./valid_output;
valid_relative_rmse = sqrt(sum((zzc).^2)/size)
valid_max_relative_error = max(zzc)

% testing_absolute_error = sum(abs(testing_res - test_output))/size
% valid_absolute_error = sum(abs(valid_res - valid_output))/size
valid_absolute_rmse = sqrt(sum(abs(valid_res - valid_output).^2)/size)
valid_max_absolute_error = max(abs(valid_res - valid_output))
 
% testing_percent_absolute_error = testing_absolute_rmse ./ (max(testing_res) - min(testing_res))
% valid_percent_absolute_error = valid_absolute_rmse ./ (max(valid_res) - min(valid_res))
% testing_percent_relative_error = testing_relative_rmse ./ (max(testing_res) - min(testing_res))
% valid_percent_relative_error = valid_relative_rmse ./ (max(testing_res) - min(testing_res))

%% Debug
% test_test = abs(testing_res - test_output);
% test_valid = abs(valid_res - valid_output);
% for i = 1:5
%     test_max_index(i) = find(testing_max_absolute_error(i) - test_test(1:end,i) <= 0.001,1);
%     valid_max_index(i) = find(valid_max_absolute_error(i) - test_valid(1:end,i) <= 0.001,1);
% end
% %Print
% test_max_index
% valid_max_index