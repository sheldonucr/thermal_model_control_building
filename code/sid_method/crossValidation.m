%% Initialization
clear all;
close all;
clc;

% Variables
% csv_location = 'data/only consider ambient temperature (5zone model).csv';
% csv_location = 'data/consider all the factors (5 zone model).csv';
csv_location = 'data/Consider all facotrs(2 zone model).csv';

csvwrite_file = strcat('cv_results_',datestr(clock,30),'.csv');

order = 40;
%ordersweep = 5:5:90;
ordersweep = 4:2:90;

numSplits = 12; % Note: Cannot be zero.

%% Load Data
fulldata = csvread(csv_location,1,1);
[fd_size, fd_points] = size(fulldata);

time = zeros(fd_size,1);
for i = 0:fd_size-1
    time(i+1,1) = mod(i,24)+1;
end

time2 = (1:fd_size)';
fulldata = [fulldata time time2];
[fd_size, fd_points] = size(fulldata);
clear time time2;

% Cross valid sweep paraemeters calculation
numShift = floor((fd_size/24)/numSplits)*24;
mp = floor((fd_size/24)/2)*24; % midpoint

figure(1);
zRes = [];

%% Begin cross validation.
for i = 0:numSplits-1
    
    fprintf('Processing for iteration: %d\n',i);
    
    % Split the data.
    tmp = i*numShift;
    if tmp > mp % Training is split.
        valid = fulldata(floor(mod(tmp+1,mp)):mp+floor(mod(tmp,mp)),1:end);
        train = fulldata([mp+floor(mod(tmp+1,mp)):mp+floor(mod(tmp+1,mp))+(mp-floor(mod(tmp,mp))-1) 1:floor(mod(tmp,mp))],1:end);
        
        figure(1);
        hold on;
        xV = [floor(mod(tmp+1,mp)) mp+floor(mod(tmp,mp))];
        y = [i i];
        plot(xV, y, 'blue');
        
        xT1 = [mp+floor(mod(tmp,mp)) fd_size];
        xT2 = [1 floor(mod(tmp,mp))];
        plot(xT1,y,'red');
        plot(xT2,y,'red');
        clear xV y xT1 xT2;
        
        hold off;
        
    else % Validation is split.
        train = fulldata(floor(mod(tmp+1,mp)):mp+tmp,1:end);
        valid = fulldata([mp+tmp+1:mp+tmp+1+(mp-tmp-1) 1:tmp],1:end);
        
        figure(1);
        hold on;
        xT = [floor(mod(tmp+1,mp)) mp+tmp];
        y = [i i];
        plot(xT,y,'red');
        
        xV1 = [mp+tmp+1 fd_size];
        xV2 = [1 tmp];
        plot(xV1,y,'blue');
        plot(xV2,y,'blue');
        clear xV1 xV2 xT y;
        
        if i == 0
            legend('Training','Validation');
            xlabel('Time (hour)');
            ylabel('Iteration number');
            title('Sweeping data');
        end
        
        hold off;
    end
    
    % --- Begin work on splitted data here.
    
    % Parse Training data
    test_time = train(1:end,end-1);
    test_time2 = train(1:end,end);
    
    % test_input = [test_time, train(1:end,1), train(1:end,7:26)]; % 5zone,only considering ambient temperature
    % test_input = [test_time, train(1:end,1), train(1:end,7:56)]; % 5zone,conder all the factors
    test_input = [test_time, train(1:end,1), train(1:end,4:13)];% 2zone, Consider all facotrs
    
    test_output = train(1:end,2:6);
    
    % Parse Validation data
    valid_time = valid(1:end,end-1);
    valid_time2 = valid(1:end,end);
    
    % valid_input = [valid_time, valid(1:end,1), valid(1:end,7:26)];% 5zone,only considering ambient temperature
    % valid_input = [valid_time, valid(1:end,1), valid(1:end,7:56)];% 5zone,consider all the factors
    valid_input = [valid_time, valid(1:end,1), valid(1:end,4:13)];% 2zone, Consider all facotrs

    valid_output = valid(1:end,2:6);
    
    % Order Sweep
    for j=1:length(ordersweep)
        order = ordersweep(j);
        fprintf('Processing for order: %d\n',order);

        % Run simulation & Training Phase
        fprintf('Simulating for iteration,order: %d, %d\n',i,order);
        SSM = my_subspace(test_input, test_output, test_time2, order);
        testing_res = sim(SSM,test_input);

        % Validation Phase
        fprintf('Validating for iteration: %d\n',i);
        valid_res = sim(SSM,valid_input);

        % Calculate Error and save results.
        zzc = abs(testing_res - test_output)./test_output;
        testing_relative_rmse = sqrt(sum((zzc).^2)/mp);
        testing_max_relative_error = max(zzc);
        testing_absolute_rmse = sqrt(sum(abs(testing_res - test_output).^2)/mp);
        testing_max_absolute_error = max(abs(testing_res - test_output));

        zzc = abs(valid_res-valid_output)./valid_output;
        valid_relative_rmse = sqrt(sum((zzc).^2)/mp);
        valid_max_relative_error = max(zzc);
        valid_absolute_rmse = sqrt(sum(abs(valid_res - valid_output).^2)/mp);
        valid_max_absolute_error = max(abs(valid_res - valid_output));
        clear zzc;

        gTest = gfit2(test_output,testing_res,'all','v');
        gValid = gfit2(valid_output,valid_res,'all','v');

        zTmp = [i order testing_relative_rmse testing_max_relative_error testing_absolute_rmse testing_max_absolute_error valid_relative_rmse valid_max_relative_error valid_absolute_rmse valid_max_absolute_error gTest gValid];
        dlmwrite(csvwrite_file,zTmp,'-append');
        zRes = [zRes ; zTmp];
    
    end
end

%% Finalize
fprintf('Done executing script.\n');