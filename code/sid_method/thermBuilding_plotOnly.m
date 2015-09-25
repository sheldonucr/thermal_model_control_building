clc;
close all;

%% Fixing plot to output only first space.
rooms=1;
%% Generating scales
time = zeros(size,1);
for i = 0:size-1
    time(i+1,1) = mod(i,24)+1;
end

titles = ['Zone 1: Testing Phase - Simulation Results',
          'Zone 2: Testing Phase - Simulation Results',
          'Zone 3: Testing Phase - Simulation Results',
          'Zone 4: Testing Phase - Simulation Results',
          'Zone 5: Testing Phase - Simulation Results' ];

titles2 = ['Zone 1: Validation Phase - Simulation Results',
          'Zone 2: Validation Phase - Simulation Results',
          'Zone 3: Validation Phase - Simulation Results',
          'Zone 4: Validation Phase - Simulation Results',
          'Zone 5: Validation Phase - Simulation Results' ];

titles3 = ['Zone 1: EnergyPlus Results',
          'Zone 2: EnergyPlus Results',
          'Zone 3: EnergyPlus Results',
          'Zone 4: EnergyPlus Results',
          'Zone 5: EnergyPlus Results' ];

%% Trimming Set data to plot only 500.
plotsize = 500;
testing_res = testing_res(1:plotsize,1:end);
test_output = test_output(1:plotsize,1:end);
valid_res = valid_res(1:plotsize,1:end);
valid_output = valid_output(1:plotsize,1:end);
      
%% Testing Phase
%testing_res = sim(SSM,test_input);
figure(1);
%sim(SSM,test_input);
for i=1:rooms
   subplot(rooms,1,i);
   plot(1:plotsize,testing_res(1:end,i), 'red');
   hold on;
   plot(1:plotsize,test_output(1:end,i),':blue');
   %plot(1:672,out_temp+(25-median(out_temp)),'green');
   legend('Estimated','EnergyPlus');
   xlabel('Time (Hours)');
   ylabel('Temperature (C)');
   title(titles(i,1:end));
   hold off;
end

set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% Validation Phase
%valid_res = sim(SSM,valid_input);
figure(2);
%sim(SSM,valid_input);
for i=1:rooms
   subplot(rooms,1,i);
   plot(1:plotsize,valid_res(1:end,i),'red');
   hold on;
   plot(1:plotsize,valid_output(1:end,i),':blue');
   %plot(1:672,out_temp,'green');
   legend('Estimated','EnergyPlus');
   xlabel('Time (Hours)');
   ylabel('Temperature (C)');
   title(titles2(i,1:end));
   hold off;
end

set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% Plot for energyPlus section
% figure(3);
% for i=1:rooms
%    subplot(rooms,1,i);
%    plot(1:plotsize,test_output(1:end,i),'blue');
%    %plot(1:672,out_temp+(25-median(out_temp)),'green');
%    legend('EnergyPlus');
%    xlabel('Time (Hours)');
%    ylabel('Temperature (C)');
%    title(titles3(i,1:end));
%    hold off;
% end