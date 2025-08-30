clc
clear
close all

%% AMS Generator general
% Method: collect the year component only from the dates and then find the
% annual maximum flow from that. this will rectify the problem of incomplete
% years

File_Path="E:\Mayank\4_New_Zealand\Complete_year\cleaned_9150_original_cleaned.xlsx";
[~, ~, DATA] = xlsread(File_Path);

 dates =DATA(2:end, 1); % Skip the first row (heading)
 dateValues = datetime(dates, 'InputFormat', 'dd-MM-yyyy'); % Adjust the date format
 uniqueYears = unique(year(dateValues));
 streamflow=(DATA(2:end,3));
% Convert streamflow cell array elements to double
streamflow_double = zeros(size(streamflow));
for i = 1:numel(streamflow)
    % Convert each element to double
    if isnumeric(streamflow{i})
        streamflow_double(i) = streamflow{i}; % If already numeric, just copy
    else
        streamflow_double(i) = str2double(streamflow{i}); % Convert string to double
    end
end
 streamflow=streamflow_double;


% Initialize arrays to store peak flow values and their durations
peakFlowValues = zeros(size(uniqueYears));
peakFlowDurations = zeros(size(uniqueYears));
floodVolumes = zeros(size(uniqueYears));
for i = 1:length(uniqueYears)
    yearIndex = year(dateValues) == uniqueYears(i);
    yearStreamflow = streamflow(yearIndex);
    
    % Find peak flow and its index
    [peakFlow, peakIndex] = max(yearStreamflow);
    
    % Find adjacent left local minima
    leftMinIndex = peakIndex - 1;
    % while leftMinIndex > 1 && yearStreamflow(leftMinIndex) < yearStreamflow(leftMinIndex+1)
    %     leftMinIndex = leftMinIndex - 1;
    % end
      while leftMinIndex > 1 && yearStreamflow(leftMinIndex)> 0
        leftMinIndex = leftMinIndex - 1;
    end
    % % Find adjacent right local minima
    rightMinIndex = peakIndex + 1;
    % while rightMinIndex < length(yearStreamflow) && yearStreamflow(rightMinIndex) < yearStreamflow(rightMinIndex-1)
    %     rightMinIndex = rightMinIndex + 1;
    % end
    while rightMinIndex < length(yearStreamflow) && yearStreamflow(rightMinIndex)>0
        rightMinIndex = rightMinIndex + 1;
    end
    if(rightMinIndex>365)
        rightMinIndex=365;
    end
    % Find the start and end dates of the peak flow event
    startDate = dateValues(yearIndex);
    startDate = startDate(leftMinIndex+1);
    endDate = dateValues(yearIndex);
    endDate = endDate(rightMinIndex-1);
    
    % Calculate duration of peak flow event
    duration = (endDate - startDate);
    duration1=(duration)/days(1);
    % Store peak flow value and its duration
    peakFlowValues(i) = peakFlow;
    peakFlowDurations(i) = days(duration)+1; % Convert duration to days
     % Calculate flood volume
     if(leftMinIndex==0)
         leftMinIndex=1;
     end
     floodVolume = sum(yearStreamflow(leftMinIndex:rightMinIndex)) - (yearStreamflow(rightMinIndex) - yearStreamflow(leftMinIndex)) *         (duration1 / 2);
    
    % Store flood volume
    floodVolumes(i) = floodVolume;
end
Flood_characterstics = [uniqueYears, peakFlowValues, peakFlowDurations, floodVolumes];

