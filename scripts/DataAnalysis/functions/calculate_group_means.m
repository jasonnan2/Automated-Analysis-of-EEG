function [tbldata,newtable] = calculate_group_means(tbldata, groups, groupnames)
% Calculate the mean for specified groups of columns in the table data.
    % Inputs:
    %   - tbldata: The input table data.
    %   - groups: A cell array of cell arrays, where each inner cell array contains column keywords.
    %   - groupnames: A cell array of strings, corresponding to the names of the calculated means.
    %
    % Example usage:
    %   tbldata = calculate_group_means(tbldata, groups, groupnames);
    %
    % For each group, calculate the mean and add a new column with the group name.
    newtable=tbldata(:,1:2); % getting subID and group labels
    for g = 1:length(groups)
        cols = contains(tbldata.Properties.VariableNames, groups{g});
        tbldata.(groupnames{g}) = mean(tbldata{:, cols}, 2);
        newtable.(groupnames{g})= mean(tbldata{:, cols}, 2);
    end
    
end