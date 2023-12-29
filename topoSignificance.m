%% Function to calculate which electrodes are significant btw two groups
function pvals=topoSignificance(s1,s2)

    % s1 and s2 are data matrices size C channels x N subjects
    % key is a two element cell array with markers for non-sig and sig
    % new_chan is the chan_locs
    % Perform t-tests
    
    s1=squeeze(s1);
    s2=squeeze(s2);
    for chan = 1:size(s1,1)
        % Get non-NaN indices
        [~,pvals(chan)] = ttest2(s1(chan,:), s2(chan,:));
        if isnan(pvals(chan))
            disp('asdfawrsasdbaweasbas')
            waitforbuttonpress
        end
    end

end
        
