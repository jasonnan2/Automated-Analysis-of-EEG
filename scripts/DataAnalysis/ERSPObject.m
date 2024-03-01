classdef ERSPObject < DataAnalysis
    properties
        DATA
        erspResults
    end
    
    methods
        function obj = ERSPObject(DATA, info, baselineTime, timeRange)
            % Call superclass constructor
            obj = obj@DataAnalysis(info, baselineTime, timeRange);
            obj.DATA=DATA;
        end
        
        function newObj=createERP(obj,freqRanges)
            newObj=struct(); newObj.info=obj.info;
            N=length(obj.info.groupNames);
            properties=obj.info.variables;
            freq_list=fieldnames(freqRanges);
            for p=1:length(properties)
                for n=1:N
                    group=obj.info.groupNames{n};
                    scalpData=[];
                    for i = 1:length(freq_list)
                        freqidx = find(obj.info.freqAxis >= freqRanges.(freq_list{i})(1)& obj.info.freqAxis < freqRanges.(freq_list{i})(2));
                        scalpData(i,:,:,:)=nanmean(obj.DATA.(group).(properties{p})(freqidx,:,:,:),1);
                    end
                    newObj.DATA.(group).(properties{p})=scalpData;
                    newObj.DATA.(group).groupName=group;
                    newObj.DATA.(group).subList=obj.DATA.(group).subList;
                    if isfield(obj.DATA.(group),'missingSubs')
                        newObj.DATA.(group).missingSubs=obj.DATA.(group).missingSubs;
                    end
                end
            end
            newObj.info=rmfield(newObj.info,'freqAxis');
            newObj.info.freq_list=freq_list;
        end
        
        function plotERSP(obj,varargin)
            
            %   [...] = plotERSP(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
            %   parameters and their values.  Valid parameters are the following:
            %
            %        Parameter  Value
            %         'vars2plot'    cell array of variables in the object
            %                        to plot. Default is to plot all
            %                        variables
            %                  
            %         'chans2plot'   cell array of channels in the object
            %                        to plot. Default is to plot all
            %                        frequencies
            %                   
            %         'times2plot'   cell array of time ranges in the object
            %                        to plot. Default is to plot all
            %                        time ranges
            %         'combinations' 2D matrix defining which groups to run
            %                        comparisons on.
            %                        ex: [1,2] will do groups 1 vs 2
            %                            [1,2;1,3] is 1 vs 2 and 1 vs 3
            %         'cRange'        [-mx mx] range for colorlimit. default
            %                       is entire range of data points 
            N = length(obj.info.groupNames);
            pnames = {'vars2plot','chans2plot','times2plot','combinations','cRange'};
            dflts  = {obj.info.variables,{obj.info.chanlocs.labels}, fieldnames(obj.info.timeRange),nchoosek(1:N,2),nan};
            [vars2plot,chans2plot,times2plot,combinations,cRange] = parseArgs(pnames,dflts,varargin{:});
            [yTicks, yTickLabels] = logFreq2Ticks(obj.info.freqAxis);
            
            for p=1:length(vars2plot)
                property=vars2plot{p}; % get property name
                for comb = 1:size(combinations, 1)

                    group1=obj.info.groupNames{combinations(comb, 1)};
                    group2=obj.info.groupNames{combinations(comb, 2)};

                    for t=1:length(times2plot)
                        
                        timeName=times2plot{t};
                        timeIdx = obj.info.timeIDX.(timeName);
                        
                        for c=1:length(chans2plot)
                            figure
                            figCount=0;
                            chan=chans2plot{c};
                            cIDX=find(strcmp(chan,{obj.info.chanlocs.labels}));
                            
                            s1 = squeeze((obj.getGroupData(group1,property,'all',timeName,cIDX)));
                            s2 = squeeze((obj.getGroupData(group2,property,'all',timeName,cIDX)));
                            
                            if strcmp(obj.info.experimentalDesign,'paired')
                                [~,p]=ttest(s1,s2,'dim',3);
                            else
                                [~,p]=ttest2(s1,s2,'dim',3);
                            end
       
                            plotdata=squeeze(nanmean(s2,3)-nanmean(s1,3));
                            plotdata(p>0.05)=nan;
                            
                            if isnan(cRange)
                                mx = max(abs(prctile(nonzeros(plotdata),[10 90])));
                                mn = -mx;
                                cRange=[mn mx];
                            end
                            
                            figCount = figCount + 1;
                            subplot(t,3,figCount) % each row is a time, columns are significant mask, then individual plots
                            imagesc(obj.info.timeAxis(timeIdx(1):timeIdx(2)),(obj.info.freqAxis),squeeze(nanmean(s1,3)))
                            set(gca,'YDir','normal','FontSize',16,'YTick',yTicks,'YTickLabels',yTickLabels,'yLim',[2 30],'CLim',cRange);
                            ylabel('Freq (Hz)')
                            title(group1)
                            colorbar
                            
                            figCount = figCount + 1;
                            subplot(t,3,figCount) % each row is a time, columns are significant mask, then individual plots
                            imagesc(obj.info.timeAxis(timeIdx(1):timeIdx(2)),(obj.info.freqAxis),squeeze(nanmean(s2,3)))
                            set(gca,'YDir','normal','FontSize',16,'YTick',yTicks,'YTickLabels',yTickLabels,'yLim',[2 30],'CLim',cRange);
                            title(group2)
                            xlabel('Time (mS)')
                            colorbar
                            xt=xticks;
                            xtl=string(xticks);

                            figCount = figCount + 1;
                            subplot(t,3,figCount) % each row is a time, columns are significant mask, then individual plots
                            imagesc(obj.info.timeAxis(timeIdx(1):timeIdx(2)),(obj.info.freqAxis),plotdata)
                            [nr,nc] = size(plotdata);
                            pcolor([plotdata nan(nr,1); nan(1,nc+1)]);
                            shading flat;
                            colorbar
                            set(gca,'XTick',find(ismember(obj.info.timeAxis(timeIdx(1):timeIdx(2)),xt)), 'XTickLabels',xtl)
                            set(gca,'YDir','normal','FontSize',16,'YTick',yTicks,'YTickLabels',yTickLabels,'yLim',[2 30],'CLim',cRange);
                            title(group2+"-"+group1)

                            sgtitle(property+" "+chan)
                 
                        end
                    end
                end
            end
        end
        function plotERSPgroups(obj,varargin)
            
            %   [...] = plotERSP(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies additional
            %   parameters and their values.  Valid parameters are the following:
            %
            %        Parameter  Value
            %         'vars2plot'    cell array of variables in the object
            %                        to plot. Default is to plot all
            %                        variables
            %                  
            %         'chans2plot'   cell array of channels in the object
            %                        to plot. Default is to plot all
            %                        frequencies
            %                   
            %         'times2plot'   cell array of time ranges in the object
            %                        to plot. Default is to plot all
            %                        time ranges
            %         'groups2plot'   cell array of groupnames in the object
            %                        to plot. Default is to plot all
            %                        groups
            %         'cRange'        [-mx mx] range for colorlimit. default
            %                       is entire range of data points 
            N = length(obj.info.groupNames);
            pnames = {'vars2plot','chans2plot','times2plot','groups2plot','cRange'};
            dflts  = {obj.info.variables,{obj.info.chanlocs.labels}, fieldnames(obj.info.timeRange),obj.info.groupNames,nan};
            [vars2plot,chans2plot,times2plot,groups2plot,cRange] = parseArgs(pnames,dflts,varargin{:});
            [yTicks, yTickLabels] = logFreq2Ticks(obj.info.freqAxis);
            
            for p=1:length(vars2plot)
                property=vars2plot{p}; % get property name
                for t=1:length(times2plot)
                    timeName=times2plot{t};
                    timeIdx = obj.info.timeIDX.(timeName);
                    figure
                    figCount=0;
                    for c=1:length(chans2plot)
                        chan=chans2plot{c};
                        cIDX=find(strcmp(chan,{obj.info.chanlocs.labels}));                        
                        for n=1:length(groups2plot)
                            group=groups2plot{n};
                            s = squeeze((obj.getGroupData(group,property,'all',timeName,cIDX)));
                            figCount = figCount + 1;
                            subplot(length(chans2plot),3,figCount) % each row is a time, columns are significant mask, then individual plots
                            imagesc(obj.info.timeAxis(timeIdx(1):timeIdx(2)),(obj.info.freqAxis),squeeze(nanmean(s,3)))
                            set(gca,'YDir','normal','FontSize',16,'YTick',yTicks,'YTickLabels',yTickLabels,'yLim',[2 30],'CLim',cRange);
                            if n==1
                                %ylabel('Freq (Hz)')
                                ylabel(chan)
                            end
                            if c==1
                                title(group)
                            end
                            
                            colorbar
                        end
                        sgtitle(chan,'fontweight','bold','fontsize',16)
                    end
                end
            end 
        end
        %--------------END OF CLASS METHODS-----------------------%
    end
    
    
    
    %--------------END OF CLASS-----------------------%
end

%%
function [yTicks, yTickLabels] = logFreq2Ticks(freqAxis)
YTicks = [2 3 5 8 10 13 20 40];
YTicks(YTicks>max(freqAxis) | YTicks<min(freqAxis)) = [];
n = length(YTicks);
yTicks = zeros(n,1);
yTickLabels = cell(n,1);
for k=1:n
    [~, yTicks(k)] = min(abs(YTicks(k)-freqAxis));
    yTickLabels{k} = num2str(yTicks(k));
end
end
