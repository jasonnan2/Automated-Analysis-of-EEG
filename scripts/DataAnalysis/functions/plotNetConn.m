function plotNetConn(netConn, p_mat,netwrk,varargin )
% plotNetConn Visualizes a network connectivity matrix with significance overlay
%
%   plotNetConn(netConn, p_mat, netwrk)
%
%   This function plots a network connectivity matrix with overlaid markers
%   indicating statistically significant connections (p < 0.05).
%
%   Inputs:
%     netConn   - Square matrix of network connectivity values [N x N]
%     p_mat     - Corresponding p-value matrix [N x N]
%     netwrk    - Struct array with field 'name', providing network labels
%   Optional Parameters:
%     'MarkerSize' - Size of the significance markers (default: 12)
%     'LineWidth'  - Line width of the markers (default: 2)
%     'FontSize'   - Font size for axis tick labels (default: 12)
%     'Colormap'   - Colormap for the heatmap (default: parula)
%   Visualization Details:
%     - Uses imagesc to display the connectivity matrix
%     - Highlights significant connections with red '+' markers
%     - Y-axis labels are flipped to match matrix orientation
%     - Tick labels use network names from the netwrk input
%
%   Example:
%     plotNetConn(netConnMatrix, pValuesMatrix, networkStruct);

p = inputParser;
addParameter(p, 'MarkerSize', 12);
addParameter(p, 'LineWidth', 2);
addParameter(p, 'FontSize', 12);
addParameter(p, 'Colormap', parula);
parse(p, varargin{:});
opts = p.Results;



axis tight
p_mask=flipud(p_mat<0.05);
hold on
imagesc(flipud(netConn))
[xax,yax]=meshgrid(1:size(p_mask,2),1:size(p_mask,1));
plot(xax(p_mask), yax(p_mask),'r+','MarkerSize',12,'linewidth',2)
colorbar()
hold off

ax=gca;
ax.FontSize=12;
ax.YTick=[1:length(netwrk)];
ax.YTickLabel=flipud({netwrk.name}')';
ax.XTick=[1:length(netwrk)];
ax.XTickLabel={netwrk.name};
colorbar;
end