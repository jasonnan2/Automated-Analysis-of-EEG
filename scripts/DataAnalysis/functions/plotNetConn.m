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
%     'Colormap'   - Colormap for the heatmap (default: red-blue-gradient)
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

cn=256;
bluePart = [linspace(0.4, 0.9, cn)', linspace(0.4, 0.9, cn)', ones(cn, 1)];
redPart  = [ones(cn, 1), linspace(0.9, 0.3, cn)', linspace(0.9, 0.3, cn)'];
defaultCmap = [bluePart;[1,1,1]; redPart];
% defaultCmap = 'parula';
addParameter(p, 'Colormap', defaultCmap);
parse(p, varargin{:});
opts = p.Results;

axis tight
p_mask=flipud(p_mat<0.05);
hold on
imagesc(flipud(netConn))

max_val = max(abs(netConn(:)));
if max_val > 0
    clim([-max_val, max_val]);
end
colormap(gca, opts.Colormap);
[xax,yax]=meshgrid(1:size(p_mask,2),1:size(p_mask,1));
plot(xax(p_mask), yax(p_mask),'k+','MarkerSize',12,'linewidth',2)
colorbar()
hold off

ax=gca;
ax.FontSize=12;
ax.YTick=[1:length(netwrk)];
ax.YTickLabel=flipud({netwrk.name}')';
ax.XTick=[1:length(netwrk)];
ax.XTickLabel={netwrk.name};

grid on;
ax.Layer = 'top';
ax.XTick = 0.5 : 1 : size(netConn, 2) + 0.5;
ax.YTick = 0.5 : 1 : size(netConn, 1) + 0.5;
ax.GridColor = [0, 0, 0];       % Black lines (use [0.5 0.5 0.5] for gray)
ax.GridAlpha = 0.3;             % Transparency (0 = invisible, 1 = solid)
ax.GridLineStyle = '-';         % Solid line style

colorbar;
end