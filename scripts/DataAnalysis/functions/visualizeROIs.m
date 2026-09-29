function visualizeROIs(hmFile, rois2plot)
% VISUALIZEROIS Renders a 3D visualization of specified Regions of Interest (ROIs).
%
%   VISUALIZEROIS(hmFile, rois2plot) loads a head model file and highlights
%   the designated ROI indices on the cortical surface.
%
%   Inputs:
%       hmFile    - String or char array. Path to the head model file.
%       rois2plot - Numeric vector. Indices of the ROIs to highlight (1 to 68).
%
%   Example:
%       visualizeROIs('headmodel.mat', [5, 12, 34]);
% Plot ROI 
hm = headModel.loadFromFile(hmFile);
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';                            
plotdata = zeros(68,1);
plotdata(rois2plot) = 1;
X=T'*plotdata;
plot68roi(hm,X , 1,{''})
end
