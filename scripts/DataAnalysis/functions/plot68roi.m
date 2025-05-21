function fig = plot68roi(hm, X, colorlim,labelnames)

%X(X<0) = 2;

[hm.fvLeft,hm.fvRight, hm.leftH, hm.rightH] = geometricTools.splitBrainHemispheres(hm.cortex);

fig = figure('Color',[1 1 1]);
fig.Position(3:4) = [660   800];

if length(colorlim)==2
    mn = colorlim(1);
    mx = colorlim(2);
    
elseif colorlim
    top = max(abs([max(X(:)) min(X(:))]));
    mx = top;
    mn = -top;
else
    mx = prctile(abs(X(:)), 100);
    mn = -mx; %prctile(abs(X(:)), 10);
    if mn == mx && mn == 0
        mx = prctile(abs(X(:)), 100);
        mn = -mx;
        if mn == mx && mn == 0
            mn = -1;
            mx = 1;
        end
    end
end

if mx<=mn
    mx=mn+0.0001;
end

AxesH = axes('Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off', ...
    'YLimMode', 'manual', 'YLim',  [0, 1], ...
    'XTick',    [],       'YTick', [], ...
    'NextPlot', 'add', ...
    'HitTest',  'off');


const = 2;
ncols = size(X,2)+const;


%%


for i = 1:ncols-const

ax = axes('Position', [ (1/ncols)*i*1-0.1    0.800    0.3000    0.1500]); % subplot(511);
patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X(:,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([-90 90])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');

nsteps=2049;
cmap = colormap(ax,bipolar(nsteps,0.99));%[customCMap1;customCMap2]; %colormap(bipolar(512, 0.99)); %cmap =
colorsteps=linspace(mn,mx,nsteps);
cmap(find(abs(colorsteps)==min(abs(colorsteps))),:)=[1 1 1]; % set closest point to 0 to white


colormap(ax, cmap);

if i == ncols-const
    cb=colorbar;

    cb.Position = cb.Position + 1e-10;
end


% title('top view');

ax = axes('Position', [ (1/ncols)*i*1-0.08    0.715    0.2600    0.0500]); %subplot(512);
patch('vertices',hm.fvLeft.vertices,'faces',hm.fvLeft.faces,'FaceVertexCData',X(hm.leftH,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([-180 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');

colormap(ax, cmap);

% 
ax = axes('Position',[ (1/ncols)*i*1-0.08    0.590    0.2600    0.0500]); %subplot(513);
patch('vertices',hm.fvRight.vertices,'faces',hm.fvRight.faces,'FaceVertexCData',X(hm.rightH,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);

set(ax,'Clim',[mn mx]);
view([0 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');

colormap(ax, cmap);


%title('Right')

ax = axes('Position',[ (1/ncols)*i*1-0.03    0.4596   0.1400    0.05]); %subplot(514);
patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X(:,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);

set(ax,'Clim',[mn mx]);
view([90 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');

colormap(ax, cmap);


% title('front view');

ax = axes('Position',[ (1/ncols)*i*1-0.035    0.314    0.1400    0.05]); %subplot(515);
patch('vertices',hm.cortex.vertices,'faces',hm.cortex.faces,'FaceVertexCData',X(:,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([-90 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');

colormap(ax, cmap);

ax = axes('Position',[ (1/ncols)*i*1-0.1    0.2114    0.2800    0.05]); %subplot(515);
patch('vertices',hm.fvLeft.vertices,'faces',hm.fvLeft.faces,'FaceVertexCData',X(hm.leftH,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([0 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');

colormap(ax, cmap);
% 
ax = axes('Position',[ (1/ncols)*i*1-0.1    0.0914    0.2800    0.05]); %subplot(515);
patch('vertices',hm.fvRight.vertices,'faces',hm.fvRight.faces,'FaceVertexCData',X(hm.rightH,i),...
    'FaceColor','interp','FaceLighting','phong','LineStyle','none','FaceAlpha',1,'SpecularColorReflectance',0,...
    'SpecularExponent',25,'SpecularStrength',0.25,'Parent',ax);
set(ax,'Clim',[mn mx]);
view([-180 0])
axis(ax,'equal','vis3d','tight');
camlight(0,180)
camlight(0,0)
axis(ax,'on');

colormap(ax, cmap);

set(findall(fig,'type','axes'),'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[],'Box','off', 'visible','off')
% fig.Position(3:4) = [405   755];
text((1/ncols)*i*1,0.97,labelnames{i},'units','normalized','Parent', AxesH);


end


text(0.06,0.9,'top','units','normalized','Parent', AxesH);
text(0.06,0.75,'left','units','normalized','Parent', AxesH);
text(0.06,0.61,'right','units','normalized','Parent', AxesH);
text(0.06,0.49,'front','units','normalized','Parent', AxesH);
text(0.06,0.37,'back','units','normalized','Parent', AxesH);
text(0.03,0.24,'medial-left','units','normalized','Parent', AxesH);
text(0.03,0.11,'medial-right','units','normalized','Parent', AxesH);

