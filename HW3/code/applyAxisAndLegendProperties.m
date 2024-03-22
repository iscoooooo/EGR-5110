function applyAxisAndLegendProperties(figHandle)
% Iterate over all axes in the figure and adjust properties
allAxes = findall(figHandle, 'type', 'axes');
for ax = allAxes'
    set(ax, ...
        'Color', 'k', ...
        'XColor', 'w', ...
        'YColor', 'w', ...
        'GridColor','w', ...
        'MinorGridColor','w');
    ax.Title.Color = 'w'; % Set title color to white
end

% Apply axes-specific commands
box(ax, 'on'),

% Adjust legend properties
legends = findobj(figHandle, 'Type', 'Legend');
for lg = legends'
    set(lg, 'Color', 'k', 'TextColor', 'w', 'EdgeColor', 'w');
end

end