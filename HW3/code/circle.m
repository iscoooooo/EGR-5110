function p = circle(x, y, radius, faceColor, edgeColor)
% Plots a circle with specified center, radius, and color.

% Calculate the position of the rectangle that will be drawn as a circle
% The position is defined as [xLeft yBottom width height]
position = [x-radius, y-radius, 2*radius, 2*radius];

% Create the circle with specified properties
rectangle('Position', position, ...
    'Curvature', [1 1], ...
    'FaceColor', faceColor, ...
    'EdgeColor', edgeColor);

% Create a proxy object for the circle in the legend
p = plot(NaN,NaN,edgeColor, 'Marker', 'o', 'MarkerFaceColor', faceColor, ...
    'LineStyle', 'none');

end