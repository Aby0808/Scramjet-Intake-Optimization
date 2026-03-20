function [x_int, y_int] = FUN_find_intersection(x1, y1, m1, x2, y2, m2)
% findIntersection - Find the intersection point of two lines
%
% Syntax:
%   [x_int, y_int] = findIntersection(x1, y1, m1, x2, y2, m2)
%
% Inputs:
%   x1, y1 - Coordinates of a point on the first line
%   m1     - Slope of the first line
%   x2, y2 - Coordinates of a point on the second line
%   m2     - Slope of the second line
%
% Outputs:
%   x_int  - x-coordinate of the intersection point
%   y_int  - y-coordinate of the intersection point

% Check if the slopes are equal (parallel lines)
if abs(m1 - m2) < 1e-10
    error('The lines are parallel and do not intersect.');
end

% Line equations:
% Line 1: y = m1 * (x - x1) + y1
% Line 2: y = m2 * (x - x2) + y2

% Solve for the intersection point
x_int = (m1 * x1 - m2 * x2 + y2 - y1) / (m1 - m2);
y_int = m1 * (x_int - x1) + y1;
end
