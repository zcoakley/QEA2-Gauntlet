close all;

neatov2.connect("192.168.16.115");
% combined_map = [origin_0]; %; origin_135; one_1dot5_225; zero_1_60

[theta_clean, r_clean] = scan();
combined_map = scan_to_global(theta_clean, r_clean, 0, 0, 0);
pause(1)
plot_global(combined_map, "Global Plot");

[map, grad] = map_maker(combined_map, .2, 1.25, .125, true);
v = .1;
w_dot = .1;
lambda = .0005;
delta = 1;
lambda_i = lambda;
r = [0, 0];
theta = [0];
stop_threshold = 0.01;
total_drive_time = 0;
total_drive_dist = 0;

[xvals, yvals] = meshgrid(-3:0.05:3,-3:0.05:3);
[part_x, part_y] = grad(xvals, yvals);
plot_global(combined_map, "Global Plot");
quiv(xvals, yvals, part_x, part_y, 0, 0)
% neatov2.setFlatlandContours(xvals, yvals, zvals);
% 
% neatov2.setPositionAndOrientation(r(1), r(2), theta(1))
% neatov2.plotSim();

i = 1;
while (i == 1 || lambda_i*norm(grad(r(i, 1), r(i, 2))) > stop_threshold) %&& i < 20
    
    % if i == 6
    %     [theta_clean, r_clean] = scan();
    %     combined_map = scan_to_global(theta_clean, r_clean, r(i, 1), r(i, 2), theta(i));
    %     [map, grad] = map_maker(newPoints, xc, yc, radius, false);
    %     pause(1)
    %     figure
    %     plot_global(newPoints, "Global Plot without circle");
    %     hold on
    % 
    %     [xvals, yvals] = meshgrid(-3:0.05:3,-3:0.05:3);
    %     [part_x, part_y] = grad(xvals, yvals);
    %     quiv(xvals, yvals, part_x, part_y, r(i, 1), r(i, 2))
    % end
    
    if i == 12
        [theta_clean, r_clean] = scan();
        combined_map = scan_to_global(theta_clean, r_clean, r(i, 1), r(i, 2), theta(i));
        possible_circ = possible_circles(combined_map, .1, .05, 4, 10);
        [xc, yc, radius] = findClusterSmallestArea(possible_circ);
        newPoints = remove_circ(combined_map, xc, yc, radius+.05);
        [map, grad] = map_maker(newPoints, xc, yc, radius, true);
        pause(1)
        figure
        plot_global(newPoints, "Global Plot without circle");
        hold on
        circle(xc, yc, radius)
         
        [xvals, yvals] = meshgrid(-3:0.05:3,-3:0.05:3);
        [part_x, part_y] = grad(xvals, yvals);
        quiv(xvals, yvals, part_x, part_y, r(i, 1), r(i, 2))    

    end
    [partial_x, partial_y] = grad(r(i, 1), r(i, 2));
    
    r(i + 1, :) = r(i, :) - lambda_i.*[partial_x, partial_y];

    theta(i + 1) = pi + atan2(partial_y, partial_x);

    [time, Vl, Vr] = calculate_turn(theta(i), theta(i + 1), w_dot);
    drive(time, Vl, Vr)
    total_drive_time = total_drive_time + time;
    total_drive_dist = total_drive_dist + time * (Vl + Vr) ./ 2;
    
    [time, Vl, Vr] = calculate_straight(r(i, :), r(i + 1, :), v);
    drive(time, Vl, Vr)
    total_drive_time = total_drive_time + time;
    total_drive_dist = total_drive_dist + time * (Vl + Vr) ./ 2;

    i = i + 1;
end

total_drive_time
total_drive_dist

function [time, Vl, Vr] = calculate_turn(theta_i, theta_f, w)

    % angle = mod(theta_f - theta_i + pi, 2*pi) - pi;
    angle = theta_i - theta_f;
    if angle < 0
        w = -w;
    end
    time = angle/w;
    Vl = .245*.5*w;
    Vr = -.245*.5*w;

    % time = time + .

end

function [time, Vl, Vr] = calculate_straight(ri, rf, v)
    
    distance = abs(norm(rf-ri));
    time = abs(distance/v);
    Vl = v;
    Vr = v;
    
end

function drive(drivetime, Vl, Vr)
% neatov2.plotSim();
tic %%start your timer in Matlab
t=toc; %initiate t as the time since you started
while t<drivetime
    neatov2.setVelocities(Vl, Vr);
    pause(.01); %you can add a short delay so we aren't constantly changing the velocities. 
    % neatov2.plotSim();
    t=toc; %t update t
end
neatov2.setVelocities(0, 0);

end

function [fun, grad] = map_maker(x_y_mat, circ_x, circ_y, rad, add_ball)

    function val = map(x, y) 
        val = 0;

        for i = 1:size(x_y_mat, 1)
            val = val - log(sqrt(((x-x_y_mat(i, 1)).^2 + (y-x_y_mat(i, 2)).^2)));
        end

        if add_ball
            for theta = 0:0.1:2*pi
                a = circ_x+rad*cos(theta);
                b = circ_y+rad*sin(theta);
                val = val + log(sqrt((x-a).^2 + (y-b).^2));
            end
        end

    end

    function [partial_x, partial_y] = gradients(x, y) 
        partial_x = 0;
        partial_y = 0;

        for i = 1:size(x_y_mat, 1)
            partial_x = partial_x - (x - x_y_mat(i, 1))./((x-x_y_mat(i, 1)).^2+(y-x_y_mat(i, 2)).^2);
            partial_y = partial_y - (y - x_y_mat(i, 2))./((x-x_y_mat(i, 1)).^2+(y-x_y_mat(i, 2)).^2);
        end


        if add_ball
            for theta = 0:0.1:2*pi
                a = circ_x+rad*cos(theta);
                b = circ_y+rad*sin(theta);
                partial_x = partial_x + (x - a)./((x-a).^2+(y-b).^2);
                partial_y = partial_y + (y - b)./((x-a).^2+(y-b).^2);
            end
        end

    end

    fun = @map;
    grad = @gradients;

end

function quiv(X, Y, U, V, startx, starty)

    % Calculate the magnitude of each vector
    magnitudes = sqrt(U.^2 + V.^2);
    
    % Avoid division by zero
    % magnitudes(magnitudes == 0) = 1;
    
    % Normalize the vectors
    U_normalized = U ./ magnitudes;
    V_normalized = V ./ magnitudes;
    
    % Plot the vectors
    % quiver(X, Y, U_normalized, V_normalized, '');
    streamline(X, Y, - U_normalized, - V_normalized, startx, starty); % The '0' argument prevents scaling
    axis equal; % Keeps the aspect ratio so that circles look like circles
    title('Unit Vectors Showing Direction Only');
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
end

function circ_removed = remove_circ(data, cent_x, cent_y, rad)
    x = data(:, 1);
    y = data(:, 2);
    x_minus = (x - cent_x).^2;
    y_minus = (y - cent_y).^2;
    distances = sqrt(x_minus + y_minus);
    circ_removed = data(distances > rad, :);
end

function [theta_clean, r_clean] = scan()
    pause(3)
    sensors = neatov2.receive();
    pause(3)
    sensors = neatov2.receive();
    angles = sensors.thetasInRadians;
    ranges = sensors.ranges;
    theta_clean = angles(ranges ~= 0);
    r_clean = ranges(ranges ~= 0);
end


function global_coords = scan_to_global(angles, ranges, x, y, phi)
    d = .1;
    global_coords = [];
    for i = 1:length(angles)
        new_x = ranges(i)*cos(angles(i)+phi)-d*cos(phi) + x;
        new_y = ranges(i)*sin(angles(i)+phi)-d*sin(phi) + y;
        global_coords = [global_coords; new_x, new_y];
    end
end

function plot_global(global_coords, title_str)
    figure
    plot(global_coords(:, 1), global_coords(:, 2), '.')
    xlabel("x-axis")
    ylabel("y-axis")
    title(title_str)
end

function [xc, yc, r] = findClusterSmallestArea(data)
    % This function finds the three largest clusters using k-means clustering and returns the points in the cluster that occupies the smallest area,
    % as well as the indices of these points in the original dataset.
    % Inputs:
    %   data - A matrix where each row represents a data point and each column represents a feature.
    %   k - The number of clusters to identify (should be at least 3).
    %
    % Outputs:
    %   smallestAreaCluster - The index of the cluster with the smallest area.
    %   smallestArea - The area of this cluster.
    %   clusterPoints - A matrix containing the coordinates of the points in the smallest area cluster.
    %   pointIndices - The indices of the points in the original data matrix that belong to the smallest area cluster.

    % Apply k-means clustering to the data
    [idx, ~] = kmeans(data(:, 1:2), fix(size(data,1)/3));

    % Calculate the size of each cluster
    clusterCounts = accumarray(idx, 1);

    % Find indices of the three largest clusters
    [~, sortedIndices] = sort(clusterCounts, 'descend');
    topClusters = sortedIndices(1:fix(size(data,1)/9), end);

    % Initialize variables to find the cluster with the smallest area
    smallestArea = inf;
    smallestAreaCluster = NaN;

    % Loop through the top three clusters to find the one with the smallest area
    for i = 1:length(topClusters)
        clusterIndex = topClusters(i);
        pointsInCluster = data(idx == clusterIndex, :);
        indicesInCluster = find(idx == clusterIndex);

        % Calculate the bounding box area of this cluster
        width = max(pointsInCluster(:,1)) - min(pointsInCluster(:,1));
        height = max(pointsInCluster(:,2)) - min(pointsInCluster(:,2));
        area = width * height;

        % Update if this cluster has the smallest area
        if area < smallestArea
            smallestArea = area;
            smallestAreaCluster = clusterIndex;
            clusterPoints = pointsInCluster;
            pointIndices = indicesInCluster;
        end
    end
    means = mean(data(pointIndices, :), 1);
    xc = means(1);
    yc = means(2);
    r = means(3);

end

function [circles] = possible_circles(points, radius, radius_error, num_points, points_within)

    circles = [];
    for i=1:(length(points)-num_points)
        points_can = points(i:i+num_points, :);
        [xc_can, yc_can, r_can] = circ_finder(points_can);
        if abs(r_can-radius) > radius_error
            continue
        end
        circle(xc_can, yc_can, r_can);
        circles = [circles; xc_can, yc_can, r_can];
        % if enough_within_dist(points_can, radius, xc_can, yc_can, points_within, .1) == false
        %     continue
        % end
    end

end

function circle(x,y,r)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit);
end

function [xc, yc, r] = circ_finder(points)

    A = [2.*points(:, 1), 2.*points(:, 2), -1*ones(length(points),1)];
    b = points(:, 1).^2 + points(:, 2).^2;
    w = A\b;
    xc = w(1);
    yc = w(2);
    r = sqrt(xc^2 + yc^2 - w(3));

end