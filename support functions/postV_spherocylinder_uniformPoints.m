% This function is a supporting function for post-segmentation using
% LCuts to construct a spherocylinder model.
%
% For more information, please refer to the following papers:
% Zhang M, Zhang J, Wang Y, Wang J, Achimovich AM, Acton ST and 
% Gahlmann A. Non-invasive single-cell morphometry in living bacterial 
% biofilms. Nature communications, 11(1), pp.1-13, 2020.
%
% This code is written by Mingxing Zhang, University of Virginia.

function [model_points] = postV_spherocylinder_uniformPoints(l,r,grid_size)
%% generate grid points with l and r
x = -r:grid_size:r;
y = -r:grid_size:r;
z_c = -l/2:grid_size:l/2; 
z_up = l/2:grid_size:l/2+r;
z_down = -l/2-r:grid_size:-l/2;
radius = r;

%% generate the cylinder first
[points_rows_c,points_columns_c,points_pages_c] = meshgrid(x,y,z_c);

mesh_size_cylinder = size(points_rows_c);
xyz_page_cylinder = [];

for i = 1:mesh_size_cylinder(3)
    x_coord = points_rows_c(:,:,i);
    y_coord = points_columns_c(:,:,i);
    xy = [];
    logic_Points_In_circle = (points_rows_c(:,:,i)).^2 + (points_columns_c(:,:,i)).^2 <= radius.^2;
    for j = 1:length(logic_Points_In_circle(:,1))
        x_row = x_coord(j,logic_Points_In_circle(j,:));
        y_row = y_coord(j,logic_Points_In_circle(j,:));
        xy_temp = [x_row' y_row'];
        if ~isempty(xy_temp)
            xy = [xy; xy_temp];
        end
    end
    if ~isempty(xy)
        z = ones(length(xy),1) * points_pages_c(1,1,i);
        xyz = [xy z];
        xyz_page_cylinder = [xyz_page_cylinder; xyz];
    end
end
clear xy_temp xy xyz z

%% Create the sphere centered at (0 0 l/2) and (0 0 -l/2)
[points_rows_up,points_columns_up,points_pages_up] = meshgrid(x,y,z_up);
[points_rows_down,points_columns_down,points_pages_down] = meshgrid(x,y,z_down);

centerX = 0;
centerY = 0;
centerZ_up = l/2;
centerZ_down = -1*l/2;

sphereVoxels_up = (points_rows_up - centerY).^2 ...
    + (points_columns_up - centerX).^2 + (points_pages_up - centerZ_up).^2 <= radius.^2;

sphereVoxels_down = (points_rows_down - centerY).^2 ...
    + (points_columns_down - centerX).^2 + (points_pages_down - centerZ_down).^2 <= radius.^2;

mesh_size_up = size(sphereVoxels_up);
xyz_page_up = [];

logic_Points_In_circle = (points_rows_up).^2 + (points_columns_up).^2 + (points_pages_up - centerZ_up).^2 <= radius.^2;
for i = 1:mesh_size_up(3)
    x_coord = points_rows_up(:,:,i);
    y_coord = points_columns_up(:,:,i);
    xy = [];
    logic_Points_In_circle_temp = logic_Points_In_circle(:,:,i);
    for j = 1:length(logic_Points_In_circle_temp(:,1))
        x_row = x_coord(j,logic_Points_In_circle_temp(j,:));
        y_row = y_coord(j,logic_Points_In_circle_temp(j,:));
        xy_temp = [x_row' y_row'];
        if ~isempty(xy_temp)
            xy = [xy; xy_temp];
        end
    end
    if ~isempty(xy)
        z = ones(length(xy(:,1)),1) * points_pages_up(1,1,i);
        xyz = [xy z];
        xyz_page_up = [xyz_page_up; xyz];
    end
end
clear xy_temp xy xyz z
mesh_size_down = size(sphereVoxels_down);
xyz_page_down = [];
logic_Points_In_circle = (points_rows_down).^2 + (points_columns_down).^2 + (points_pages_down - centerZ_down).^2 <= radius.^2;
for i = 1:mesh_size_down(3)
    x_coord = points_rows_down(:,:,i);
    y_coord = points_columns_down(:,:,i);
    xy = [];
    logic_Points_In_circle_temp = logic_Points_In_circle(:,:,i);
    for j = 1:length(logic_Points_In_circle_temp(:,1))
        x_row = x_coord(j,logic_Points_In_circle_temp(j,:));
        y_row = y_coord(j,logic_Points_In_circle_temp(j,:));
        xy_temp = [x_row' y_row'];
        if ~isempty(xy_temp)
            xy = [xy; xy_temp];
        end
    end
    if ~isempty(xy)
        z = ones(length(xy(:,1)),1) * points_pages_down(1,1,i);
        xyz = [xy z];
        xyz_page_down = [xyz_page_down; xyz];
    end
end
clear xy_temp xy xyz z
model_points = [xyz_page_up; xyz_page_cylinder; xyz_page_down];

end