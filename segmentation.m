clc;
clear;
List = dir("C:\Users\25297\Desktop\LiDAR\data\WREF\*.las");
k =length(List);
for i=1:1:k
%pc = LASread('C:\Users\25297\Desktop\LiDAR1\HARV_011.las');
path = List(i).folder;
name = List(i).name;
%path1 = path + ' \ '+ name;
path1 = strcat(path,'\');
path1 = strcat(path1,name);

path2 = strcat(path,'1\');
path2 = strcat(path2,name);

pc = LASread(path1);
cellSize = 0.8;
[models, refmat] = elevationModels([pc.record.x, pc.record.y, pc.record.z], ...
    pc.record.classification, ...
    'classTerrain', [2], ...%表示地形
    'classSurface', [1], ...%表示点
    'cellSize', cellSize, ...
    'interpolation', 'idw', ...
    'searchRadius', inf, ...
    'weightFunction', @(d) d^-3, ...
    'smoothingFilter', fspecial('gaussian', [3, 3], 0.8), ...
    'outputModels', {'terrain', 'surface', 'height'}, ...
    'fig', true, ...
    'verbose', true);

% export the Digital Terrain Model (DTM) to an ARC/INFO ASCII grid file
ASCwrite('zh_2014_a_dtm.asc', ...
    models.terrain.values, ...
    refmat, ...
    'precision', 2, ...
    'noData', -99999, ...
    'verbose', true);

% export the Digital Surface Model (DSM) to an ARC/INFO ASCII grid file
ASCwrite('zh_2014_a_dsm.asc', ...
    models.surface.values, ...
    refmat, ...
    'precision', 2, ...
    'noData', -99999, ...
    'verbose', true);

   



% export the Digital Height Model (DHM) to an ARC/INFO ASCII grid file
ASCwrite('zh_2014_a_dhm.asc', ...
    models.height.values, ...
    refmat, ...
    'precision', 2, ...
    'noData', -99999, ...
    'verbose', true);

[peaks_crh, ~] = canopyPeaks(models.height.values, ...
    refmat, ...
    'minTreeHeight', 2, ...
    'searchRadius', @(h) 0.28 * h^0.59, ...
    'fig', true, ...
    'verbose', true);

    [label_2d, colors] = treeWatershed(models.height.values, ...
    'markers', peaks_crh(:,1:2), ...
    'minHeight', 1, ...
    'removeBorder', true, ...
    'fig', true, ...
    'verbose', true);

metrics_2d = regionprops(label_2d, models.height.values, ...
    'Area', 'Centroid', 'MaxIntensity');

idxl_veg = ismember(pc.record.classification, [1]);

% convert map coordinates (x,y) to image coordinates (column, row)
RC = [pc.record.x - refmat(3,1), pc.record.y - refmat(3,2)] / refmat(1:2,:);
RC(:,1) = round(RC(:,1)); % row
RC(:,2) = round(RC(:,2)); % column
ind = sub2ind(size(label_2d), RC(:,1), RC(:,2));

% transfer the label
label_3d = label_2d(ind);
label_3d(~idxl_veg) = 0;
[label_3d(label_3d ~= 0), ~] = grp2idx(label_3d(label_3d ~= 0));

% transfer the color index
color_3d = colors(ind);
color_3d(~idxl_veg) = 1;



% define a colormap
cmap = [0, 0, 0;
    166,206,227;
    31,120,180;
    178,223,138;
    51,160,44;
    251,154,153;
    227,26,28;
    253,191,111;
    255,127,0;
    202,178,214;
    106,61,154;
    255,255,153;
    177,89,40] ./ 255;

%     figure('Color', [1,1,1])
%     scatter3(pc.record.x(idxl_veg), ...
%         pc.record.y(idxl_veg), ...
%         pc.record.z(idxl_veg), ...
%         6, ...
%         color_3d(idxl_veg), ...
%         'Marker', '.')
%     axis equal tight
%     colormap(cmap)
%     caxis([1, size(cmap,1)])
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
    
% compute metrics
[metrics_3d, fmt, idxl_scalar] = treeMetrics(label_3d, ...
    [pc.record.x, pc.record.y, pc.record.z], ...
    pc.record.intensity, ...
    pc.record.return_number, ...
    pc.record.number_of_returns, ...
    nan(length(pc.record.x), 3), ...
    models.terrain.values, ...
    refmat, ...
    'metrics', {'UUID', 'XPos', 'YPos', 'ZPos', 'H','BBOX2D', 'XCVH2D', 'YCVH2D', 'CVH2DArea', 'IQ50'}, ...
    'intensityScaling', true, ...
    'alphaMin', 1.5, ...
    'verbose', true);

% list field names
sfields = fieldnames(metrics_3d);

% IMPORTANT: adjust the path to the output CSV file (otherwise it will be created in the current folder)
fid = fopen('zh_2014_a_seg_metrics.csv', 'w+'); % open file
fprintf(fid, '%s\n', strjoin(sfields(idxl_scalar), ',')); % write header line
C = struct2cell(rmfield(metrics_3d, sfields(~idxl_scalar))); % convert structure to cell array
fprintf(fid, [strjoin(fmt(idxl_scalar), ','), '\n'], C{:}); % write cell array to CSV file
fclose(fid); % close file


%将彩色点云导出到las 1.4文件中
% duplicate the source file
r = pc;

% rescale the RGB colors to 16 bit range and add them to the point record
rgb = uint16(cmap(color_3d,:) * 65535);
r.record.red = rgb(:,1); 
r.record.green = rgb(:,2);
r.record.blue = rgb(:,3);

% add the "label" field to the point record (as an uint32 field)
r.record.label = uint32(label_3d);

r.record.blue = label_3d;

% add the "label" uint32 field metadata in the variable length records
% check the ASPRS LAS 1.4 specification for details about the meaning of the fields
% https://www.asprs.org/a/society/committees/standards/LAS_1_4_r13.pdf
vlr = struct;
vlr.value.reserved = 0;
vlr.value.data_type = 5;
vlr.value.options.no_data_bit = 0;
vlr.value.options.min_bit = 0;
vlr.value.options.max_bit = 0;
vlr.value.options.scale_bit = 0;
vlr.value.options.offset_bit = 0;
vlr.value.name = 'label';
vlr.value.unused = 0;
vlr.value.no_data = [0; 0; 0];
vlr.value.min = [0; 0; 0];
vlr.value.max = [0; 0; 0];
vlr.value.scale = [0; 0; 0];
vlr.value.offset = [0; 0; 0];
vlr.value.description = 'LABEL';

vlr.reserved = 43707;
vlr.user_id = 'LASF_Spec';
vlr.record_id = 4;
vlr.record_length_after_header = length(vlr.value) * 192;
vlr.description = 'Extra bytes';

% append the new VLR to the existing VLR
if isfield(r, 'variable_length_records')

    r.variable_length_records(length(r.variable_length_records)+1) = vlr;

else

    r.variable_length_records = vlr;

end

% if necessary, adapt the output record format to add the RGB channel
switch pc.header.point_data_format_id

    case 1 % 1 -> 3

        recordFormat = 3;

    case 4 % 4 -> 5

        recordFormat = 5;

    case 6 % 6 -> 7

        recordFormat = 7;

    case 9 % 9 -> 10

        recordFormat = 10;

    otherwise % 2,3,5,7,8,10

        recordFormat = pc.header.point_data_format_id;

end

% write the LAS 1.4 file
% IMPORTANT: adjust the path to the output LAS file
LASwrite(r, ...
    path2, ...
    'version', 12, ...
    'guid', lower(strcat(dec2hex(randi(16,32,1)-1)')), ...
    'systemID', 'SEGMENTATION', ...
    'recordFormat', recordFormat, ...
    'verbose', true);

% you can optionally read the exported file and check it has the
% RGB color and label records
% IMPORTANT: adjust the path to the input LAS file

r2 = LASread(path2);


%     idxl_sample = (label_3d == 6);
%     figure
%     scatter3(pc.record.x(idxl_sample), ...
%         pc.record.y(idxl_sample), ...
%         pc.record.z(idxl_sample), ...
%         12, ...
%         pc.record.intensity(idxl_sample), ...
%         'Marker', '.');
%     colorbar
%     axis equal tight
%     title('Return intensity')
%     xlabel('x')
%     ylabel('y')
%     ylabel('z')
end
    

