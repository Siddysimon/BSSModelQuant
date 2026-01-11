headList = dir('BB_*'); 
subjects = append('TIME',{'020','041','045', '052', '057', '072', '081','082','085'});
plotNames = ["AbsCondChange_fit__vs_field_change", "ROI2elect_dist_fit_vs_field_change", "ROI2elect_dist_abs_rel_cylinder_Cond_fit_vs_field_change", "Closest_elect_Rel_cylinder_Cond_fit_vs_field_change", "Rel_Cylinder_Cond_fit_vs_field_change", "Closest_elect_cylinder_Cond_fit_vs_field_change", "Cylinder_Cond_fit_vs_field_change", "AbsCondChange_vs_field_change", "Closest_elect_Rel_cylinder_Cond_vs_field_change", "Closest_elect_cylinder_Cond_vs_field_change", "Cylinder_Cond_vs_field_change", "ROI2elect_dist_Rel_cylinder_Cond_vs_field_change", "ROI2elect_dist_cylinder_Cond_vs_field_change", "Rel_Cylinder_Cond_vs_field_change"];
datetimeRun = strrep(strrep(strrep(datestr(datetime('now')), ' ', '_'), '-', '_'), ':', '_');
titleRun = append('Combined_', datetimeRun);
mkdir(titleRun);
% Define marker shape based on index
marker_shapes = {'o', '+', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
for j =1:length(plotNames)
    fig_handles = cell(length(headList), 1);
    % Initialize figure and axes
    combined_figs{j} = figure;
    for i=1:numel(headList)
      plotFile=join([convertStringsToChars(headList(i).name), '/', convertStringsToChars(plotNames(j)), convertStringsToChars(string(subjects(i))), '.fig']);
      strrep(plotFile,' ','');
      fig_handles{i} = openfig(plotFile);

      % Get handle to the current axes
      ax = gca;

      % Get handle to the current scatter plot
      h = findobj(ax, 'Type', 'scatter');

      % Get x and y data of the scatter plot
      x_data = get(h, 'XData');
      if length(x_data) == 0
          fail = 1;
      end
      y_data = get(h, 'YData');
      legend_data = ax.Children.CData;

      % Get x-axis and y-axis labels
      xlabel_str = ax.XLabel.String;
      ylabel_str = ax.YLabel.String;
      colorbarHandle = findobj(gcf, 'Type', 'colorbar');
      labelHandle = get(colorbarHandle, 'YLabel');
      colorbarLabel = get(labelHandle, 'String');
      
        
      close(fig_handles{i})

      % Plot the data from the current .fig file
      marker_shape = marker_shapes{mod(i, size(marker_shapes, 2))};

      if ~isa(x_data(1,1),'double')
          %% TODO: Understand why there are missing marker shapes. 
          x_data = cell2mat(x_data);
          y_data = cell2mat(y_data);
      end
      finalFig = scatter(x_data', y_data', [], legend_data, 'Marker', marker_shape);

      hold on

    end
    hold off

    % Add legend with the legend entries
    %legend(subjects);
    if ~isempty(colorbarLabel)
        c = colorbar;
        c.Label.String = colorbarLabel;
        c.Label.Rotation = 270;
    end

    % Add labels and title
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    try
        axis(gca, [min(x_data), max(x_data), min(y_data), max(y_data)])
        fail = 0;
    catch
        fail = 1;
    end
    title(xlabel_str + " vs " + ylabel_str);
    saveas(finalFig,append(titleRun, '/', convertStringsToChars(plotNames(j))))
    saveas(finalFig,append(titleRun, '/', convertStringsToChars(plotNames(j)), '.png'))
    close all


 
    

end