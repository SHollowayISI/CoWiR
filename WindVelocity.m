%% Housekeeping

% clear variables
close all
addpath(genpath('Velocity Figures'));
addpath(genpath('MAT Files'));
addpath(genpath('Figures'));


for loop = 1:length(files)
    
    %% Unpack Data
    
    % Unpack data from Figures
    %{
    filename = files{loop};
    
    for chan = 1:3
        
        filepath = sprintf([filename '_ch%d_velocity.fig'], chan);
        
        velgraph = openfig(filepath);
        
        axobjs = velgraph.Children;
        dataobjs = axobjs(2).Children;
        
        range(chan,:) = dataobjs(1).XData;
        vels(chan,:) = dataobjs(1).YData;
        
%         close all
        
    end
    %}
    
    % Unpack data from file
    
    for chan = 1:3
        filename = ['MAT Files\' seq_name, '_ch', sprintf('%d', chan), '.mat'];
        
        vars_in = load(filename);
        vels(chan,:) = vars_in.mean_vel;
        range(chan,:) = vars_in.range_axis(1,1:size(vels,2));
    end
        
    
    
    % figure;
    % plot(range(1,:), vels(1,:), range(2,:), vels(2,:), range(3,:), vels(3,:))
    
    
    %% Calculate Wind Direction
    
    % a = 1.38;
    a = 3.2;
    b = 30;
    
    % wind_mat = ...
    %     [[-sind(a)*sind(60), sind(a)*cosd(60), cosd(a)]; ...
    %      [sind(a), 0, cosd(a)]; ...
    %      [-sind(a)*sind(60), -sind(a)*cosd(60), cosd(a)]];
    
    % vels = [vels(3,:); vels(2,:); vels(1,:)];
    
    % wind_vel = ((wind_mat)\vels);
    
    wind_vel(1,:) = (2/(3*sind(a)))*(vels(2,:) - 0.5*vels(1,:) - 0.5*vels(3,:));
    wind_vel(2,:) = (1/(sqrt(3)*sind(a)))*(vels(1,:)+vels(3,:));
    wind_vel(3,:) = (1/(3*cosd(a)))*(vels(1,:)+vels(2,:)+vels(3,:));
    
    % plot(wind_vel)
    
    %% Visualization
    
%     close all
    
    ymax = 40;
    
    % X, Y, Z-Direction Subplots
    %
    figure('Name', 'Directional_Velocity_Subplot')
    
    plot_title = ['X-Direction'; 'Y-Direction'; 'Z-Direction'];
    
    for ind = 1:3
        subplot(3,1,ind)
        plot(range(ind,:), wind_vel(ind,:), 'LineWidth', 2);
        title(plot_title(ind,:))
        
        grid on;
        set(gca, 'YLim', [-ymax, ymax], 'FontWeight', 'bold');
        xlabel('Range [m]');
        ylabel('Estimated Velocity [m/s]');
        
    end
    %}
    
    % X, Y, Z-Direction Subplots
    %
    figure('Name', 'Directional_Velocity')
    
    plot_title = ['X-Direction'; 'Y-Direction'; 'Z-Direction'];
    plot_color = ['r', 'g', 'b'];
    
    for ind = 1:3
        %     subplot(3,1,ind)
        plot(range(ind,:), wind_vel(ind,:), 'LineWidth', 2, ...
            'DisplayName', plot_title(ind,:), ...
            'Color', plot_color(ind));
        hold on
        %     title(plot_title(ind,:))
        
        grid on;
        set(gca, 'YLim', [-ymax, ymax], 'FontWeight', 'bold');
        xlabel('Range [m]');
        ylabel('Estimated Velocity [m/s]');
        
    end
    legend
    %}
    
    % Magnitude Plot
    %
    figure('Name', 'Magnitude')
    
    mag_vel = squeeze(sqrt(sum(wind_vel.^2, 1)));
    
    plot(range(1,:), mag_vel, 'LineWidth', 2);
    title('Wind Velocity Magnitude')
    grid on;
    set(gca, 'YLim', [0, ymax], 'FontWeight', 'bold');
    xlabel('Range [m]');
    ylabel('Estimated Velocity [m/s]');
    %}
    
    
    % Tangent Plot
    %
    figure('Name', 'Tangent')
    
    tangent_vel = atand(vels(2,4:end)./vels(3,4:end));
    
    plot(range(1,4:end), tangent_vel, 'LineWidth', 2);
    title('Wind Direction')
    grid on;
    set(gca, 'YLim', [-90, 90], 'FontWeight', 'bold');
    xlabel('Range [m]');
    ylabel('Direction [degree]');
    %}
    
    %% Save plots to file
    %
    filepath = ['Figures\Wind_Velocity\', ...
        seq_name, '\'];
    if ~exist(filepath, 'dir')
        mkdir(filepath)
    end
    
    FolderName = filepath;   % Your destination folder
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = get(FigHandle, 'Name');
        savefig(FigHandle, fullfile(FolderName, [seq_name, '_', FigName, '.fig']));
        saveas(FigHandle, fullfile(FolderName, [seq_name, '_', FigName, '.png']));
    end
    
    close all
    %}
    
end
