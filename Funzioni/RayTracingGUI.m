function RayTracingGUI(Sat, rays)
    % GUI per visualizzare dati 3D e interagire con superfici e raggi riflessi

    % Creazione della finestra principale
    fig = figure('Name', 'Ray Tracing 3D Viewer', 'Position', [100, 100, 1000, 600]);

    % Creazione di una struttura per salvare lo stato della GUI
    guiData.SelectedSurface = [];
    guiData.SelectedRays = [];
    guidata(fig, guiData);

    % Pannello per la visualizzazione 3D
    ax = axes('Parent', fig, 'Position', [0.05, 0.1, 0.6, 0.8]);
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    view(ax, 3);

    % Visualizzazione iniziale del satellite
    p = patch('Faces', Sat.Faces, 'Vertices', Sat.Vertex, 'Parent', ax);
    p.FaceColor = [0.7, 0.7, 0.7];
    p.EdgeColor = 'k';
    p.FaceAlpha = 0.8;
    title(ax, '3D Satellite Viewer');

    % Pannello per i controlli
    ctrlPanel = uipanel('Title', 'Controls', 'FontSize', 12, ...
                        'Position', [0.7, 0.1, 0.28, 0.8]);

    % Pulsante per selezionare una superficie
    uicontrol(ctrlPanel, 'Style', 'pushbutton', 'String', 'Select Surface', ...
              'Position', [10, 350, 150, 30], ...
              'Callback', @(~, ~) selectSurface(ax, Sat, fig));

    % Checkbox per attivare/disattivare i raggi riflessi
    chkReflectedRays = uicontrol(ctrlPanel, 'Style', 'checkbox', 'String', 'Show Reflected Rays', ...
                                 'Position', [10, 300, 150, 30], ...
                                 'Callback', @(src, ~) toggleRays(ax, rays, src.Value, fig));

    % Slider per trasparenza delle superfici
    uicontrol(ctrlPanel, 'Style', 'text', 'String', 'Surface Transparency:', ...
              'Position', [10, 250, 150, 20], 'HorizontalAlignment', 'left');
    sliderTransparency = uicontrol(ctrlPanel, 'Style', 'slider', ...
                                   'Position', [10, 220, 150, 20], ...
                                   'Min', 0, 'Max', 1, 'Value', 0.8, ...
                                   'Callback', @(src, ~) adjustTransparency(p, src.Value));

    % Pulsante per reimpostare la vista
    uicontrol(ctrlPanel, 'Style', 'pushbutton', 'String', 'Reset View', ...
              'Position', [10, 180, 150, 30], ...
              'Callback', @(~, ~) resetView(ax, Sat));

    % Aggiunta di un'area per informazioni
    lblSurfaceInfo = uicontrol(ctrlPanel, 'Style', 'text', 'String', 'Selected Surface Info: None', ...
                               'Position', [10, 100, 150, 40], ...
                               'HorizontalAlignment', 'left');

    % Funzioni utili

    function selectSurface(ax, Sat, fig)
        % Funzione per selezionare una superficie cliccando
        disp('Select a surface by clicking on the mesh.');
        dcm = datacursormode(fig);
        datacursormode on;
        pause;
        info = getCursorInfo(dcm);
        if ~isempty(info)
            selectedFaceIdx = info.DataIndex;
            lblSurfaceInfo.String = sprintf('Selected Face: %d', selectedFaceIdx);
            highlightSurface(ax, Sat, selectedFaceIdx);

            % Salva lo stato nella GUI
            guiData = guidata(fig);
            guiData.SelectedSurface = selectedFaceIdx;
            guidata(fig, guiData);
        end
        datacursormode off;
    end

    function highlightSurface(ax, Sat, idx)
        % Evidenzia la superficie selezionata
        patch('Faces', Sat.Faces(idx, :), 'Vertices', Sat.Vertex, ...
              'FaceColor', 'r', 'EdgeColor', 'none', 'Parent', ax);
    end

    function toggleRays(ax, rays, show, fig)
        % Mostra o nasconde i raggi riflessi
        guiData = guidata(fig);
        if show
            guiData.SelectedRays = quiver3(ax, ...
                                           rays.ray_source_P(:, 1), rays.ray_source_P(:, 2), rays.ray_source_P(:, 3), ...
                                           rays.ray_rif_dir(:, 1), rays.ray_rif_dir(:, 2), rays.ray_rif_dir(:, 3), ...
                                           0.5, 'r');
        else
            if ~isempty(guiData.SelectedRays)
                delete(guiData.SelectedRays);
                guiData.SelectedRays = [];
            end
        end
        guidata(fig, guiData);
    end

    function adjustTransparency(p, alpha)
        % Modifica la trasparenza della mesh
        p.FaceAlpha = alpha;
    end

    function resetView(ax, Sat)
        % Reimposta la vista iniziale
        cla(ax);
        patch('Faces', Sat.Faces, 'Vertices', Sat.Vertex, ...
              'Parent', ax, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'k');
        view(ax, 3);
        grid(ax, 'on');
    end
end
