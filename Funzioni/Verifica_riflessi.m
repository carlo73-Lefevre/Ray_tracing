% verifica sorgenti di superficii illuminate dai riflessi
% FOC_3_003
% Box_Wing
% Specchi


% Cannon_fine
% Mirror

clear;clc;

path = 'D:\OneDrive - INAF - Istituto Nazionale di Astrofisica\G4S\Matlab_ray';

if ~exist('DAT')
    DAT = matfile('Cannon_fine.mat');
    ray_b = DAT.ray_b;
    Sat_b = DAT.Sat_b;
    ACC   = DAT.ACC;
end

video = 'on';

if strcmp(video,'on')
    v = VideoWriter('specchi2.mp4','MPEG-4');
    open(v);
    f = figure('visible','off');
else
    figure(1);
end

for s1=1:length(ACC.beta) % beta

    for   s2 = 1:length(ACC.alf_x) %alfa

        ray_sun      = ray_b(s1,1,s2);
        ray_sour     = ray_b(s1,2,s2);
        ray_rif      = ray_b(s1,3,s2);
%%
        clf;
        p5 = patch('Faces',Sat_b(s1,1).Faces,'Vertices',Sat_b(s1,1).Vertex);%,'EdgeColor','none'
        p5.FaceColor = [0.7 0.7 0.7];
        hold on;
        p5.FaceAlpha = 0.9;
        axis equal;hold on;grid on;
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        %     xlim([-1.2 1.8])
%         ylim([-3 3])
        %     zlim([0.4 1.4])
        view(-130,30);
%         zoom(6)

%             p6 = patch('Faces',Sat_b(s1,s2).Faces(ray_sour.face_source,:),'Vertices',Sat_b(s1,s2).Vertex,'EdgeColor','none');
%             p6.FaceColor = 'r';
%             p6.FaceAlpha = 0.9; 

             [ACC_S1,ACC_N1] = SRP(Sat_b(s1,1),ray_sun,ray_sour,ray_sun.dir(1,:));
%              plot_quiv([-0.5 0 0], ray_sun.dir(1,:),0.2)

% %              plot_tri(ray_b(1,1).p_source+[7.5 0 0]);%,ray_b(1,1).dir,8,'b')
%              plot_quiv(ray_b(1,1).p_source+[7.5 0 0],ray_b(1,1).dir,12,'b')

%              plot_quiv(ray_sour.p_source,ray_sour.dir,3,'r')
view(0,90);
       %%
        % faccie di arrivo
        p7 = patch('Faces',Sat_b(s1,1).Faces(ray_b(1,2).face_source,:),'Vertices',Sat_b(s1,1).Vertex);
        p7.FaceColor = 'c';
        p7.FaceAlpha = 0.9;


        %% Seconda riflessione
        figure(2);
         clf;
        p5 = patch('Faces',Sat_b(s1,1).Faces,'Vertices',Sat_b(s1,1).Vertex);%,'EdgeColor','none'
        p5.FaceColor = [0.7 0.7 0.7];
        hold on;
%         p5.FaceAlpha = 0.9;
        axis equal;hold on;grid on;
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        view(-130,30);
%         zoom(6)

%             p6 = patch('Faces',Sat_b(s1,s2).Faces(ray_sour.face_source,:),'Vertices',Sat_b(s1,s2).Vertex,'EdgeColor','none');
%             p6.FaceColor = 'r';
%             p6.FaceAlpha = 0.9; 

%              plot_quiv([-0.5 0 0], ray_sun.dir(1,:),0.2)

% %              plot_tri(ray_b(1,1).p_source+[7.5 0 0]);%,ray_b(1,1).dir,8,'b')
%              plot_quiv(ray_b(1,1).p_source+[7.5 0 0],ray_b(1,1).dir,12,'b')

%              plot_quiv(ray_sour.p_source,ray_sour.dir,3,'r')

        p7 = patch('Faces',Sat_b(s1,1).Faces(ray_b(1,3).face_source,:),'Vertices',Sat_b(s1,1).Vertex);
        p7.FaceColor = 'm';
%         p7.FaceAlpha = 0.9;

        %%

%%        %             plot_quiv(ray_sour.p_source,ray_sour.dir,1)

        % plotta raggi riflessi 2
        for ii = 1:length(ray_rif.face_from)

            idx = (ray_sour.face_source == ray_rif.face_from(ii));


            line([ray_sour.p_source(idx,1) ray_rif.p_source(ii,1)],...
                [ray_sour.p_source(idx,2) ray_rif.p_source(ii,2)],...
                [ray_sour.p_source(idx,3) ray_rif.p_source(ii,3)],'Color','red')

        end

        if strcmp(video,'on')
            set(gca,'nextplot','replacechildren');
            frame = getframe(gcf);
            writeVideo(v,frame);
        end

    end
end
if strcmp(video,'on')
    close(v);
end