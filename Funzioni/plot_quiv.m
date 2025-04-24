function q = plot_quiv(C,V,varargin)
   % C punto sorgente
   % V direzione 
c = 1;
% for ii =1:length(V(:,2))
%     if abs(V(ii,2)) > 0.1
%         VV(c,:)  = V(ii,:);
%         CC(c,:) = C(ii,:);
%         c = c+1;
%     end
% end

% V =VV;
% C =CC;


   if size(varargin)>0

        q = quiver3(C(:,1),C(:,2),C(:,3),V(:,1),V(:,2),V(:,3),varargin{:});
           q.ShowArrowHead = 'on';
           q.Marker = '.';
           

   else
        q = quiver3(C(:,1),C(:,2),C(:,3),V(:,1),V(:,2),V(:,3));
           q.ShowArrowHead = 'off';
           q.Marker = '.';

   end
   if size(varargin)>2
        q = quiver3(C(:,1),C(:,2),C(:,3),V(:,1),V(:,2),V(:,3),varargin{1},string(varargin{2}));
        q.ShowArrowHead = 'off';
        q.Marker = '.';
   end

end
%patch di 3 punti

