function Dentro = ver_inner(A,B,C,P,varargin)
% A B C  vertici triangolo
% P punto da valutare se interno

ba = B-A;
bc = B-C;
ca = C-A;

Area_sect = norm(cross(ba,ca))/2 ;  %calc_area(ba,bc,ca);
Area_1    = norm(cross(ba,A-P))/2;  %calc_area(ba,A-P,B-P);
Area_2    = norm(cross(bc,B-P))/2;  %calc_area(bc,B-P,C-P);
Area_3    = norm(cross(ca,C-P))/2;  %calc_area(ca,C-P,A-P);

if ~isempty(varargin)
    
    if ((Area_1+Area_2+Area_3)- Area_sect) <= varargin{1,1}%|| ...
        %floor(Area_sect-(Area_1+Area_2+Area_3)) == 0
        Dentro = true;
    else
        Dentro = false;
    end

else

    if ((Area_1+Area_2+Area_3)- Area_sect) <= 1e-6%|| ...
        %floor(Area_sect-(Area_1+Area_2+Area_3)) == 0
        Dentro = true;
    else
        Dentro = false;
    end
end
end