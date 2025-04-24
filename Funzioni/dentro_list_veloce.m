function IDX_Facce_ill = dentro_list_veloce(A,B,C,P0,IDX_Facce_ill,cont,varargin)

    if cont >= length(IDX_Facce_ill)
        return
    end

    P = repmat(P0(cont,:),size(A,1), 1);

    ba = B-A;
    bc = B-C;
    ca = C-A;

    Area_sect = 0.5*norm(cross(ba,ca),2);
    Area_1    = 0.5*norm(cross(ba,A-P),2);
    Area_2    = 0.5*norm(cross(bc,B-P),2);
    Area_3    = 0.5*norm(cross(ca,C-P),2);

    if ~isempty(varargin)
        Dentro_idx = all((Area_1+Area_2+Area_3)- Area_sect <= varargin{1,1});
    else
        Dentro_idx = all((Area_1+Area_2+Area_3)- Area_sect <= 1e-6);
    end

    l1 = find(Dentro_idx);
    IDX_Facce_ill(l1) = 0;
    cont = cont +1;

    IDX_Facce_ill = dentro_list_veloce(A,B,C,P0,IDX_Facce_ill,cont,1e-11);

end
