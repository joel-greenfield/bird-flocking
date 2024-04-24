function dvdt_ali = alignment(~,~,~,v_del,C_ali,~,~,~,n_n,M)
    dvdt_ali = zeros(size(v_del));
    for i = 1:M
        dvdt_ali = dvdt_ali+(v_del(n_n(:,i),:)-v_del);
    end
    dvdt_ali = C_ali/M * dvdt_ali;
end