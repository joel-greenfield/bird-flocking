%% Attraction force
function dvdt_att = attraction(~,~,x_del,~,~,C_att,~,~,n_n,M)
    dvdt_att = zeros(size(x_del));
    for i = 1:M
        dvdt_att = dvdt_att+(x_del(n_n(:,i),:)-x_del);
    end
    dvdt_att = C_att * dvdt_att;
end