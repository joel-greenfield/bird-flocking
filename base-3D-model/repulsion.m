function dvdt_rep = repulsion(~,~,x_del,~,~,~,C_rep,ep,n_n,M)
    dvdt_rep = zeros(size(x_del));
    for i = 1:M
        num = x_del(n_n(:,i),:)-x_del;
        den = sqrt(sum((x_del(n_n(:,i),:)-x_del).^2,2))+ep;
        dvdt_rep = dvdt_rep+num./den;
    end
    dvdt_rep = -C_rep * dvdt_rep;
end