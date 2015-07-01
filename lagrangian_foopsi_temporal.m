function [c,ld] = lagrangian_foopsi_temporal(y,a,thr,G,ld_in)

% solves min sum(G*c) 
%       s.t.: G*c >= 0, ||a(j)*c - y(j)|| <= thr(j)
T = size(G,1);
myfun = @(Ald) lagrangian_temporal_grad(Ald,y,thr);
options = optimset('GradObj','On','Display','Off','Algorithm','interior-point','TolX',1e-3);
if nargin < 5
    ld_in = 10*ones(length(a),1);
end
ld = fmincon(myfun,ld_in,[],[],[],[],zeros(length(a),1),[],[],options);

    function [f,grad] = lagrangian_temporal_grad(Al,y,thr)
        v = G'*ones(T,1);
        options2 = optimset('Display','Off','Algorithm','interior-point-convex');
        c = quadprog(2*norm(a)^2*speye(T),-2*y'*a+v,-G,zeros(T,1),[],[],[],[],[],options2);

        f = v'*c; 
        for i = 1:length(a)
            f = f + Al(i)*(norm(a(i)*c'-y(i,:))^2 - thr(i));
        end
        grad = zeros(length(a),1);
        for i = 1:length(a)
            grad = grad + norm(a(i)*c'-y(i,:))^2 - thr(i);
        end
        f = -f;
        grad = -grad;
    end

end