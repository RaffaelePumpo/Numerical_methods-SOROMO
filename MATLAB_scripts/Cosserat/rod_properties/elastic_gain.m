function   dy = elastic_gain(X,~,Const)


L= Const.L;


B = Const.B;

H = diag([Const.GI,Const.EI,Const.EI,Const.EA,Const.GA,Const.GA]);
D = Const.mu*H;



Ha  = B'*H*B;
Da  = B'*D*B;

Phi = getPhi(X,Const);

Krr_prime = Phi'*Ha*Phi;
Drr_prime = Phi'*Da*Phi;

Krr_prime_vec = reshape(Krr_prime,[Const.dim_base*Const.dim_base, 1]);
Drr_prime_vec = reshape(Drr_prime,[Const.dim_base*Const.dim_base, 1]);


% sorties
dy  = L*[Krr_prime_vec;
    Drr_prime_vec];

end