function h = Field_Hour_Values(r,N,cc,nl,Ccrown_c, ...
    Rn,H,QE,G,Ts,EG,ELitter,ESN,ESN_In,EWAT,EICE,EIn_urb,EIn_rock, ...
    Lk,Rh,Rd,T_H,T_L,EIn_H,EIn_L,An_H,An_L,Rdark_H,Rdark_L,O,Tdp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% One hour's values of every square, from row r of the 2-D frame's hourly
%% arrays (hours x k x squares), for Field_Hourly_Means. Ccrown_c is the
%% squares' crown cover that hour (squares x cc); nl the active-layer soil
%% layers kept. The x-y comparison calls it on a full record with r the
%% hour itself, so both build the identical values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h.Rn = one(Rn); h.H = one(H); h.QE = one(QE); h.G = one(G); h.Ts = one(Ts);
h.EG = one(EG); h.ELitter = one(ELitter); h.ESN = one(ESN);
h.ESN_In = one(ESN_In); h.EWAT = one(EWAT); h.EICE = one(EICE);
h.EIn_urb = one(EIn_urb); h.EIn_rock = one(EIn_rock);
h.Lk = one(Lk); h.Rh = one(Rh); h.Rd = one(Rd);
h.T_H = slots(T_H); h.T_L = slots(T_L); h.EIn_H = slots(EIn_H); h.EIn_L = slots(EIn_L);
h.An_H = slots(An_H); h.An_L = slots(An_L);
h.Rdark_H = slots(Rdark_H); h.Rdark_L = slots(Rdark_L);
h.Ccrown = Ccrown_c.';
h.O = reshape(O(r,1:nl,:),nl,N);
h.Tdp = reshape(Tdp(r,1:nl,:),nl,N);

    function x = one(X)
        x = reshape(X(r,1,:),1,N);
    end
    function x = slots(X)
        x = reshape(X(r,:,:),cc,N);
    end
end
