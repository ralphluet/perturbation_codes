function denHaanLoop(x,distrSS, XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC,DCD,IDCD)

    JDminus = copy(distrSS)

    m_SS          = dropdims(sum(distrSS,dims=(2,3)),dims=(2,3))

    k_SS          = dropdims(sum(distrSS,dims=(1,3)),dims=(1,3))

    h_SS          = dropdims(sum(distrSS,dims=(1,2)),dims=(1,2))

    Brealized_RC=ones(nlag).*exp(XSS[indexes.BSS]);
    Krealized_RC=ones(nlag).*exp(XSS[indexes.KSS]);

    RBminus=copy(exp(XSS[indexes.RBSS]))

    xHAAN      = zeros(size(LOMstate,1), 1)

    PG1=pinv(Γ[1])
    PG2=pinv(Γ[2])

    for t=1:nlag-1#mpar.T(1,1)-1
        xHAAN[indexes.Z] = x[indexes.Z,t]
        IRFX = (MX*xHAAN)
        xHAANPrime = LOMstate * xHAAN
        IRFXPrime = MX*xHAANPrime

        QuantDiff, RBminus, JDminus  = denHaan(JDminus, RBminus, IRFX,IRFXPrime ,XSS,m_par,n_par,indexes,Γ,compressionIndexes,DC,IDC,DCD,IDCD)

        xHAAN = LOMstate * xHAAN

        Brealized_RC[t+1]=dot(JDminus,n_par.mesh_m)
        Krealized_RC[t+1]=dot(JDminus,n_par.mesh_k)

        temp          = dropdims(sum(JDminus,dims=(2,3)),dims=(2,3))
        # xHAAN[indexes.distr_m] = PG1*(temp[1:end] - m_SS[1:end])
        xHAAN[indexes.distr_m] = (temp[1:end-1] - m_SS[1:end-1])
        temp          = dropdims(sum(JDminus,dims=(1,3)),dims=(1,3))
        # xHAAN[indexes.distr_k] = PG2*(temp[1:end] - k_SS[1:end])
        xHAAN[indexes.distr_k] = (temp[1:end-1] - k_SS[1:end-1])

        # temp          = dropdims(sum(JDminus,dims=(1,2)),dims=(1,2))
        # xHAAN[indexes.distr_y] = (cum_h[1:end-1] - cum_h_SS[1:end-1])

        xHAAN[indexes.RB] = (log(RBminus) - XSS[indexes.RBSS])

        println(RBminus)


    end

    return Krealized_RC, Brealized_RC, RBminus, xHAAN

end
