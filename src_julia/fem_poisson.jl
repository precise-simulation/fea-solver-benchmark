function fem_poisson( n = 16 )
# Solves uxx + uyy = 1 on a unit square with Q1 finite elements.
# Returns the computed solution vector in u, and the timimngs vector where:
#
# t_grid   = timings[1] time for grid generation
# t_ptr    = timings[2] time to calculate matrix pointers
# t_asm    = timings[3] time for system matrix assembly
# t_rhs    = timings[4] time for right hand side assembly
# t_bdr    = timings[5] time to set boundary conditions
# t_sparse = timings[6] time to convert to sparse matrix format
# t_sol    = timings[7] time for solution

# Copyright 2013-2018 Precise Simulation, Ltd.


    # Grid
    gc_enable(false)
    tic()
    nx = n
    ny = n
    grid = rectgrid( nx, ny );
    neq = size( grid[1], 2 )
    t_grid = toq()
    gc_enable(true)
    gc()


    # CSR matrix pointers
    gc_enable(false)
    tic()
    kcol, kld, na = ap7( grid[2] )
    t_ptr = toq()
    gc_enable(true)
    gc()


    # Assemble system matrix
    gc_enable(false)
    tic()
    avals = assemblea( grid[1], grid[2], kcol, kld, na, neq )
    t_asm = toq()
    gc_enable(true)
    gc()


    # Assemble right hand side
    gc_enable(false)
    tic()
    f = assemblef( grid[1], grid[2], kcol, kld, na, neq )
    t_rhs = toq()
    gc_enable(true)
    gc()


    # Boundary conditions
    tic()
    set_boundary_conditions( nx, ny, avals, kld, f )
    t_bdr = toq()
    gc()


    # Sparsification
    tic()
    rowptr, colptr = convert_csr_to_triplet( avals, kcol, kld )
    A = sparse( rowptr, colptr, avals, neq, neq )
    t_sparse = toq()
    gc()


    # Sparse mv
    tic()
    for i = 1:100
        tmp = A*f
    end
    t_sparse_mv = toq()/100.0
    gc()


    # Solve system
    tic()
    u = A\f
    t_sol = toq()


    return u, [ t_grid, t_ptr, t_asm, t_rhs, t_bdr, t_sparse, t_sparse_mv, t_sol ]

end


function set_boundary_conditions( nx, ny, avals, kld, f )

    nx = nx + 1
    ny = ny + 1
    bind = [ collect(1:nx)'                 collect(nx:nx:nx*ny)'
             collect(nx*ny:-1:nx*ny-nx+1)'  collect(nx*ny-nx+1:-nx:1)' ]
    bind = sort( unique( bind ) )

    # Set matrix rows to zero with diagonal entry to one
    for ib = 1:length(bind)
        ieq = bind[ib]
        for ia = kld[ieq]:kld[ieq+1]-1
            avals[ia] = 0.0
        end
        avals[kld[ieq]] = 1.0
    end

    # Set rhs to zero
    f[bind] = 0.0

end


function rectgrid( n_cx = 10, n_cy = 0, xp = [0 1;0 1] )


    x = []
    y = []
    if n_cx <= 0
       n_cx = 10
    end
    if n_cy <= 0
       n_cy = n_cx
    end
    if length(n_cx)>1
        xp[1] = min(n_cx)
        xp[3] = max(n_cx)
        x = sort(n_cx)
        n_cx = length(n_cx) - 1
    end
    if length(n_cy)>1
        xp[2] = min(n_cy)
        xp[4] = max(n_cy)
        y = sort(n_cy)
        n_cy = length(n_cy) - 1
    end

    lx   = xp[3] - xp[1]   # Length in x-direction.
    ly   = xp[4] - xp[2]   # Length in y-direction.
    x0   = xp[1]           # x-coordinate of lower left corner.
    y0   = xp[2]           # y-coordinate of lower left corner.
    n_c  = n_cx*n_cy
    n_px = n_cx + 1
    n_py = n_cy + 1
    n_p  = n_px*n_py


    # Grid point vertex coordinates (grid points ordered
    # in rows from left to right and top to bottom).
    if isempty(x)
        x = linspace( x0, x0+lx, n_px )
    end
    if isempty(y)
        y = linspace( y0, y0+ly, n_py )
    end
    yy = Array{Float64}(n_p)
    for j = 1:n_py
        yy[(j-1)*n_px+1:j*n_px] = y[j]
    end
    xx = repmat( x, n_py, 1 )
    p = [ xx yy ]'

    # Cell cell connectivities (grid points ordered counterclockwise for
    # each cell starting with the bottom left vertex as the first node).
    c = Array{Int64}(4,n_c)
    ic = 0
    for j = 1:n_py-1
        for i = 1:n_px-1
            ic += 1
            c[1,ic] = i + (j-1)*n_px
            c[2,ic] = i + 1 + (j-1)*n_px
            c[3,ic] = i + 1 + j*n_px
            c[4,ic] = i + j*n_px
        end
    end

    # Grid cell adjacencies.
    a = Array{Int64}(4,n_c)
    ic = 0
    for j = 1:n_cy
        for i = 1:n_cx
            ic += 1
            a[1,ic] = ic - n_cx
            a[2,ic] = ic + 1
            a[3,ic] = ic + n_cx
            a[4,ic] = ic - 1
        end
        a[2,ic] = 0
        a[4,ic-n_cx+1] = 0
    end
    a[1,1:n_cx] = 0
    a[3,n_c-n_cx+1:n_c] = 0


    # Boundary information.
    ox = ones( 1, n_cx )
    oy = ones( 1, n_cy )
    b = [ (1:n_cx)'      (n_cx:n_cx:n_cx*n_cy)'  (n_c:-1:n_c-n_cx+1)'  (n_c-n_cx+1:-n_cx:1)' ;
          ox             2*oy                    3*ox                  4*oy                  ;
          ox             2*oy                    3*ox                  4*oy                  ;
          zeros(1,n_cx)  oy                      zeros(1,n_cx)        -oy                    ;
         -ox             zeros(1,n_cy)           ox                    zeros(1,n_cy)         ]

    # Subdomain indexes.
    s = ones(Int64,1,n_c)

    # Output.
    grid = Array{Any}(5)

    grid[1] = p
    grid[2] = c
    grid[3] = a
    grid[4] = b
    grid[5] = s

    return grid

end


function convert_csr_to_triplet( avals, kcol, kld )

    na  = kld[end] - 1
    neq = length( kld ) - 1

    rowptr = zeros( Int64, na )
    colptr = zeros( Int64, na )
    cnt = 0
    for irow = 1:neq
        for j = kld[irow]:kld[irow+1]-1
            jcol = kcol[j]
            cnt  = cnt + 1
            rowptr[cnt] = irow
            colptr[cnt] = jcol
        end
    end

    return rowptr, colptr

end


function ap7( kvert::Array{Int64,2} )

    nel = size(kvert,2)
    neq = maximum(kvert)
    namax = 15*neq

    kcol = Array{Int64}(namax)
    kld  = zeros(Int64,neq+1)

    na    = neq
    kcol1 = Array{Int64}(namax)
    kcol1[1:neq] = 1:neq
    kind  = Array{Int64}(namax)
    kind[1:neq]  = 0

    kdfg = zeros(Int64,4)
    kdfl = zeros(Int64,4)

    bsymm  = false
    idofe1 = 1
    idfl   = 4   #ndfl[ieltyp]   # Determine number of degrees of freedom per element

    jvg = zeros(Int64,8)
    jvl = zeros(Int64,8)

    for iel = 1:nel   # 100

        # kdfg returns the global degrees of freedom in increasing order
        kdfg, kdfl = ndfgl( iel, 0, kvert, kdfg, kdfl, jvg, jvl )

        # loop over local number of degrees of freedom
        for idofe = 1:idfl   # 110
            irow = kdfg[idofe]

            if bsymm
                if irow == neq
                    continue
                end
                idofe1 = idofe
            end

            # Loop over off-diagonal elements
            # only upper triangular part in symmetric case
            for jdofe = idofe1:idfl   # 120
                if idofe == jdofe
                    continue
                end
                jcol = kdfg[jdofe]
                # JCOL is the global number of d.o.f. to be inserted into row IROW

                # Look whether entry (IROW,JCOL) is already provided
                ipos = irow
                # @label lbl_121
                while true
                    if kind[ipos] == 0   # 121

                        # insert new entry
                        na = na + 1
                        kcol1[na] = jcol
                        kind[ipos] = na
                        kind[na] = 0
                        break

                    else

                        ipos = kind[ipos]
                        # if kcol1[ipos] != jcol
                        #     @goto lbl_121
                        # end

                    end
                end

            end   # 120

        end   # 110

    end   # 100

    # Collect entries on KCOL1 separately for each row
    na = 0
    for ieq = 1:neq   # 200
        na = na + 1
        kld[ieq] = na
        kcol[na] = ieq
        ipos = ieq

        while kind[ipos] != 0   # 201
            na = na + 1
            ipos = kind[ipos]
            kcol[na] = kcol1[ipos]
        end

    end   # 200
    kld[neq+1] = na + 1

    # Sort off-diagonal entries on KCOL separately for each row
    for ieq = 1:neq   # 300

        # @label lbl_301
        bsort = false
        while !bsort
            bsort = true
            for icol = kld[ieq]+1:kld[ieq+1]-2
                if kcol[icol] > kcol[icol+1]
                    ihelp = kcol[icol]
                    kcol[icol] = kcol[icol+1]
                    kcol[icol+1] = ihelp
                    bsort = false
                end
            end
        end
        # if !bsort
            # @goto lbl_301
        # end

    end   # 300
    deleteat!(kcol,na+1:length(kcol))

    return kcol, kld, na

end


function assemblea( dcorvg::Array{Float64,2}, kvert::Array{Int64,2},
                    kcola::Array{Int64,1},    klda::Array{Int64,1}, na::Int64, neq::Int64 )

    # initializations
    a = zeros(Float64,na)

    aux  = 1.0  # coefficient for bilinear form
    bder = [false true true]   # basis function evaluations

    # related to grid cell
    nel  = size(kvert,2)
    nve  = 4
    kve  = zeros(Int64,nve)
    dx   = zeros(Float64,nve)
    dy   = zeros(Float64,nve)

    # for elem e011
    djac   = zeros(Float64,2,2)
    dbas   = zeros(Float64,4,3)
    dhelp  = zeros(Float64,4,2)
    idfl   = 4
    kentry = zeros(Int64,idfl,idfl)
    dentry = zeros(Float64,idfl,idfl)
    kdfg   = zeros(Int64,idfl)
    kdfl   = collect(1:idfl)

    # 2x2 Gauss rule   # call cb2q(icub)
    ncubp  = 4
    c      = 0.577350269189626
    dxi    = [ -c -c; c -c; -c c; c c ]
    domega = ones(Float64,4);

    jvg = zeros(Int64,8)
    jvl = zeros(Int64,8)

    # loop over all cells
    for iel = 1:nel

        # get local to global dof mapping
        kdfg, foo1 = ndfgl( iel, 1, kvert, kdfg, kdfl, jvg, jvl )

        # determine entry positions in matrix
        iaj = 0
        for irowl = 1:idfl
            irow = kdfg[irowl]
            ia   = klda[irow]           # pointer to diagonal entry in a(irow,irow), and start of row irow
            kentry[irowl,irowl] = ia
            dentry[irowl,irowl] = 0.0
            for jcoll = 1:idfl
                if jcoll != irowl             # diagonal already processed
                    jcol = kdfg[jcoll]
                    for iaj = (ia+1):na       # loop over row entries
                        if kcola[iaj] == jcol
                            break
                        end
                    end
                    ia = ia + 1   # possibly not necessary if kdfl is sorted?
                    kentry[irowl,jcoll] = iaj
                    dentry[irowl,jcoll] = 0.0
                end
            end
        end

        # extract vertex coordinates
        for ive = 1:nve
          jp       = kvert[ive,iel]
          kve[ive] = jp
          dx[ive]  = dcorvg[1,jp]
          dy[ive]  = dcorvg[2,jp]
        end

        dj1 = 0.5*( -dx[1] - dx[2] + dx[3] + dx[4] )
        dj2 = 0.5*(  dx[1] - dx[2] + dx[3] - dx[4] )
        dj3 = 0.5*( -dy[1] + dy[2] - dy[3] + dy[4] )
        dj4 = 0.5*( -dy[1] + dy[2] + dy[3] - dy[4] )

        # loop over all cubature points
        for icubp = 1:ncubp

            xi1 = dxi[icubp,1]
            xi2 = dxi[icubp,2]

            # jacobian of the bilinear mapping onto the reference element
            djac[1,1] = 0.5*(dx[2]-dx[1]+dj2) + 0.5*dj2*xi2
            djac[1,2] = 0.5*dj1 + 0.5*dj2*xi1
            djac[2,1] = 0.5*dj4 - 0.5*dj3*xi2
            djac[2,2] = 0.5*(dy[3]-dy[1]-dj4) - 0.5*dj3*xi1
            detj      = djac[1,1]*djac[2,2] - djac[1,2]*djac[2,1]
            om        = domega[icubp]*detj

            # evaluate basis function
            dbas = e011( xi1, xi2, djac, bder, dbas, dhelp, detj )

            # summing up over all pairs of multiindices
            for jdofe = 1:idfl
                jdofeh = kdfl[jdofe]
                hbasj2 = dbas[jdofeh,2]
                hbasj3 = dbas[jdofeh,3]

                for idofe = 1:idfl
                    if idofe == jdofe
                        ah     = aux*(hbasj2^2 + hbasj3^2)
                        # ah     = aux*dbas[jdofeh,1]^2   # mass matrix
                    else
                        idofeh = kdfl[idofe]
                        hbasi2 = dbas[idofeh,2]
                        hbasi3 = dbas[idofeh,3]
                        ah     = aux*(hbasj2*hbasi2 + hbasj3*hbasi3)
                        # ah     = aux*dbas[jdofeh,1]*dbas[idofeh,1]   # mass matrix
                    end

                    dentry[idofe,jdofe] = dentry[idofe,jdofe] + om*ah
                end
            end

        end   # end loop over cubature points


        for jdofe = 1:idfl
            for idofe = 1:idfl
                ia    = kentry[jdofe,idofe]
                a[ia] = a[ia] + dentry[jdofe,idofe]
            end
        end

    end   # end loop over cells

    return a

end


function assemblef( dcorvg::Array{Float64,2}, kvert::Array{Int64,2}, kcola::Array{Int64,1}, klda::Array{Int64,1}, na::Int64, neq::Int64 )

    # initializations
    f = zeros(Float64,neq,1)

    aux  = 1.0  # coefficient for bilinear form
    bder = [true false false]   # basis function evaluations

    # related to grid cell
    nel  = size(kvert,2)
    nve  = 4
    kve  = zeros(Int64,nve)
    dx   = zeros(Float64,nve)
    dy   = zeros(Float64,nve)

    # for elem e011
    djac   = zeros(Float64,2,2)
    dbas   = zeros(Float64,4,3)
    dhelp  = zeros(Float64,4,2)
    idfl   = 4
    kdfg   = zeros(Int64,idfl)
    kdfl   = collect(1:idfl)

    # 2x2 Gauss rule   # call cb2q(icub)
    ncubp  = 4
    c      = 0.577350269189626
    dxi    = [ c c; -c c; -c -c; c -c ]
    domega = ones(Float64,4);

    jvg = zeros(Int64,8)
    jvl = zeros(Int64,8)

    # loop over all cells
    for iel = 1:nel

        # get local to global dof mapping
        kdfg, foo1 = ndfgl( iel, 1, kvert, kdfg, kdfl, jvg, jvl )

        # extract vertex coordinates
        for ive=1:nve
          jp       = kvert[ive,iel]
          kve[ive] = jp
          dx[ive]  = dcorvg[1,jp]
          dy[ive]  = dcorvg[2,jp]
        end

        dj1 = 0.5*( -dx[1] - dx[2] + dx[3] + dx[4] )
        dj2 = 0.5*(  dx[1] - dx[2] + dx[3] - dx[4] )
        dj3 = 0.5*( -dy[1] + dy[2] - dy[3] + dy[4] )
        dj4 = 0.5*( -dy[1] + dy[2] + dy[3] - dy[4] )

        # loop over all cubature points
        for icubp = 1:ncubp

            xi1 = dxi[icubp,1]
            xi2 = dxi[icubp,2]

            # jacobian of the bilinear mapping onto the reference element
            djac[1,1] = 0.5*(dx[2]-dx[1]+dj2) + 0.5*dj2*xi2
            djac[1,2] = 0.5*dj1 + 0.5*dj2*xi1
            djac[2,1] = 0.5*dj4 - 0.5*dj3*xi2
            djac[2,2] = 0.5*(dy[3]-dy[1]-dj4) - 0.5*dj3*xi1
            detj      = djac[1,1]*djac[2,2] - djac[1,2]*djac[2,1]
            om        = domega[icubp]*detj

            # evaluate basis function
            dbas = e011( xi1, xi2, djac, bder, dbas, dhelp, detj )

            # summing up over all pairs of multiindices
            for jdofe = 1:idfl
                ieq = kdfg[jdofe]
                f[ieq] = f[ieq] + aux*dbas[kdfl[jdofe],1]*om
            end

        end   # end loop over cubature points

    end   # end loop over cells

    return f

end


function ndfgl( iel::Int64, ipar::Int64, kvert::Array{Int64,2}, kdfg::Array{Int64}, kdfl::Array{Int64}, jvg::Array{Int64}, jvl::Array{Int64} )

    nve = size(kvert,1)
    # jvg = zeros(Int64,8)
    # jvl = zeros(Int64,8)

    for ive = 1:nve
        jvg[ive] = kvert[ive,iel]
        jvl[ive] = ive
    end
    nke = nve

    if ipar > 0
        jvg, jvl = ngls( jvg, jvl, nke )
    end

    for ike = 1:nke
        kdfg[ike] = jvg[ike]
    end
    if ipar == 1
        for  ike = 1:nke
            kdfl[ike] = jvl[ike]
        end
    end

    return kdfg, kdfl

end


# Bubble sort of the arrays KV1 and KV2
function ngls( kv1::Array{Int64}, kv2::Array{Int64}, idim::Int64 )

    bmore = true
    while bmore
        bmore = false
        for icomp = 1:(idim-1)
            if kv1[icomp] > kv1[icomp+1]
                jaux1 = kv1[icomp]
                jaux2 = kv2[icomp]
                kv1[icomp] = kv1[icomp+1]
                kv2[icomp] = kv2[icomp+1]
                kv1[icomp+1] = jaux1
                kv2[icomp+1] = jaux2
                bmore = true
            end
        end
    end

    return kv1, kv2

end


function e011( xi1::Float64, xi2::Float64, djac::Array{Float64,2}, bder::Array{Bool}, dbas::Array{Float64,2}, dhelp::Array{Float64,2}, detj::Float64 )

    const q4 = 0.25

    if( bder[1] )  # function values
        dbas[1,1] = q4*(1.0 - xi1)*(1.0 - xi2)
        dbas[2,1] = q4*(1.0 + xi1)*(1.0 - xi2)
        dbas[3,1] = q4*(1.0 + xi1)*(1.0 + xi2)
        dbas[4,1] = q4*(1.0 - xi1)*(1.0 + xi2)
    end

    if( bder[2] || bder[3] )   # first order derivatives
        dhelp[1,1] = -(1.0 - xi2)
        dhelp[2,1] =  (1.0 - xi2)
        dhelp[3,1] =  (1.0 + xi2)
        dhelp[4,1] = -(1.0 + xi2)
        dhelp[1,2] = -(1.0 - xi1)
        dhelp[2,2] = -(1.0 + xi1)
        dhelp[3,2] =  (1.0 + xi1)
        dhelp[4,2] =  (1.0 - xi1)
        xj = q4/detj

        if( bder[2] )
            dbas[1,2] =  xj*( djac[2,2]*dhelp[1,1] - djac[2,1]*dhelp[1,2] )
            dbas[2,2] =  xj*( djac[2,2]*dhelp[2,1] - djac[2,1]*dhelp[2,2] )
            dbas[3,2] =  xj*( djac[2,2]*dhelp[3,1] - djac[2,1]*dhelp[3,2] )
            dbas[4,2] =  xj*( djac[2,2]*dhelp[4,1] - djac[2,1]*dhelp[4,2] )
        end

        if( bder[3] )
            dbas[1,3] = -xj*( djac[1,2]*dhelp[1,1] - djac[1,1]*dhelp[1,2] )
            dbas[2,3] = -xj*( djac[1,2]*dhelp[2,1] - djac[1,1]*dhelp[2,2] )
            dbas[3,3] = -xj*( djac[1,2]*dhelp[3,1] - djac[1,1]*dhelp[3,2] )
            dbas[4,3] = -xj*( djac[1,2]*dhelp[4,1] - djac[1,1]*dhelp[4,2] )
        end

    end

    return dbas

end
