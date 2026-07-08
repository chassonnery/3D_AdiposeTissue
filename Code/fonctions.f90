module fonctions

implicit none

contains


subroutine init_random_seed(yesorno)
! Initialize a pseudo-random number sequence, with an option for repeatability
	implicit none
	integer, intent(in)                :: yesorno
	integer                            :: i, n, clock
	integer, dimension(:), allocatable :: seed
	
	call Random_seed(size = n)
	allocate(seed(n))

	if (yesorno.eq.0) then
		call System_clock(COUNT=clock)
		seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	else
		seed = 51748307 * (/ (i - 1, i = 1, n) /)
	end if
	
	call Random_seed(PUT = seed)
	deallocate(seed)
end subroutine



subroutine pick_and_reject(unif) 
! génère un vecteur aléatoire choisi selon une loi uniforme dans la sphère unité
	implicit none
	double precision, dimension(3), intent(out) :: unif
	integer :: IsInUnitSphere
	
	IsInUnitSphere = 0
	do while(IsInUnitSphere == 0)
		call random_number(unif(1))
		call random_number(unif(2))
		call random_number(unif(3))
		unif(:) = (unif(:) - 0.5)*2
	
		if (norm(unif)<=1.0 .and. norm(unif)>=0.5) then
			IsInUnitSphere = 1
		end if
		unif(:) = unif(:)*1./norm(unif)
	end do
end subroutine pick_and_reject



function norm(X)
	implicit none
	double precision, dimension(3), intent(in) :: X
	double precision                           :: norm

	norm = sqrt(X(1)**2+X(2)**2+X(3)**2)
end function



subroutine insemination_unif(Xcell,Rcell,dt,time,Pmod,Pdom)
! Randomly inseminate a new cell in the domain :
! 		1- Check if a cell must be inseminated now according to a Poisson process of constant frequency
! 		2- Draw insemination position according to a uniform law in the spacial domain
use definitions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell
	double precision, dimension(:), intent(inout)   :: Rcell
	double precision, intent(in)                    :: dt, time
	type(ParametersMod), intent(inout)              :: Pmod
	type(ParametersDom), intent(in)                 :: Pdom
	double precision :: alea
	
	call random_number(alea)
	if ((alea<1-exp(-Pmod%freq_ens*dt)).and.(Pmod%Ncell.lt.Pmod%Nmax)) then
		Pmod%Ncell = Pmod%Ncell + 1
	
		call random_number(Xcell(Pmod%Ncell,1))
		call random_number(Xcell(Pmod%Ncell,2))
		call random_number(Xcell(Pmod%Ncell,3))
		Xcell(Pmod%Ncell,1) = (2d0*Xcell(Pmod%Ncell,1) - 1)*Pdom%xmaxR
		Xcell(Pmod%Ncell,2) = (2d0*Xcell(Pmod%Ncell,2) - 1)*Pdom%ymaxR
		Xcell(Pmod%Ncell,3) = (2d0*Xcell(Pmod%Ncell,3) - 1)*Pdom%zmaxR
		Rcell(Pmod%Ncell) = Pmod%Rmin
		
		!write(*,'(a,f0.4,a,i0)') "New insemination at time ", time, ", number of cells is now : ", Pmod%Ncell
	end if
end subroutine



subroutine insemination_biased_simple(Xcell,Rcell,dt,time,Pmod,Pdom)
! Randomly inseminate a new cell in the domain :
! 		1- Check if a cell must be inseminated now according to a Poisson process of constant frequency
! 		2- Draw insemination position according to a law proportional to the local cell density power bias
use definitions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell
	double precision, dimension(:), intent(inout)   :: Rcell
	double precision, intent(in)                    :: dt, time
	type(ParametersMod), intent(inout)              :: Pmod
	type(ParametersDom), intent(in)                 :: Pdom
	double precision, dimension(:), allocatable     :: CellDensity
	double precision                                :: alea_time, alea_pos
	double precision                                :: normalisateur, Xnew(3), local_proba
	integer                                         :: success, it, it_max
	
	call random_number(alea_time)
	if ((alea_time<1-exp(-Pmod%freq_ens*dt)).and.(Pmod%Ncell.lt.Pmod%Nmax)) then
		
		allocate(CellDensity(1:(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1)*(2*Pdom%NzPIC+1)))
		call computeCellDensity(CellDensity,Xcell,Rcell,Pmod,Pdom)
		normalisateur = maxval(CellDensity**Pmod%bias_ens)
		
		it = 0
		success = 0
		it_max = 1000
		do while ((it.lt.it_max).and.(success.eq.0))
			call random_number(Xnew(1))
			call random_number(Xnew(2))
			call random_number(Xnew(3))
			Xnew(1) = (2d0*Xnew(1) - 1)*Pdom%xmaxR
			Xnew(2) = (2d0*Xnew(2) - 1)*Pdom%ymaxR
			Xnew(3) = (2d0*Xnew(3) - 1)*Pdom%zmaxR
			
			local_proba = interp3D(CellDensity,Xnew,Pdom)**Pmod%bias_ens / normalisateur
			call random_number(alea_pos)
			if (alea_pos.lt.local_proba) then
				success = 1
				Pmod%Ncell = Pmod%Ncell + 1
				Xcell(Pmod%Ncell,:) = Xnew
				Rcell(Pmod%Ncell) = Pmod%Rmin
				!write(*,'(a,f0.4,a,i0)') "New insemination at time ", time, ", number of cells is now : ", Pmod%Ncell
			end if
			
			it = it+1
		end do
		
		!if (success.eq.0) write(*,'("Warning : insemination failed at time ",f0.4)') time 
		deallocate(CellDensity)
	end if
end subroutine



subroutine insemination_biased(Xcell,Rcell,dt,time,Pmod,Pdom)
! Randomly inseminate new cell(s) in the domain :
!		Draw random position 
!		Check if a cell must be inseminated now at that position according to a law proportional to the temporal insemination frequency and the local cell density power b, divided by it_max
!		Repeat the process it_max times to be sure the domain is correctly sampled
use definitions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell
	double precision, dimension(:), intent(inout)   :: Rcell
	double precision, intent(in)                    :: dt, time
	type(ParametersMod), intent(inout)              :: Pmod
	type(ParametersDom), intent(in)                 :: Pdom
	double precision, dimension(:), allocatable     :: CellDensity
	double precision                                :: alea
	double precision                                :: normalisateur, Xnew(3), local_proba
	integer                                         :: it, it_max
	
	it_max = 100 * 8*Pdom%NxPIC*Pdom%NyPIC*Pdom%NzPIC
	
	allocate(CellDensity(1:(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1)*(2*Pdom%NzPIC+1)))
	call computeCellDensity(CellDensity,Xcell,Rcell,Pmod,Pdom)
	normalisateur = integral3D(CellDensity,Pmod%bias_ens,it_max,Pdom)/(8d0*Pdom%xmaxR*Pdom%ymaxR*Pdom%zmaxR)
	
	do it=1,it_max
		call random_number(Xnew(1))
		call random_number(Xnew(2))
		call random_number(Xnew(3))
		Xnew(1) = (2d0*Xnew(1) - 1)*Pdom%xmaxR
		Xnew(2) = (2d0*Xnew(2) - 1)*Pdom%ymaxR
		Xnew(3) = (2d0*Xnew(3) - 1)*Pdom%zmaxR
		
		local_proba = (1d0/it_max) * Pmod%freq_ens * (interp3D(CellDensity,Xnew,Pdom)**Pmod%bias_ens)/normalisateur
		call random_number(alea)
		if ((alea<1-exp(-local_proba*dt)).and.(Pmod%Ncell.lt.Pmod%Nmax)) then
			Pmod%Ncell = Pmod%Ncell + 1
			Xcell(Pmod%Ncell,:) = Xnew
			Rcell(Pmod%Ncell) = Pmod%Rmin
			!write(*,'(a,f0.4,a,i0)') "New insemination at time ", time, ", number of cells is now : ", Pmod%Ncell
			
			if (Pmod%Ncell.eq.Pmod%Nmax) exit
		end if
	end do
	
	deallocate(CellDensity)
end subroutine



subroutine computeCellDensity(CellDensity,Xcell,Rcell,Pmod,Pdom)
! Compute the density of cells at each node of a cuboid grid
use definitions
	implicit none
	double precision, dimension(:), intent(out)  :: CellDensity
	double precision, dimension(:,:), intent(in) :: Xcell
	double precision, dimension(:), intent(in)   :: Rcell
	type(ParametersMod), intent(in) :: Pmod
	type(ParametersDom), intent(in) :: Pdom
	double precision                :: Xnode(3), weight, Vc
	integer :: i,ix,iy,iz,lx,ly,lz,l
	double precision, parameter :: pi=3.141592653589793238
	
	CellDensity(:) = 0.0
	do i=1,Pmod%Ncell
		! Compute cell volume
		Vc = (4.0/3)*pi*Rcell(i)**3
		
		! find indexes of the "bottom left" corner of the box containing agent i
		ix = floor( (Xcell(i,1)+Pdom%xmaxR)/Pdom%dxPIC ) ! use fonction floor for the truncature
		iy = floor( (Xcell(i,2)+Pdom%ymaxR)/Pdom%dyPIC )
		iz = floor( (Xcell(i,3)+Pdom%zmaxR)/Pdom%dzPIC )
		
		! Iterate over all corners of that box (= all nodes of the map surrounding the cell)
		do lx = 0,1
			do ly = 0,1
				do lz = 0,1
					! Weight of the cell over this node = volume of the cuboid space between the cell and the opposite node
					Xnode(1) = (ix+1-lx)*Pdom%dxPIC - Pdom%xmaxR
					Xnode(2) = (iy+1-ly)*Pdom%dyPIC - Pdom%ymaxR
					Xnode(3) = (iz+1-lz)*Pdom%dzPIC - Pdom%zmaxR
					weight = abs((Xnode(1)-Xcell(i,1)) * (Xnode(2)-Xcell(i,2)) * (Xnode(3)-Xcell(i,3)))
					
					! Add cell contribution to density map
					l = 1 + ix + lx + (iy+ly)*(2*Pdom%NxPIC+1) + (iz+lz)*(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1)
					CellDensity(l) = CellDensity(l) + Vc*weight/((Pdom%dxPIC*Pdom%dyPIC*Pdom%dzPIC)**2)
				end do
			end do
		end do
	end do
	
	call HandleDensityBoundaryCondition(CellDensity,Pmod,Pdom)
end subroutine



subroutine computeFiberDensity(FibDensity,Xfib,omega,Pmod,Pdom,Box)
! Compute the density of fibers at each node of a cuboid grid. 
! Each fiber is divided into Ndiv sections of equal length Ddiv, with Ndiv = floor(Lf/(2*Rf)) + 1 so that the sections are as close as possible to spheres of radius Rf.
use definitions
	implicit none
	double precision, dimension(:), intent(out)  :: FibDensity
	double precision, dimension(:,:), intent(in) :: Xfib, omega
	type(ParametersMod), intent(in) :: Pmod
	type(ParametersDom), intent(in) :: Pdom
	type(ParametersBox), intent(in) :: Box
	double precision                :: Xnode(3), Xtemp(3), weight, Ddiv, Vfib_div
	integer :: Ndiv
	integer :: i,j,k,ix,iy,iz,lx,ly,lz,l
	double precision, parameter :: pi=3.141592653589793238
	
	Ndiv = floor(Pmod%Lf/(2*Pmod%Rf)) + 1
	Ddiv = Pmod%Lf/(Ndiv - 1) ! ≥ 2*Pmod%Rf. Should be close to 2*Pmod%Rf for a good approximation
	Vfib_div = ( (4.0/3)*pi*Pmod%Rf**3 + pi*Pmod%Lf*Pmod%Rf**2 )/Ndiv
	
	FibDensity(:) = 0.0
	do i=1,Pmod%Nfib
		do j=0,Ndiv-1
			! Find position of the jth section of fiber i
			Xtemp = Xfib(i,:) + (j*Ddiv - Pmod%Lf/2d0)*omega(i,:)
			
			call findBox(Xtemp,k,Pdom)
			if (Box%BoxType(k,1).eq.1) then ! apply periodic boundary condition to this section of the fiber
				Xtemp(:) = Xtemp(:) - Box%BoxType(k,2:4)
			elseif (Box%BoxType(k,1).gt.1) then
				if ((abs(Xtemp(1)).gt.Pdom%xmaxR+Pdom%dxPIC).or.(abs(Xtemp(2)).gt.Pdom%ymaxR+Pdom%dyPIC) &
				.or.(abs(Xtemp(3)).gt.Pdom%zmaxR+Pdom%dzPIC)) then
					cycle ! ignore sections that are FAR outside the domain
				end if
			end if
			
			! find indexes of the "bottom left" corner of the box containing section j of fiber i
			! /!\ to accomodate for "dangling" sections located outside the domain, the indexes ix,iy,iz range from -1 to 2*Pdom%NxPIC+1 instead of 0 to 2*Pdom%NxPIC
			ix = floor( (Xtemp(1)+Pdom%xmaxR+Pdom%dxPIC)/Pdom%dxPIC ) - 1 ! use fonction floor for the truncature
			iy = floor( (Xtemp(2)+Pdom%ymaxR+Pdom%dyPIC)/Pdom%dyPIC ) - 1
			iz = floor( (Xtemp(3)+Pdom%zmaxR+Pdom%dzPIC)/Pdom%dzPIC ) - 1
			
			! Iterate over all corners of that box (= all nodes of the map surrounding the concerned section)
			do lx = 0,1
				do ly = 0,1
					do lz = 0,1
						if ((ix+lx.ge.0).and.(ix+lx.le.2*Pdom%NxPIC).and.(iy+ly.ge.0).and.(iy+ly.le.2*Pdom%NyPIC) &
						.and.(iz+lz.ge.0).and.(iz+lz.le.2*Pdom%NzPIC)) then
							! Weight of the cell over this node = volume of the space between the section center and the opposite node
							Xnode(1) = (ix+1-lx)*Pdom%dxPIC - Pdom%xmaxR
							Xnode(2) = (iy+1-ly)*Pdom%dyPIC - Pdom%ymaxR
							Xnode(3) = (iz+1-lz)*Pdom%dzPIC - Pdom%zmaxR
							weight = abs((Xnode(1)-Xtemp(1)) * (Xnode(2)-Xtemp(2)) * (Xnode(3)-Xtemp(3)))
							
							! Add contribution to density map
							l = 1 + ix + lx + (iy+ly)*(2*Pdom%NxPIC+1) + (iz+lz)*(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1)
							FibDensity(l) = FibDensity(l) + Vfib_div*weight/((Pdom%dxPIC*Pdom%dyPIC*Pdom%dzPIC)**2)
						end if
					end do
				end do
			end do
		end do
	end do
	
	call HandleDensityBoundaryCondition(FibDensity,Pmod,Pdom)
end subroutine



subroutine HandleDensityBoundaryCondition(Density,Pmod,Pdom)
! Apply periodic boundary coundition to the borders of a density map
use definitions
	implicit none
	double precision, dimension(:), intent(out)  :: Density
	type(ParametersMod), intent(in) :: Pmod
	type(ParametersDom), intent(in) :: Pdom
	double precision                :: sum
	integer :: ix,iy,iz,l
	
	if ( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) ) then
		do iy = 0,2*Pdom%NyPIC
			do iz = 0,2*Pdom%NzPIC
				l = 1 + iy*(2*Pdom%NxPIC+1) + iz*(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1)
				sum = Density(l+0) + Density(l+2*Pdom%NxPIC)
				Density(l+0) = sum
				Density(l+2*Pdom%NxPIC) = sum
			end do
		end do
	end if
	
	if ( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) ) then
		do ix = 0,2*Pdom%NxPIC
			do iz = 0,2*Pdom%NzPIC
				l = 1 + ix + iz*(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1)
				sum = Density(l+0) + Density(l+2*Pdom%NyPIC*(2*Pdom%NxPIC+1))
				Density(l+0) = sum
				Density(l+2*Pdom%NyPIC*(2*Pdom%NxPIC+1)) = sum
			end do
		end do
	end if

	if ( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) ) then
		do ix = 0,2*Pdom%NxPIC
			do iy = 0,2*Pdom%NyPIC
				l = 1 + ix + iy*(2*Pdom%NxPIC+1)
				sum = Density(l+0) + Density(l+2*Pdom%NzPIC*(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1))
				Density(l+0) = sum
				Density(l+2*Pdom%NzPIC*(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1)) = sum
			end do
		end do
	end if
end subroutine



function interp3D(DensityMap,X,Pdom)
! Compute the approximate density at point X using a trilinear interpolation of the 3D grid DensityMap
use definitions
	implicit none
	double precision                           :: interp3D
	double precision, dimension(:), intent(in) :: DensityMap
	double precision, dimension(3), intent(in) :: X
	type(ParametersDom), intent(in) :: Pdom
	double precision                :: Xnode(3), weight
	integer :: ix,iy,iz,lx,ly,lz,l
	
	! Check that point X is in range
	if (abs(X(1)).gt.Pdom%xmaxR) then
		write(*,'(a)') "Trying to interpolate density at a point outside the domain ! (function inter3D, axis X)"
		STOP
	end if
	if (abs(X(2)).gt.Pdom%ymaxR) then
		write(*,'(a)') "Trying to interpolate density at a point outside the domain ! (function inter3D, axis Y)"
		STOP
	end if
	if (abs(X(3)).gt.Pdom%zmaxR) then
		write(*,'(a)') "Trying to interpolate density at a point outside the domain ! (function inter3D, axis Z)"
		STOP
	end if
	
	! find indexes of the "bottom left" corner of the box containing point X
	ix = floor( (X(1)+Pdom%xmaxR)/Pdom%dxPIC ) ! use fonction floor for the truncature
	iy = floor( (X(2)+Pdom%ymaxR)/Pdom%dyPIC )
	iz = floor( (X(3)+Pdom%zmaxR)/Pdom%dzPIC )
			
	! Iterate over all corners of that box (= all nodes of the map surrounding point X)
	interp3D = 0.0
	do lx = 0,1
		do ly = 0,1
			do lz = 0,1
				! Weight of this node over point X = volume of the space between point X and the opposite node
				Xnode(1) = (ix+1-lx)*Pdom%dxPIC - Pdom%xmaxR
				Xnode(2) = (iy+1-ly)*Pdom%dyPIC - Pdom%ymaxR
				Xnode(3) = (iz+1-lz)*Pdom%dzPIC - Pdom%zmaxR
				weight = abs((Xnode(1)-X(1)) * (Xnode(2)-X(2)) * (Xnode(3)-X(3)))
				
				! Add node contribution to density at point X
				l = 1 + ix + lx + (iy+ly)*(2*Pdom%NxPIC+1) + (iz+lz)*(2*Pdom%NxPIC+1)*(2*Pdom%NyPIC+1)
				interp3D = interp3D + DensityMap(l)*weight/(Pdom%dxPIC*Pdom%dyPIC*Pdom%dzPIC)
			end do
		end do
	end do
end function



function integral3D(Map,power,it_max,Pdom)
! Compute the sum (approximated integrale) of Map while accounting for boundary conditions
use definitions
	implicit none
	double precision                           :: integral3D
	double precision, dimension(:), intent(in) :: Map
	double precision, intent(in)               :: power
	integer, intent(in)                        :: it_max
	type(ParametersDom), intent(in) :: Pdom
	double precision :: X(3), Nperunit, dx, dy, dz
	integer :: i, j, k, Nx, Ny, Nz
	
	Nperunit = max(Pdom%xmaxR,Pdom%ymaxR,Pdom%zmaxR) * ( it_max / (Pdom%xmaxR*Pdom%ymaxR*Pdom%zmaxR) )**(1/3d0)
	Nx = int( Nperunit*Pdom%xmaxR/max(Pdom%xmaxR,Pdom%ymaxR,Pdom%zmaxR) )
	Ny = int( Nperunit*Pdom%ymaxR/max(Pdom%xmaxR,Pdom%ymaxR,Pdom%zmaxR) )
	Nz = int( Nperunit*Pdom%zmaxR/max(Pdom%xmaxR,Pdom%ymaxR,Pdom%zmaxR) )
	dx = 2d0*Pdom%xmaxR/Nx
	dy = 2d0*Pdom%ymaxR/Ny
	dz = 2d0*Pdom%zmaxR/Nz
	
	integral3D = 0
	do i=0,Nx-1
		do j=0,Ny-1
			do k=0,Nz-1
				X(1) = (i+0.5)*dx - Pdom%xmaxR 
				X(2) = (j+0.5)*dy - Pdom%ymaxR 
				X(3) = (k+0.5)*dz - Pdom%zmaxR
				integral3D = integral3D + interp3D(Map,X,Pdom)**power *dx*dy*dz
			end do
		end do
	end do
end function



subroutine find_overlap_fibfib(Xi,Xk,wi,wk,Pmod,dik,lik,lki)
! Find the closest points segments and rays of two fibers
! http://geomalgorithms.com/a07-_distance.html
use definitions
	implicit none
	double precision, dimension(3), intent(in) :: Xi, Xk, wi, wk
	type(ParametersMod), intent(in)            :: Pmod
	double precision                           :: prod,s,t
	double precision, dimension(3)             :: x
	double precision, intent(out)              :: dik,lik,lki

	x = Xi-Xk
	prod = (dot_product(wi,wi)*dot_product(wk,wk) - dot_product(wi,wk)**2)
	
	if (prod < 0.001) then ! lines almost parallel
		s = 0
		t = dot_product(wk,x)/dot_product(wk,wk)
		!! alternative : s and t at equidistance from the two center points Xi and Xk
		! s = dot_product(wi,-x)/2d0
		! t = dot_product(wk, x)/2d0
	else
		s = (dot_product(wi,wk)*dot_product(wk,x) - dot_product(wi,x))/prod
		t = (dot_product(wk,x) - dot_product(wi,wk)*dot_product(wi,x))/prod
		
	endif
	
	! Cut-off to range [-Lf/2,Lf/2]
	if (abs(s) > Pmod%Lf/2d0) then
		s = sign(Pmod%Lf/2d0,s) ! sign(a,b) = value of a with the sign of b
		t = dot_product(wk,x) + s*dot_product(wk,wi)
	end if
	
	if (abs(t) > Pmod%Lf/2d0) then
		t = sign(Pmod%Lf/2d0,t)
		s = - dot_product(wi,x) + t*dot_product(wi,wk)
		if (abs(s) > Pmod%Lf/2d0) then
			s = sign(Pmod%Lf/2d0,s)
		end if
	end if
	
	dik = norm(Xi + s*wi - (Xk + t*wk))
	lik = s
	lki = t
end subroutine find_overlap_fibfib



subroutine find_overlap_cellfib(Xc,Xf,wf,Pmod,d,t)
! Find point on fiber i closest to the center of cell k
! http://geomalgorithms.com/a07-_distance.html
use definitions
    implicit none
	double precision, dimension(3), intent(in) :: Xc, Xf, wf
	type(ParametersMod), intent(in)            :: Pmod
    double precision, intent(out)              :: d,t

	t = - dot_product(wf,Xf-Xc)/dot_product(wf,wf)
	
	if (t <= -Pmod%Lf/2d0) then
	    t = -Pmod%Lf/2d0
	elseif (t >= Pmod%Lf/2d0) then
	    t = Pmod%Lf/2d0
	endif
	
	d = norm(Xf + t*wf - Xc)
end subroutine find_overlap_cellfib



subroutine grad_cellcell(grad,Xi,Xk,Ri,Rk,Pmod)
use definitions
	implicit none
	double precision, dimension(3), intent(inout) :: grad
	double precision, dimension(3), intent(in) :: Xi,Xk
	double precision, intent(in)               :: Ri,Rk
	type(ParametersMod), intent(in)            :: Pmod
	double precision                           :: dik,Reff
	double precision, dimension(3)             :: Vik

	
	dik = norm(Xi - Xk)
	if (dik<Ri+Rk) then ! cells i and k overlap
		if (dik<0.001) then
			call pick_and_reject(Vik)
		else
			Vik = (Xi - Xk)/dik
		end if
		
		! Hertzian contact force
		Reff = Ri*Rk*1d0/(Ri + Rk)
		grad(:) = grad(:) + Pmod%alpha_repCC * (4d0/3) * Reff**(0.5) * (Ri+Rk-dik)**(1.5) * Vik(:)
	end if	
end subroutine grad_cellcell



subroutine grad_cellfib(grad,Xi,Xk,Ri,omegak,Pmod)
! Contrary to the function grad_cellcell and grad_fibfib, which are antisymetric with respect to arguments i and k (they return the force/torque exerted by agent k on agent i), in this function the arguments i and k CAN NOT be reverted. 
! This function always return the force and torque exerted by the cell i on the fiber k. To obtain the force exerted by the fiber k on the cell i, call grad_cellfib(gradtemp,Xi,Xk,Ri,omegak,Pmod) then -gradtemp(1:3)
use definitions
	implicit none
	double precision, dimension(6), intent(inout) :: grad
	double precision, dimension(3), intent(in) :: Xi,Xk,omegak
	double precision, intent(in)               :: Ri
	type(ParametersMod), intent(in)            :: Pmod
	double precision                           :: dik,lki,Reff
	double precision, dimension(3)             :: Vki,torque_rep

	
	call find_overlap_cellfib(Xi,Xk,omegak,Pmod,dik,lki)

	if (dik < Ri+Pmod%Rf) then ! cell i and fiber k overlap
		if (dik<0.001) then
			call pick_and_reject(Vki)
		else
			Vki = (Xk + lki*omegak - Xi)/dik
		end if
		
		torque_rep(1) = lki*( Vki(1)*(omegak(2)**2 + omegak(3)**2) - Vki(2)*omegak(1)*omegak(2) - Vki(3)*omegak(1)*omegak(3) )  
		torque_rep(2) = lki*( Vki(2)*(omegak(1)**2 + omegak(3)**2) - Vki(1)*omegak(1)*omegak(2) - Vki(3)*omegak(2)*omegak(3) ) 
		torque_rep(3) = lki*( Vki(3)*(omegak(1)**2 + omegak(2)**2) - Vki(1)*omegak(1)*omegak(3) - Vki(2)*omegak(2)*omegak(3) ) 
		
		! Hertzian contact force
		Reff = Ri*Pmod%Rf*1d0/(Ri + Pmod%Rf)
		grad(1:3) = grad(1:3) + Pmod%alpha_repCF * (4d0/3) * Reff**(0.5) * (Ri+Pmod%Rf-dik)**(1.5) * Vki(:)
		grad(4:6) = grad(4:6) + Pmod%alpha_repCF * (4d0/3) * Reff**(0.5) * (Ri+Pmod%Rf-dik)**(1.5) * torque_rep(:)
	end if
end subroutine grad_cellfib



subroutine grad_fibfib(grad,Xi,Xk,omegai,omegak,i,k,Pmod,link,distpoints)
use definitions
	implicit none
	double precision, dimension(6), intent(inout) :: grad
	integer, dimension(:,:), intent(in)           :: link
	double precision, dimension(:,:), intent(in)  :: distpoints
	double precision, dimension(3), intent(in)    :: Xi,Xk,omegai,omegak
	type(ParametersMod), intent(in)               :: Pmod
	integer, intent(in)                           :: i,k
	double precision                              :: dik, lik, lki, c
	double precision, dimension(3)                :: Vik,vpoints,torque_rep,torque_rappel,torque_align,om,v
	double precision, dimension(3,3)              :: R,Rot
	double precision, parameter                   :: pi=3.141592653589793238

	
	call find_overlap_fibfib(Xi,Xk,omegai,omegak,Pmod,dik,lik,lki)
	! reference : Interparticle torques suppress motility-induced phase separation for rodlike particles, Damme, Rondenburg, 2019, arXiv

	if ( (link(i,k)==0).or.(Pmod%cumul_replink.eq.0) ) then
		if (dik<2*Pmod%Rf) then ! i and k overlap
			vpoints = Xi + lik*omegai - (Xk + lki*omegak)
			if (norm(Vpoints)<0.001) then
				call pick_and_reject(Vik)
!				write(*,'(a,i0,a,i0)') "Warning : degenerated interaction between fibers ", i, " and ", k
			else
				Vik = vpoints/norm(vpoints)
			end if
			
			torque_rep(1) = lik*( Vik(1)*(omegai(2)**2 + omegai(3)**2) - Vik(2)*omegai(1)*omegai(2) - Vik(3)*omegai(1)*omegai(3) )  
			torque_rep(2) = lik*( Vik(2)*(omegai(1)**2 + omegai(3)**2) - Vik(1)*omegai(1)*omegai(2) - Vik(3)*omegai(2)*omegai(3) ) 
			torque_rep(3) = lik*( Vik(3)*(omegai(1)**2 + omegai(2)**2) - Vik(1)*omegai(1)*omegai(3) - Vik(2)*omegai(2)*omegai(3) ) 
			
			! Hertzian contact force between crossed cylinders of equal radius
			grad(1:3) = grad(1:3) + Pmod%alpha_repFF * (4d0/3) * Pmod%Rf**(0.5) * (2.0*Pmod%Rf-dik)**(1.5) * Vik(:)
			grad(4:6) = grad(4:6) + Pmod%alpha_repFF * (4d0/3) * Pmod%Rf**(0.5) * (2.0*Pmod%Rf-dik)**(1.5) * torque_rep(:)
		end if
	end if
		
	if (link(i,k)==1) then
		! Elastic restoring force following Hooke law
		vpoints = Xi + distpoints(i,k)*omegai - (Xk + distpoints(k,i)*omegak)
		Vik = (Pmod%d0 - norm(vpoints))*vpoints/norm(vpoints)
		
		torque_rappel(1) = distpoints(i,k)*( Vik(1)*(omegai(2)**2 + omegai(3)**2) &
		                   - Vik(2)*omegai(1)*omegai(2) - Vik(3)*omegai(1)*omegai(3) )  
		torque_rappel(2) = distpoints(i,k)*( Vik(2)*(omegai(1)**2 + omegai(3)**2) &
		                   - Vik(1)*omegai(1)*omegai(2) - Vik(3)*omegai(2)*omegai(3) ) 
		torque_rappel(3) = distpoints(i,k)*( Vik(3)*(omegai(1)**2 + omegai(2)**2) &
		                   - Vik(1)*omegai(1)*omegai(3) - Vik(2)*omegai(2)*omegai(3) ) 
		
		! NEMATIC alignment force
		c = dot_product(omegai,omegak)
		if (acos(c) < pi/2d0) then
			v = omegak
		else
			v = -omegak
		endif
		c = dot_product(omegai,v)
		
		om(1) = omegai(2)*v(3) - omegai(3)*v(2)
		om(2) = omegai(3)*v(1) - omegai(1)*v(3)
		om(3) = omegai(1)*v(2) - omegai(2)*v(1)
		Rot=0
		if (abs(1.0 - c**2)>0.01) then
			R = 0
			R(1,2) = -om(3)
			R(1,3) =  om(2)
			R(2,1) =  om(3)
			R(2,3) = -om(1)
			R(3,1) = -om(2)
			R(3,2) =  om(1)
			Rot = R + (1-c)/(norm(om)**2) * matmul(R,R)
			! Rot+Identity = rotation of axis c/||c|| and angle arccos(c)
			! So (Rot+Identity).omegai = v and the variation needed to bring omegai colinear to v is (v-omegai) = Rot.omegai
		endif
		
		torque_align(1) = Rot(1,1)*omegai(1) + Rot(1,2)*omegai(2) + Rot(1,3)*omegai(3)
		torque_align(2) = Rot(2,1)*omegai(1) + Rot(2,2)*omegai(2) + Rot(2,3)*omegai(3)
		torque_align(3) = Rot(3,1)*omegai(1) + Rot(3,2)*omegai(2) + Rot(3,3)*omegai(3)
		
		grad(1) = grad(1) + Pmod%alpha_rappel*Vik(1)
		grad(2) = grad(2) + Pmod%alpha_rappel*Vik(2)
		grad(3) = grad(3) + Pmod%alpha_rappel*Vik(3)
		grad(4) = grad(4) + Pmod%alpha_rappel*torque_rappel(1) + Pmod%alpha_align*torque_align(1)
		grad(5) = grad(5) + Pmod%alpha_rappel*torque_rappel(2) + Pmod%alpha_align*torque_align(2)
		grad(6) = grad(6) + Pmod%alpha_rappel*torque_rappel(3) + Pmod%alpha_align*torque_align(3)
	endif
end subroutine grad_fibfib



subroutine computeBoxType(Pdom,Pmod,Box)
use definitions
    implicit none
	type(ParametersDom),intent(in)    :: Pdom
	type(ParametersMod),intent(in)    :: Pmod
	type(ParametersBox),intent(inout) :: Box
	integer                           :: ix,iy,iz,l
		
	Box%BoxType(:,:) = 0.0
	! BoxType(i,1) = 0 --> real box
	!				 1 --> periodic ghost box
	!				 2 --> must-be-empty box (mixed border condition)
	!				 3 --> dirichlet box
	
	! Types are computed by increasing order of priority : if a box is on the edge between a dirichlet (3) or must-be-empty (2) border and a periodic border (1), then it is of type periodic (1). If a box is on the edge between a must-be-empty (2) and a dirichlet border (3), then it is of type must-be-empty (2).
	
	do ix = 0,2*Pdom%Nx+1
		do iy = 0,2*Pdom%Ny+1
			do iz = 0,2*Pdom%Nz+1
				l = 1 + ix + 2*iy*(Pdom%Nx+1) + 4*iz*(Pdom%Nx+1)*(Pdom%Ny+1)
				
				! Dirichlet boxes
				if ( (ix.eq.0).and.(Pmod%DirLim(1).eq.1) ) then
					Box%BoxType(l,1) = 3
				elseif ( (ix.eq.2*Pdom%Nx+1).and.(Pmod%DirLim(2).eq.1) ) then
					Box%BoxType(l,1) = 3
				end if
				
				if ( (iy.eq.0).and.(Pmod%DirLim(3).eq.1) ) then
					Box%BoxType(l,1) = 3
				elseif ( (iy.eq.2*Pdom%Ny+1).and.(Pmod%DirLim(4).eq.1) ) then
					Box%BoxType(l,1) = 3
				end if
				
				if ( (iz.eq.0).and.(Pmod%DirLim(5).eq.1) ) then
					Box%BoxType(l,1) = 3
				elseif ( (iz.eq.2*Pdom%Nz+1).and.(Pmod%DirLim(6).eq.1) ) then
					Box%BoxType(l,1) = 3
				end if

				! Must-be-empty boxes
				if ( (ix.eq.0).and.(Pmod%DirLim(1).eq.0) ) then
					Box%BoxType(l,1) = 2
				elseif ( (ix.eq.2*Pdom%Nx+1).and.(Pmod%DirLim(2).eq.0) ) then
					Box%BoxType(l,1) = 2
				end if
				if ( (iy.eq.0).and.(Pmod%DirLim(3).eq.0) ) then
					Box%BoxType(l,1) = 2
				elseif ( (iy.eq.2*Pdom%Ny+1).and.(Pmod%DirLim(4).eq.0) ) then
					Box%BoxType(l,1) = 2
				end if
				if ( (iz.eq.0).and.(Pmod%DirLim(5).eq.0) ) then
					Box%BoxType(l,1) = 2
				elseif ( (iz.eq.2*Pdom%Nz+1).and.(Pmod%DirLim(6).eq.0) ) then
					Box%BoxType(l,1) = 2
				end if
				
				! Periodic boxes
				if ( (ix.eq.0).and.( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,2) = -2*Pdom%xmaxR
				else if ( (ix.eq.2*Pdom%Nx+1).and.( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,2) = 2*Pdom%xmaxR
				end if
				
				if ( (iy.eq.0).and.( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,3) = -2*Pdom%ymaxR
				else if ( (iy.eq.2*Pdom%Ny+1).and.( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,3) = 2*Pdom%ymaxR
				end if
				
				if ( (iz.eq.0).and.( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,4) = -2*Pdom%zmaxR
				else if ( (iz.eq.2*Pdom%Nz+1).and.( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) ) ) then
					Box%BoxType(l,1) = 1
					Box%BoxType(l,4) = 2*Pdom%zmaxR
				end if
				
			end do
		end do
	end do
end subroutine computeBoxType



subroutine findBox(X,k,Pdom)
! Computes the index k of the box in which point X is located
! NB: this function INCLUDES THE CL (box 1 = left bottom with CL)
use definitions
	implicit none
	double precision, dimension(:), intent(in) :: X
	integer, intent(out)                       :: k
	type(ParametersDom), intent(in)            :: Pdom
	integer                                    :: ix,iy,iz
	
	ix = floor( (X(1)+Pdom%xmax)/Pdom%dx ) ! use fonction floor for the truncature
	iy = floor( (X(2)+Pdom%ymax)/Pdom%dy )
	iz = floor( (X(3)+Pdom%zmax)/Pdom%dz )
	
	k = 1 + ix + 2*iy*(Pdom%Nx+1) + 4*iz*(Pdom%Nx+1)*(Pdom%Ny+1)
end subroutine findBox



subroutine UpdateVerlet(i,X,Pdom,Box)
use definitions
	implicit none
	integer, intent(in)                       :: i
	double precision, dimension(3),intent(in) :: X
	type(ParametersDom),intent(in)            :: Pdom
	type(ParametersBox),intent(inout)         :: Box
	integer                                   :: k

	call findBox(X,k,Pdom)
	
	Box%AgentBox(i) = k
	
	if (Box%FirstAgent(k).eq.0) then ! The box is still empty, begin it with i
		Box%FirstAgent(k) = i 
		Box%LastAgent(k)  = i
	else ! add index i to the end of verletlist and change last agent
		Box%VerletList(Box%LastAgent(k)) = i
		Box%LastAgent(k) = i
	end if
end subroutine UpdateVerlet



subroutine cleanInsideDomain(Pdom,Box)
! Function which clean the array FirstAgent and LastAgent for all boxes *inside* the domain 
use definitions
	implicit none
	type(ParametersDom), intent(in)    :: Pdom
	type(ParametersBox), intent(inout) :: Box
	integer                            :: k,ix,iy,iz
	
	do ix = 1,2*Pdom%Nx
		do iy = 1,2*Pdom%Ny
			do iz = 1,2*Pdom%Nz
				k = 1 + ix + 2*iy*(Pdom%Nx+1) + 4*iz*(Pdom%Nx+1)*(Pdom%Ny+1)
				Box%FirstAgent(k) = 0
				Box%LastAgent(k) = 0
			end do
		end do
	end do
end subroutine cleanInsideDomain



subroutine FillGhostBoxes(Pmod,Pdom,Box)
! Function which fill the "ghost boxes" on the periodic borders of the domain
use definitions
	implicit none
	type(ParametersMod), intent(in)    :: Pmod
	type(ParametersDom), intent(in)    :: Pdom
	type(ParametersBox), intent(inout) :: Box
	integer                            :: k,kreal
	integer                            :: ix,iy,iz,ixreal,iyreal,izreal
	
	do ix = 0,2*Pdom%Nx+1
		do iy = 0,2*Pdom%Ny+1
			do iz = 0,2*Pdom%Nz+1
				
				k = 1 + ix + 2*iy*(Pdom%Nx+1) + 4*iz*(Pdom%Nx+1)*(Pdom%Ny+1)
				! If this box is a "ghost" one, then its first agent is equal to the first agent of its associated "real box" or "Dirichlet box"
				if (Box%BoxType(k,1).eq.1) then
					
					if ( (ix.eq.0).and.( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) )) then 
						ixreal = 2*Pdom%Nx
					elseif ( (ix.eq.2*Pdom%Nx+1).and.( (Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.0) )) then
						ixreal = 1
					else
						ixreal = ix
					end if
					
					if ( (iy.eq.0).and.( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) )) then 
						iyreal = 2*Pdom%Ny
					elseif ( (iy.eq.2*Pdom%Ny+1).and.( (Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.0) )) then
						iyreal = 1
					else
						iyreal = iy
					end if
					
					if ( (iz.eq.0).and.( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) )) then 
						izreal = 2*Pdom%Nz
					elseif ( (iz.eq.2*Pdom%Nz+1).and.( (Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.0) )) then
						izreal = 1
					else
						izreal = iz
					end if
				
					kreal = 1 + ixreal + 2*iyreal*(Pdom%Nx+1) + 4*izreal*(Pdom%Nx+1)*(Pdom%Ny+1)
					Box%FirstAgent(k) = Box%FirstAgent(kreal)
					Box%LastAgent(k) = Box%LastAgent(kreal)
					
					if (kreal.eq.k) then
						write(*,'(a,i0,a)') "Erreur : BoxType of box ", k, " is badly defined"
						STOP
					end if
				end if
				
			end do
		end do
	end do
end subroutine FillGhostBoxes



subroutine HandleCellBoundaryCondition(X,Pmod,Pdom,Box)
use definitions
	implicit none
	double precision, dimension(3), intent(inout) :: X
	type(ParametersMod), intent(in) :: Pmod
	type(ParametersDom), intent(in) :: Pdom
	type(ParametersBox), intent(in) :: Box
	integer                         :: ix, iy, iz, k

	call findBox(X,k,Pdom)
	
	! If a cell exit the computation domain by entering a periodic ghost box : make it come back from the opposite side.
	if (Box%BoxType(k,1).eq.1) then
		X(:) = X(:) - Box%BoxType(k,2:4)
		call findBox(X,k,Pdom)
		if (Box%BoxType(k,1).eq.1) then
			write(*,'(a)') "Error : cell is still in a periodic box after applying boundary condition."
			STOP
		end if
	end if
	
	! For the cells, dirichlet boxes are the same as must-be-empty boxes.
	! If a cell exit the computation domain by entering a must-be-empty box : bring it back inside from the same side.
	if ((Box%BoxType(k,1).eq.2).or.(Box%BoxType(k,1).eq.3)) then
		ix = floor( (X(1)+Pdom%xmax)/Pdom%dx )
		iy = floor( (X(2)+Pdom%ymax)/Pdom%dy )
		iz = floor( (X(3)+Pdom%zmax)/Pdom%dz )
		
		if ( ((ix.eq.0).or.(ix.eq.2*Pdom%Nx+1)) .and. ((Pmod%DirLim(1).eq.1).or.(Pmod%DirLim(2).eq.1)) ) then
			X(1) = sign(1d0,X(1))*0.999*Pdom%xmaxR
		end if
		
		if ( ((iy.eq.0).or.(iy.eq.2*Pdom%Ny+1)) .and. ((Pmod%DirLim(3).eq.1).or.(Pmod%DirLim(4).eq.1)) ) then
			X(2) = sign(1d0,X(2))*0.999*Pdom%ymaxR
		end if
		
		if ( ((iz.eq.0).or.(iz.eq.2*Pdom%Nz+1)) .and. ((Pmod%DirLim(5).eq.1).or.(Pmod%DirLim(6).eq.1)) ) then
			X(3) = sign(1d0,X(3))*0.999*Pdom%zmaxR
		end if
		
		call findBox(X,k,Pdom)
		if ((Box%BoxType(k,1).eq.2).or.(Box%BoxType(k,1).eq.3)) then
			write(*,'(a)') "Error : cell is still in a must-be-empty or dirichlet box after applying boundary condition."
			STOP
		end if
	end if
end subroutine



subroutine HandleFibBoundaryCondition(X,IsInDirichlet,Pmod,Pdom,Box)
use definitions
	implicit none
	double precision, dimension(3), intent(inout) :: X
	integer, intent(inout)          :: IsInDirichlet
	type(ParametersMod), intent(in) :: Pmod
	type(ParametersDom), intent(in) :: Pdom
	type(ParametersBox), intent(in) :: Box
	integer                         :: ix, iy, iz, k
	
	call findBox(X,k,Pdom)
	if ( (k.lt.1).or.(k.gt.8*(Pdom%Nx+1)*(Pdom%Ny+1)*(Pdom%Nz+1)) ) then
		write(*,'(a,i0,a)') "Error : motion-step bring a fiber into the box ",k,", which does not exist !" 
		STOP
	end if
	
	! If a fiber exit the computation domain by entering a periodic ghost box : make it come back from the opposite side.
	if (Box%BoxType(k,1).eq.1) then
		X(:) = X(:) - Box%BoxType(k,2:4)
		
		call findBox(X,k,Pdom)		
		if (Box%BoxType(k,1).eq.1) then
			write(*,'(a)') "Error : fiber is still in a periodic box after applying boundary condition."
			STOP
		end if
	end if
	
	! If a fiber exit the computation domain by entering a must-be-empty box : bring it back inside from the same side.
	if (Box%BoxType(k,1).eq.2) then
		ix = floor( (X(1)+Pdom%xmax)/Pdom%dx )
		iy = floor( (X(2)+Pdom%ymax)/Pdom%dy )
		iz = floor( (X(3)+Pdom%zmax)/Pdom%dz )
		
		if ( (ix.eq.0).and.(Pmod%DirLim(1).eq.0).and.(Pmod%DirLim(2).eq.1) ) then
			X(1) = - 0.999*Pdom%xmaxR
		else if ( (ix.eq.2*Pdom%Nx+1).and.(Pmod%DirLim(1).eq.1).and.(Pmod%DirLim(2).eq.0) ) then
			X(1) = + 0.999*Pdom%xmaxR
		end if
		
		if ( (iy.eq.0).and.(Pmod%DirLim(3).eq.0).and.(Pmod%DirLim(4).eq.1) ) then
			X(2) = - 0.999*Pdom%ymaxR
		else if ( (iy.eq.2*Pdom%Ny+1).and.(Pmod%DirLim(3).eq.1).and.(Pmod%DirLim(4).eq.0) ) then
			X(2) = + 0.999*Pdom%ymaxR
		end if
		
		if ( (iz.eq.0).and.(Pmod%DirLim(5).eq.0).and.(Pmod%DirLim(6).eq.1) ) then
			X(3) = - 0.999*Pdom%zmaxR
		else if ( (iz.eq.2*Pdom%Nz+1).and.(Pmod%DirLim(5).eq.1).and.(Pmod%DirLim(6).eq.0) ) then
			X(3) = + 0.999*Pdom%zmaxR
		end if
		
		call findBox(X,k,Pdom)
		if (Box%BoxType(k,1).eq.2) then
			write(*,'(a)') "Error : fiber is still in a must-be-empty box after applying boundary condition."
			STOP
		end if
	end if
	
	! If a fiber exit the computation domain by entering a dirichlet box : fix it here.
	if (Box%BoxType(k,1).eq.3) then
		if (IsInDirichlet.eq.0) then
			IsInDirichlet = 1
		else if (IsInDirichlet.eq.1) then
			IsInDirichlet = 2
		end if
	end if
end subroutine



subroutine RandomInitialization(Xcell,Rcell,Xfib,omega,link,distpoints,Pmod,Pdom,BoxCell,BoxFib)
use definitions
	implicit none
	double precision, dimension(:,:), intent(inout) :: Xcell, Xfib, omega, distpoints
	integer, dimension(:,:), intent(inout)          :: link
	double precision, dimension(:), intent(inout)   :: Rcell
	type(ParametersMod), intent(in)                 :: Pmod
	type(ParametersDom), intent(in)                 :: Pdom
	type(ParametersBox), intent(inout)              :: BoxCell,BoxFib

	double precision, parameter :: pi=3.141592653589793238
	double precision :: alea1,alea2,alea3,chil,dik,lik,lki,Xtemp(3),ratio
	integer          :: i, k, l, l1, l2, l3, ltemp
	integer          :: ix, iy, iz, Ntemp, Ndirperbox

	!------ Cells ----------------------------------------------------------------------------------
	ratio = 4.0
	do i=1,Pmod%Ncell
		! Positions homogeneously distributed inside the whole domain
		call random_number(alea1)
		call random_number(alea2)
		call random_number(alea3)
		Xcell(i,1) = (2d0*Pdom%xmaxR*alea1 - Pdom%xmaxR)/ratio
		Xcell(i,2) = (2d0*Pdom%ymaxR*alea2 - Pdom%ymaxR)/ratio
		Xcell(i,3) = (2d0*Pdom%zmaxR*alea3 - Pdom%zmaxR)/ratio
	
		! Security in case of a change in the initialization process
		if (abs(Xcell(i,1))>Pdom%xmaxR) Xcell(i,1) = Xcell(i,1) - sign(1d0,Xfib(i,1))*2d0*Pdom%xmaxR
		if (abs(Xcell(i,2))>Pdom%ymaxR) Xcell(i,2) = Xcell(i,2) - sign(1d0,Xfib(i,2))*2d0*Pdom%ymaxR
		if (abs(Xcell(i,3))>Pdom%zmaxR) Xcell(i,3) = Xcell(i,3) - sign(1d0,Xfib(i,3))*2d0*Pdom%zmaxR

		! Each and every cell begin with the minimal radius Rmin
		Rcell(i) = Pmod%Rmin	
	end do
	
	! Attribute cells to boxes
	BoxCell%FirstAgent(:) = 0
	BoxCell%LastAgent(:)  = 0
	BoxCell%AgentBox(:)   = 0
	BoxCell%VerletList(:) = 0
	
	do i = 1,Pmod%Ncell
		call UpdateVerlet(i,Xcell(i,:),Pdom,BoxCell)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxCell)
	
	
	
	!------ Fibers ---------------------------------------------------------------------------------
	do i=1,Pmod%Nfib
		! Positions homogeneously distributed inside the whole domain
		call random_number(alea1)
		call random_number(alea2)
		call random_number(alea3)
		Xfib(i,1) = 2d0*Pdom%xmaxR*alea1 - Pdom%xmaxR
		Xfib(i,2) = 2d0*Pdom%ymaxR*alea2 - Pdom%ymaxR
		Xfib(i,3) = 2d0*Pdom%zmaxR*alea3 - Pdom%zmaxR
	
		! Security in case of a change in the initialization process
		if (abs(Xfib(i,1))>Pdom%xmaxR) Xfib(i,1) = Xfib(i,1) - sign(1d0,Xfib(i,1))*2d0*Pdom%xmaxR
		if (abs(Xfib(i,2))>Pdom%ymaxR) Xfib(i,2) = Xfib(i,2) - sign(1d0,Xfib(i,2))*2d0*Pdom%ymaxR
		if (abs(Xfib(i,3))>Pdom%zmaxR) Xfib(i,3) = Xfib(i,3) - sign(1d0,Xfib(i,3))*2d0*Pdom%zmaxR

		! Orientations homogeneously distributed in the sphere
		call pick_and_reject(omega(i,:))	
	end do
	
	
	!------ Dirichlet layer fibers -----------------------------------------------------------------
	Ndirperbox = int(Pmod%DirFillingRate * Pdom%dx*Pdom%dy*Pdom%dz)
	Ntemp = Pmod%Nfib
	
	do ix = 1,2*Pdom%Nx+2
		do iy = 1,2*Pdom%Ny+2
			do iz = 1,2*Pdom%Nz+2
				l = ix + 2*(iy-1)*(Pdom%Nx+1) + 4*(iz-1)*(Pdom%Nx+1)*(Pdom%Ny+1)

				if (BoxFib%BoxType(l,1).eq.3) then
					k = 1
					do while(k.le.Ndirperbox)
						!  Positions homogeneously distributed inside box l = (ix,iy,iz)
						call random_number(alea1)
						call random_number(alea2)
						call random_number(alea3)
						Xfib(Ntemp+k,1) = alea1*Pdom%dx + (ix-1)*Pdom%dx - Pdom%xmax
						Xfib(Ntemp+k,2) = alea2*Pdom%dy + (iy-1)*Pdom%dy - Pdom%ymax
						Xfib(Ntemp+k,3) = alea3*Pdom%dz + (iz-1)*Pdom%dz - Pdom%zmax
        				
						! Orientation parallele to the border between box l and the inside of the domain
						call random_number(alea1)
						if ( (ix.eq.1).or.(ix.eq.2*Pdom%Nx+2) ) then
							omega(Ntemp+k,1) = 0.0
							omega(Ntemp+k,2) = cos(2*pi*alea1)
							omega(Ntemp+k,3) = sin(2*pi*alea1)
						elseif ( (iy.eq.1).or.(iy.eq.2*Pdom%Ny+2) ) then
							omega(Ntemp+k,1) = cos(2*pi*alea1)
							omega(Ntemp+k,2) = 0.0
							omega(Ntemp+k,3) = sin(2*pi*alea1)
						elseif ( (iz.eq.1).or.(iz.eq.2*Pdom%Nz+2) ) then
							omega(Ntemp+k,1) = cos(2*pi*alea1)
							omega(Ntemp+k,2) = sin(2*pi*alea1)
							omega(Ntemp+k,3) = 0.0
						end if

						k = k+1
					end do
					Ntemp = Ntemp + Ndirperbox
				end if
				
			end do
		end do
	end do	
	
	
	! Attribute ALL fibers (including Dirichlet ones) to boxes
	BoxFib%FirstAgent(:) = 0
	BoxFib%LastAgent(:)  = 0
	BoxFib%AgentBox(:)   = 0
	BoxFib%VerletList(:) = 0
	
	do i = 1,Pmod%Nfib+Pmod%Ndirtot
		call UpdateVerlet(i,Xfib(i,:),Pdom,BoxFib)
	end do
	call FillGhostBoxes(Pmod,Pdom,BoxFib)
	

	!------ Links ----------------------------------------------------------------------------------
	link = 0
	distpoints = 0
	
	if ((Pmod%Nfib.ne.0).and.(Pmod%LinkType.ne.2)) then
	
		chil = Pmod%freq_link*1.0/(Pmod%freq_link + Pmod%freq_dlink)
		
		do i=1,Pmod%Nfib
			! Pour rechercher les voisins on utilise une méthode de localisation sur la grille numérique (appelée VerletMethod)
			l = BoxFib%AgentBox(i)
			do l1 = -1,1
				do l2 = -1,1
					do l3 = -1,1
						ltemp = l + l1 + 2*l2*(Pdom%Nx+1) + 4*l3*(Pdom%Nx+1)*(Pdom%Ny+1)
						k = BoxFib%FirstAgent(ltemp)
						
						do while (k /= 0)  ! i is real | k is ghost/real/dirichlet
							if ((k<=Pmod%Nfib).and.(k>i)) then ! k is a moving fiber (ghost or real) different from i
								
								if (BoxFib%BoxType(ltemp,1).eq.1) then
									Xtemp(:) = Xfib(k,:) + BoxFib%BoxType(ltemp,2:4)
								else
									Xtemp(:) = Xfib(k,:)
								end if
								
								call find_overlap_fibfib(Xfib(i,:),Xtemp,omega(i,:),omega(k,:),Pmod,dik,lik,lki)
								call random_number(alea1)
								
								if ( ((dik-2*Pmod%Rf-Pmod%eps)<0d0).and.(alea1<chil) ) then
									link(i,k) = 1
									link(k,i) = 1
									distpoints(i,k) = lik
									distpoints(k,i) = lki
								end if
								
							elseif (k>Pmod%Nfib) then ! k is a fiber from the Dirichlet layer
								
								call find_overlap_fibfib(Xfib(i,:),Xfib(k,:),omega(i,:),omega(k,:),Pmod,dik,lik,lki)
								call random_number(alea1)
								
								if ( ((dik-2*Pmod%Rf-Pmod%eps)<0d0).and.(alea1<chil) ) then
									link(i,k) = 1
									distpoints(i,k) = lik
									distpoints(k,i) = lki
								end if
								
							end if
							k = BoxFib%VerletList(k)
						end do
						
					end do
				end do
			end do
		end do
		
	end if
end subroutine RandomInitialization



end module fonctions