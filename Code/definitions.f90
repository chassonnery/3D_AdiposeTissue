module definitions

type ParametersSim
	integer            :: InitType
	double precision   :: TF, dt, Tenr
	character(len=200) :: nomdossier
end type ParametersSim

type ParametersMod
	integer           :: Nfib, Ncell, Nmax, Nseed_ens, DirLim(6), Ndirtot, cumul_replink, LinkType
	double precision  :: Rmin, Rmax, Lf, Rf, eps, d0 
	double precision  :: alpha_repCC, alpha_repCF, alpha_repFF, alpha_rappel, alpha_align
	double precision  :: freq_ens, bias_ens, Kgrowth, freq_link, freq_dlink
	double precision  :: DirFillingRate
end type ParametersMod

type ParametersDom
	integer           :: Nx,Ny,Nz ! Number of cuboid boxes in each direction for agent location (using Cell Linked List method)
	double precision  :: dx,dy,dz ! Box side-lengths
	double precision  :: xmax,ymax,zmax,xmaxR,ymaxR,zmaxR ! Half-lengths of the domain in each direction (with and without an outside layer of ghost boxes)
	integer           :: NxPIC,NyPIC,NzPIC ! Number of cuboid boxes in each direction for agent density (using Particule In Cell method)
	double precision  :: dxPIC,dyPIC,dzPIC ! Box side-lengths
end type ParametersDom


type ParametersBox
	integer, dimension(:), allocatable :: AgentBox
	integer, dimension(:), allocatable :: VerletList
	integer, dimension(:), allocatable :: FirstAgent
	integer, dimension(:), allocatable :: LastAgent
	double precision, dimension(:,:), allocatable :: BoxType
end type ParametersBox


type SaveState
	double precision, dimension(:,:,:), allocatable :: CellState, FibState, LinkState
	double precision, dimension(:,:), allocatable   :: VariousReal ! time, tcpu
	integer, dimension(:,:), allocatable            :: VariousInt ! Nlinks_int, Nlinks_ext, Ncell, kbetween
end type SaveState

contains


end module definitions