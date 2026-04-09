module calgeo
   !
   ! Module to parse TIFF files to geo.dat file for CALMET
   !

   !use datetime_module, only: datetime, timedelta, strptime
   !use utils_module
   use PROJ
   use GTIFF
   use INTERP

   implicit none

   private
   public make_geo

   include 'landuse.tbl'  !table with land use parameters

contains

subroutine make_geo(terrain_file,lulc_file, xc,yc,dx,dy,nx,ny,proj,lulc_lookup) !lu_file, lu_categories,
   implicit none
   character(len=*), intent(in) :: terrain_file, lulc_file
   character(len=8), intent(in) :: lulc_lookup(:)
   character(len=*), intent(in) :: proj
   real            , intent(in) :: xc,yc,dx,dy
   integer         , intent(in) :: nx,ny           

   type(interp_grid)    :: g        !calmet (target) grid

   real,    allocatable :: z(:,:)   !terrain height array
   integer, allocatable :: lu(:,:)  !land use category array 
   real,    allocatable :: f(:,:,:) !array with relative frequency of each land use category [ncat,nx,ny]
   real,    allocatable :: a0(:,:)  !albedo 
   real,    allocatable :: b0(:,:)  !bowen's ratio
   real,    allocatable :: z0(:,:)  !surface roughness length
   real,    allocatable ::  H(:,:)  !Heat flux 
   real,    allocatable :: Ha(:,:)  !Antropogenic heat flux
   real,    allocatable ::LAI(:,:)  !Leaf Area Index
    
   integer :: i,j,k


   print '("Creating geo.dat from terrain and land use files..")'
   !Initialize target interp_grid values
   call init_grid(g, xc,yc, dx,dy,nx,ny,proj)
   
   allocate(  z(g%nx,g%ny))
   allocate(  f(ncat,g%nx,g%ny))
   allocate( lu(g%nx,g%ny))  
   allocate( a0(g%nx,g%ny))  
   allocate( b0(g%nx,g%ny))  
   allocate( z0(g%nx,g%ny))  
   allocate(  H(g%nx,g%ny))  
   allocate( Ha(g%nx,g%ny))  
   allocate(lai(g%nx,g%ny))  

   !Topography
   call get_topo(terrain_file,g,z)

   !Land Use/Land cover cellgrid fractions (f)
   call get_lulc(lulc_file ,g,lu,lulc_lookup, f)

   !Compute surface characteristics
   do concurrent(i = 1:g%nx, j = 1:g%ny)
     lu(i,j) = cat1(maxloc(f(:,i,j),1))  
     a0(i,j) = dot_product(f(:,i,j),  a0_t(:))
     b0(i,j) = dot_product(f(:,i,j),  b0_t(:))
     z0(i,j) = dot_product(f(:,i,j),  z0_t(:))
      H(i,j) = dot_product(f(:,i,j),   H_t(:))
     Ha(i,j) = 0.0 !dot_product(f(:,i,j) * Ha_t(:))
     LAI(i,j)= dot_product(f(:,i,j), LAI_t(:))
   end do

   !Write geo.dat
   call write_geo(g,z,lu,a0,b0,z0,H,Ha,LAI)

end subroutine

subroutine init_grid(g, xc,yc,dx,dy,nx,ny,proj)
   implicit none
   type(interp_grid), intent(inout) :: g
   real             , intent(in)    :: xc,yc,dx,dy
   integer          , intent(in)    :: nx,ny           
   character(len=*) , intent(in)    :: proj

   print '("  Initialize grid parameters")'
   !general
   g%xc = dble(xc); g%yc = dble(yc)
   g%dx = dble(dx); g%dy = dble(dy)
   g%nx = nx      ; g%ny = ny
   g%proj=trim(proj)
   !boundaries
   g%xmin=dble(xc-0.5*dx*nx) 
   g%xmax=dble(xc+0.5*dx*nx)
   g%ymin=dble(yc-0.5*dy*ny)
   g%ymax=dble(yc+0.5*dy*ny)
end subroutine

subroutine get_topo(terrain_file,g,z)
  implicit none
  type(interp_grid) ,intent(in)    :: g
  character(len=*)  ,intent(in)    :: terrain_file
  real, allocatable ,intent(inout) :: z(:,:)   !calmet grid values
  type(interp_grid) :: g0
  real, allocatable,dimension(:,:) :: x0,y0,z0 !original TIFF values
  type(TIFF_FILE)   :: my_tiff

  integer :: ierr,crs!,wid,len

   print '("  Processing topography data..")'
  !Read GTIFF DEM file:
  call TIFF_Open(124,trim(terrain_file),'r', my_tiff, ierr)
  if (ierr==0) then
  
     allocate(x0(my_tiff%nx,my_tiff%ny))
     allocate(y0(my_tiff%nx,my_tiff%ny))
     allocate(z0(my_tiff%nx,my_tiff%ny))

     call GTIFF_GET_IMAGE_COORDINATES(my_tiff,x0,y0)
     call TIFF_GET_IMAGE(my_tiff, 1, z0)
     call GTIFF_GET_PROJ_STR(my_tiff, 1, g0%proj)

     !Pass values to interp_grid
     g0%nx  =my_tiff%nx      ; g0%ny=my_tiff%ny
     g0%dx  =my_tiff%scale(1); g0%dy=my_tiff%scale(2)
     g0%xmin=minval(x0)      ; g0%ymin=minval(y0)    
     g0%xmax=maxval(x0)      ; g0%ymax=maxval(y0)    
     g0%xc  =g0%xmin+0.5*g0%dx*g0%nx
     g0%yc  =g0%ymin+0.5*g0%dy*g0%ny

     call TIFF_Close(my_tiff)
  else
     print '("Failed to read TIFF file: ",A)',trim(terrain_file); stop
  endif

  !Interpolate read arrays to target grid:
  call interp_2d(g0,z0,g,z)

end subroutine

subroutine get_lulc(lulc_file, g, lu, lu_lookup, freq)
  implicit none
  character(len=*)  ,intent(in)       :: lulc_file 
  type(interp_grid) ,intent(in)       :: g                                                   
  integer, allocatable, intent(inout) :: lu(:,:)       !calmet grid values                         
  character(len=8), intent(in)        :: lu_lookup(:)
  real,    allocatable, intent(inout) :: freq(:,:,:)

  type(interp_grid)     :: g0                                                  
  real, allocatable     :: x0(:,:),y0(:,:),lu0(:,:)    !original TIFF values      
  integer, allocatable  ::ilu0(:,:)                    !original TIFF values (int)
  type(TIFF_FILE) :: my_tiff                                                              
  integer :: ierr,crs                                                                        

  !lulc.dat file
  integer   :: io8=18
  character :: ccomma=','
  integer :: i,j,k,n
      
  print '("  Processing land-use data..")'
  !Read GTIFF LULC file:
  call TIFF_Open(124,trim(lulc_file),'r', my_tiff, ierr)                                   
  if (ierr==0) then                                                                       
                                                                                          
     allocate( x0(my_tiff%nx,my_tiff%ny))                                                  
     allocate( y0(my_tiff%nx,my_tiff%ny))                                                  
     allocate(lu0(my_tiff%nx,my_tiff%ny))                                                  

     call GTIFF_GET_PROJ_STR(my_tiff, 1, g0%proj)
     call GTIFF_GET_IMAGE_COORDINATES(my_tiff,x0,y0)                                      
     call TIFF_GET_IMAGE(my_tiff, 1, lu0)
     ilu0=int(lu0)  !transfer real to integer.

     !Pass values to interp_grid                                                          
     g0%nx  =my_tiff%nx      ; g0%ny=my_tiff%ny                                       
     g0%dx  =my_tiff%scale(1); g0%dy=my_tiff%scale(2)
     g0%xmin=minval(x0)      ; g0%ymin=minval(y0)    
     g0%xmax=maxval(x0)      ; g0%ymax=maxval(y0)    
     g0%xc  =g0%xmin+0.5*g0%dx*g0%nx
     g0%yc  =g0%ymin+0.5*g0%dy*g0%ny

     call TIFF_Close(my_tiff)
  else
     print '("Failed to read TIFF file: ",A)',trim(lulc_file); stop
  endif

  print '("  Mapping land-use categories..")'
  !Map LU categories to CALPUFF LU
  call map_categories(ilu0,lu_lookup)

  !!Interpolate read arrays to target grid:
  !call interp_2d(g0,ilu0,g,lu)
  
  print '("  Computing land-use relative fraction in grid cells..")'
  call remap_freq(g0,ilu0,g,freq,cat2)

  print '("  Writing lulc.dat file..")'
  open(io8,file='lulc.dat', action='write',status='replace')
     write(io8,'(10x,999(i6))') (cat2(n),n=1,ncat)
     do j=1,g%ny                                 
     do i=1,g%nx                                                   
        write(io8,'(2i5,40(f6.3))') i,j,(freq(k,i,j),k=1,ncat)
     enddo
     enddo
  close(io8)

  deallocate(x0,y0,lu0)
end subroutine

subroutine map_categories(iarr,lookup)
     implicit none
     integer, intent(inout)       :: iarr(:,:)
     character(len=8),intent(in) :: lookup(:)
     integer :: key,val, i,j,k,pos
     character(len=8) :: left,right,pair

     print '("   Mapping land use categories..")'

     ! Loop over lookup table and apply replacements with WHERE
     do k = 1, size(lookup)
       pair= trim(lookup(k))
       pos = index(pair, ":")
       if (pos > 0) then
         left =trim( pair(:pos-1) )
         right=trim( pair(pos+1:) )

         read(left, *)  key
         read(right, *) val

         if ( key /= val) then
                 print '("     - mapping key",i5,"->","val",i5)',key, val
            where (iarr == key)
              iarr = val
            end where
         end if
       else
          !print*, "ERROR on lookup array:",lookup(k)
       end if
     end do
end subroutine


subroutine write_geo(g,z,lu,a0,b0,z0,H,Ha,LAI)
   implicit none
   type(interp_grid),      intent(in) :: g
   real,                   intent(in) :: z(:,:)
   integer,                intent(in) :: lu(:,:)
   real   ,dimension(:,:), intent(in) :: a0,b0,z0,H,Ha,LAI

   integer :: iost, io7=7
   integer :: i,j,k,n
   character(len=1) :: ccomma=","
   
   print '("   Writing geo.dat file..")'

   open(io7,file='geo.dat', action='write',status='replace')
   
      !!HEADER:
      !write(io7,*)"#GEO.DAT     3.0     Produced by CALPREP version 0.0"
      !write(io7,*)"#"//adjustl(g%proj)
      !write(io7,*)"ncols",     g%nx
      !write(io7,*)"nrows",     g%ny
      !write(io7,*)"xllcorner", g%xmin
      !write(io7,*)"yllcorner", g%ymin   
      !write(io7,*)"cellsize",  g%dx
      !write(io7,*)"NODATA_value", -9999  
      !HEADER: (original format)
      write(io7,'("GEO.DAT",/,"1",/,"Produced by makegeo v0.0")') !title
      write(io7,'("UTM")')
      write(io7,'(a10)') trim(g%proj)                               !proj 
      write(io7,'("WGS-84   01-01-2000")')
      write(io7,'(2i8,4f12.3)') g%nx,g%ny,g%xmin*1e-03,g%ymin*1e-03,g%dx*1e-03,g%dy*1e-03!grid 
      write(io7,'(a4)') 'KM'
      !LANDUSE:
      !write(io7,'(1x,f6.4,1x," - ",a70)') 1.0,'LAND USE DATA - (1=new categories)'
      write(io7,'(1x,i6,1x," - ",a70)') 1  ,'LAND USE DATA - (1=new categories)'
      write(io7,'("14  51  55 - NLU, IWAT1, IWAT2")')
      write(io7,'("10  20 -20  30  40  51  54  55  60  61  62  70  80  90")')
      !write(io7,'("#",30a)') 'land use dat' 
      do j=g%ny,1,-1
         write(io7,'(10(i7,a1))') (lu(n,j),ccomma,n=1,g%nx-1),lu(g%nx,j)
         !write(io7,'(100(i7,a1))') (lu(n,j),ccomma,n=1,g%nx-1),lu(g%nx,j)
      enddo
      !TERRAIN:
      write(io7,'(1x,f6.4,1x," - ",a70)') 1.0,'TERRAIN heights - HTFAC (Conversion to meters)'
      do j=g%ny,1,-1
         write(io7,'(100(f7.2,a1))') (z(n,j),ccomma,n=1,g%nx-1),z(g%nx,j)
      enddo
      !Z0
      write(io7,'(1x,"2",3x," - ",a70)') 'gridded z0 field'
      do j=g%ny,1,-1
         write(io7,'(100(f7.2,a1))') (z0(n,j),ccomma,n=1,g%nx-1),z0(g%nx,j)
      enddo
      !A0
      write(io7,'(1x,"2",3x," - ",a70)') 'gridded albedo field'
      do j=g%ny,1,-1
         write(io7,'(100(f7.2,a1))') (a0(n,j),ccomma,n=1,g%nx-1),a0(g%nx,j)
      enddo
      !B0
      write(io7,'(1x,"2",3x," - ",a70)') 'gridded Bowen ratio field'
      do j=g%ny,1,-1
         write(io7,'(100(f7.2,a1))') (b0(n,j),ccomma,n=1,g%nx-1),b0(g%nx,j)
      enddo
      !H
      write(io7,'(1x,"2",3x," - ",a70)') 'gridded soil heat flux parameters'
      do j=g%ny,1,-1
         write(io7,'(100(f7.2,a1))') (H(n,j),ccomma,n=1,g%nx-1),H(g%nx,j)
      enddo
      !H_antro
      write(io7,'(1x,"2",3x," - ",a70)') 'gridded anthropogenic heat flux field'
      do j=g%ny,1,-1
         write(io7,'(100(f7.2,a1))') (Ha(n,j),ccomma,n=1,g%nx-1),Ha(g%nx,j)
      enddo
      !LAI
      write(io7,'(1x,"2",3x," - ",a70)') 'gridded leaf area index field'
      do j=g%ny,1,-1
         write(io7,'(100(f7.2,a1))') (LAI(n,j),ccomma,n=1,g%nx-1),LAI(g%nx,j)
      enddo

   close(io7)

end subroutine



end module
