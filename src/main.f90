program calprep
   !
   !Prepare met and geo data for CALMET.
   !
   use prep_uppa
   use prep_surf
   use calgeo      !prep geo.dat

   implicit none
   integer, parameter :: max_files=100
   !control:
   character(19)      :: start_date, end_date
   character(256)     :: proj
   integer            :: nx,ny
   real               :: dx,dy,xc,yc
   logical            :: prep_surf=.false. ,prep_up=.false. ,prep_geo=.false.
   !surface:
   character(len=256) :: surface_files(max_files)
   integer            :: surface_nsta
   !upperair:
   character(len=256) :: upperair_files(20)
   real               :: upperair_ptop
   !geo
   character(len=256) :: terrain_file
   character(len=256) :: lulc_file
   character(len=8)   :: lulc_lookup(100)=""
   !general
   integer            :: ios
   character(len=256) :: msg

   !---read namelist variables and parameters
   namelist/control /start_date,end_date,proj,dx,dy,nx,ny,xc,yc,prep_surf,prep_up,prep_geo
   namelist/surface /surface_files,surface_nsta
   namelist/upperair/upperair_files,upperair_ptop
   namelist/geo     /terrain_file,lulc_file,lulc_lookup    
   !namelist/prog    /wrf_file,geo_file!,lulc_files,lulc_categories

   print '("CALPREP",/)'

   print '(" SETUP PHASE")'

   open(1,file='namelist.calprep')
     read(1,nml=control , iostat=ios,iomsg=msg);  if (ios /= 0) then; print*,msg;stop 'Error: Failed to read control  section on namelist.';endif!control
     read(1,nml=surface , iostat=ios,iomsg=msg);  if (ios /= 0) then; print*,msg;stop 'Error: Failed to read surface  section on namelist.';endif!surface
     read(1,nml=upperair, iostat=ios,iomsg=msg);  if (ios /= 0) then; print*,msg;stop 'Error: Failed to read upperair section on namelist.';endif!raobs   
     read(1,nml=geo     , iostat=ios,iomsg=msg);  if (ios /= 0) then; print*,msg;stop 'Error: Failed to read geo      section on namelist.';endif!geo   
   close(1)


   print '(" COMPUTATIONAL PHASE")'

   !--- Surface  Met. data (ISH) -> surf.dat
   if ( prep_surf ) call ish2surf(start_date, end_date, surface_files(1:surface_nsta))
  
   !--- Upperair Met. data (IGRA) -> up.dat
   if ( prep_up   ) call igra2up(start_date, end_date, upperair_ptop, upperair_files(:))

   !--- Topography & surface parameters  (GTIFF) -> geo.dat
   if ( prep_geo  ) call make_geo(trim(terrain_file),trim(lulc_file), xc,yc,dx,dy,nx,ny,proj, lulc_lookup) 
   
   !!--- Prognostic model data (NetCDF: geo_em & wrfout) -> 3D.dat
   !call calwrf(start_date, end_date, wrfout, geo_em,xc,yc,dx,dy,nx,ny,proj)

   print '(" TERMINATION PHASE")'

end program
