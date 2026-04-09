module prep_uppa
   !
   ! Module to parse IGRA files to UP.dat files for CALMET
   !
   use datetime_module, only: datetime, timedelta, strptime! clock,
   use utils_module

   implicit none
   integer, parameter :: mxlev=200          !max number of levels a sounding could have

   type igra_record
     real    :: p,z,t,ws!,wd,rh,dwp         !pressure, height, temperature, wind-speed
     integer :: wd                          !wind-direction
   end type

   type igra_sounding
      character(11)      :: id              !station id
      real               :: lon,lat,elev    !longitude, latitude, altitude
      integer            :: nlev            !sounding number of levels
      integer            :: istop           !level at which pressure >= pstop
      type(datetime)     :: date
      type(igra_record)  :: r(mxlev)        !array with records
      logical :: first_sounding=.true.      !flag for the first sounding to be written.
   end type

   private 
   public igra2up
   
contains

subroutine igra2up(sdate, edate, pstop, file_list) !iFile)
    implicit none
    character(*)  ,intent(in)  :: file_list(:)          !Input Surfac Met. (ISH) files list.
    !character(*)  ,intent(in)  :: iFile
    character(256)             :: iFile,oFile
    character(19) ,intent(in)  :: sdate, edate
    real          ,intent(in)  :: pStop
    type(datetime)             :: current_date,next_date,start_date, end_date
    logical                    :: file_exists=.false.   !flag file existence
    type(igra_sounding)        :: sounding              !sounding buffer
    integer                    :: io8=8, io9=9          !file unit numbers for IO 
    integer :: nn

    print '("Parsing IGRA file to up.dat ..    "/)'

    do nn=1,size(file_list)

       iFile=trim(file_list(nn))                                           ! set igra     file name
       write(oFile,'("up",i0.2,".dat")') nn                                ! set upnn.dat file name
       sounding = new_sounding()                                           ! create/reset igra sounding

       inquire(file=trim(iFile), exist=file_exists)                        !Check if input file exists.
       if ( file_exists ) then

          print '("Reading file: ",(A))', trim(iFile)
          open (io8, file=trim(iFile), action="read" , status='old'    )   !Open Input  file (IGRA file)
          open (io9, file=trim(oFile), action='write', status='replace')   !Open Output file (UPnn.DAT)

          start_date   = strptime(sdate,'%Y-%m-%d %H:%M:%S') 
          end_date     = strptime(edate,'%Y-%m-%d %H:%M:%S')
          current_date = datetime()
          next_date    = start_date
          do while (current_date <= end_date)                              !Loop over each record till end_date is reached.
             
             call Read_IGRA_sounding_header(io8, sounding)                 !Read IGRA next header, put content into sounding

             current_date = sounding%date                                  !Update current date-time
           
             if ( current_date >= start_date ) then                        !Proceed if start_date has been reached.

                call Read_IGRA_sounding(io8, sounding, pstop)              !Read IGRA sounding content, put content into sounding

                if ( sounding%first_sounding ) then                        !If this is the first sounding to be written
                   call Write_header(io9, sounding, start_date, end_date, pstop)  !Write output header
                   sounding%first_sounding=.false.                         !Future soundings will no longer be the first
                end if
                
                if ( current_date > next_date) then                        !Check if previous date is what expected (or previous)
                   print '("Warning: No soundings for this date: ",A,".")', next_date % isoformat();
                endif

                call Write_sounding(io9, sounding)                         !Write output sounding record

                next_date = current_date + timedelta(days=1)               !Define next expected date

             else                                                          !If sounding date before to start_date 
                call Skip_N_lines(io8, sounding%nlev)                      !Go to next sounding
             end if
          end do
          close(io8)                                                       !Close Input  file
          close(io9)                                                       !Close Output file

          !Execution error checks:
          if ( sounding%first_sounding ) then
             stop 'ERROR. No soundings for the dates requested found in provided file.'
          end if
       else                                                                !If file doesn't exist print error.
          print '("ERROR. File not found: ",A)',trim(iFile); !stop
       end if 

   enddo !file_list

   print '("up.dat created succesfully          ")'
   print '("-----------------------------------"/)'

end subroutine

pure elemental type(igra_record) function new_record()
   new_record % p  = 0.0
   new_record % z  = 0.0
   new_record % t  = 0.0
   new_record % ws = 0.0
   new_record % wd = 0
end function new_record
 
pure elemental type(igra_sounding) function new_sounding()
    new_sounding % id    = ''
    new_sounding % lon   = 0.
    new_sounding % lat   = 0.
    new_sounding % elev  = 0.
    new_sounding % nlev  = 0.
    new_sounding % istop = 0.
    new_sounding % date  = datetime()
    new_sounding % first_sounding = .true.
    new_sounding % r(:)  = new_record ()
end function new_sounding

subroutine read_igra_sounding_header(io, s)
     implicit none
     integer            ,intent(in)    :: io
     type(igra_sounding),intent(inout) :: s
     character(72) :: row
     !print*,"> read_igra_sounding_header."

     read(io,'(A)')  row
     
     if ( row(1:1) .ne. '#' ) then           !header line should start with '#'
        print '("Warning: header not found where expected, searching for the next header.",/,A)',row
        do while ( row(1:1) .ne. "#" )
           read(io,'(A)')  row
        end do
     endif

     ! --- Read IGRA header                  !var      position  type  description
     !                                       !----------------------------------------
     !hr     =row(1:1)                       !HEADREC   1-  1    char  '#'       
     s%id = row(2:12)                        !ID        2- 12    char  station id code
     !yyyy = row(14:17)                      !YEAR     14- 17    int   year  of release
     !mm   = row(19:20)                      !MONTH    19- 20    int   month of release
     !dd   = row(22:23)                      !DAY      22- 23    int   day   of release
     !hh   = row(25:26)                      !HOUR     25- 26    int   nominal obs hour       
     !@hhmm   = row(28:31)                   !RELTIME  28- 31    int   release time (%h%M)
     s%nlev  = atoi(row(33:36))              !NUMLEV   33- 36    int   num of levels of the sounding
     !@ p_src  = row(38:45)                  !P_SRC    38- 45    char  src of pressure levels
     !@ np_src = row(47:54)                  !NP_SRC   47- 54    char  src of non-press lvls
     s%lat = real(atoi(row(56:62)) )/1e4     !LAT      56- 62    int   latitude
     s%lon = real(atoi(row(64:71)) )/1e4     !LON      64- 71    int   longitude
     !---------------------------------------!----------------------------------------
                                             !Example: #USM00072249 2021 01 01 00 2341  195 ncdc-nws ncdc-nws  328350  -972986
     s%date=strptime(row(14:26),"%Y %m %d %H") !sounding date

end subroutine

subroutine skip_n_lines(io, n)
     implicit none
     integer            , intent(in)  :: io,n
     character(72)      :: row
     integer            :: i
     !print*,"> go_to_next_sounding."
     do i=1,n
        read(io,'(A)') row
     enddo
end subroutine

subroutine read_igra_sounding(io, s, pstop)
     implicit none
     integer             ,intent(in)    :: io
     type(igra_sounding) ,intent(inout) :: s
     real               , intent(in)    :: pstop
     character(53) :: row
     character(4)  :: elev
     integer :: i,ielev
     !print*,"> read_igra_sounding."

     s%istop=1                                      !initialize number of valid records
     do i=1, s%nlev                                 !loop over record levels

        read(io,'(A)')  row 
        ! --- Read IGRA sounding
        !Variable Position Type Description
        !----------------------------------------
        !LVLTYP1   1- 1    int  lvl indicator: (1:std pres, 2:Other pres., 3:non-pres.)                  
        !LVLTYP2   2- 2    int  lvl indicator: (1:surf.   , 2:tropop.    , 0:other) 
        !ETIME     4- 8    int  elapsed time since lunch                  
        !PRESS    10-15    int  pressure [Pa or mb*100] missing:-9999 or -8888         
        !PFLAG    16-16    char QA flag:(blank: not checked, A: tier-1, B: tier-1 & tier-2.)
        !GPH      17-21    int  geopotential height [masl], missing:-9999 or -8888
        !ZFLAG    22-22    char                   
        !TEMP     23-27    int  temperature [ºC*10], missing -9999
        !TFLAG    28-28    char                   
        !RH       29-33    int  relative humidity [%*10]                  
        !DPDP     35-39    int  dew-point [ºC*10]                  
        !WDIR     41-45    int  wind direction [deg from north (clockwise)] 
        !WSPD     47-51    int  wind speed [m/s*10]                  
        !----------------------------------------

        if ( row(1:1) /= "3" ) then                                    !do not use non-pressure levels

          !--- Check first level is at ground level
          if ( row(2:2) == "1" ) then
             !elev=trim(row(18:21))
             !ielev=atoi(elev)
             !print*,"ELEV:",ielev
             s%elev=real(atoi(row(18:21)))                   !get altitude.
             if ( s%istop > 1) print '("Warning: First level not at ground level! (date:",A,")")',s%date % isoformat();
          end if

          !--- Check Pressure value. if valid store it.
          if ( row(10:10) .ne. "-" ) then                              !check pressure isn't missing
             s%r(s%istop)%p=real(atoi(row(10:15)))/100.0;              !get pressure and convert to [mbar]
          else 
             !print '("Warning: Missing pressure at : ",I3,". on date: ",12A)',i,s%date_str;
             cycle
          endif      
          if ( s%r(s%istop)%p >= pstop ) then                          !check pressure > top pressure level (ptop)

             !--- Check Geopotential
             if ( row(17:17) .ne. "-" ) then 
                s%r(s%istop)%z=real(atoi(row(17:21)))                   !get geopotential [m]
             else 
                !print '("Warning: Missing geopotential at: ",I3,". on date: ",12A)',i,s%date_str;
                cycle
             endif

             !--- Check Temperature
             if ( (row(23:27) .ne. "-9999") .or. (row(23:27) .ne. "-8888") ) then 
                s%r(s%istop)%t=real(atoi(row(23:27)))/10.0 + 273.15     !get temp and convert to [deg Kelvin]
             else 
                !print '("Warning: Missing temperature at: ",I3,". on date: ",12A)',i,s%date_str;
                cycle
             endif

             s%r(s%istop)%wd =      atoi(row(41:45))                    !get wind dir
             s%r(s%istop)%ws = real(atoi(row(47:51)))/10.               !get wind speed and convert to [m/s]

             s%istop=s%istop + 1                                        !since is a valid record, increment iStop

          else                                                          !if press < ptop
             call skip_n_lines(io,s%nlev-i); return                     !drop the rest of the sounding
          end if
        else                                                            !if non-pressure level
            continue                                                    !go to next iteration
        end if
     enddo
end subroutine

subroutine write_header(io, s, sdate, edate, pStop)
    implicit none
    integer            , intent(in) :: io
    type(igra_sounding), intent(in) :: s
    type(datetime)     , intent(in) :: sdate, edate
    real               , intent(in) :: PSTOP
    logical      :: LHT=.true.,LTEMP=.true.,LWD=.false.,LWS=.false.

    character(16) ::  dataset,dataver
    character(64) ::  datamod
    character(16) ::  clon,clat      
    character(1)  ::  hem,mer

    ! --- Configure output variables
    data dataset/'UP.DAT'/, dataver/'2.1'/
    data datamod/'Hour Start and End Times with Seconds'/

    integer :: ibyr,ibjd,ibhr, ibsec
    integer :: ieyr,iejd,iehr, iesec
    character(len=18) :: a,b

    a=sdate%strftime("%Y %j %H %S")
    b=edate%strftime("%Y %j %H %S")
    read(a,"(i4,i4,i3,i3)") ibyr,ibjd,ibhr,ibsec
    read(b,"(i4,i4,i3,i3)") ieyr,iejd,iehr,iesec

    !Header:
    write(io,'(2a16,a64)') dataset,dataver,datamod
    write(io,*) '1'                 !n-comments   
    write(io,*) 'EPSG:4326'               !proj
    write(io,'("LL",/,"WGS-84  10-10-2002"/,"DEG")')
    !Start / End time
    write(io,'(a3,sp,i3.2,ss,i2.2)')"UTC",0,0            !time zone (allways UTC)
    write(io,'(1x,8i5,f5.0,2i5)') ibyr,ibjd,ibhr,ibsec,ieyr,iejd,iehr,iesec, pStop,3,2 
    !FORMAT(   1X,8I5,F5.0,2I5)
    !Field flags:
    write(io,'(1x,4(4x,l1))') lht,ltemp,lwd,lws
    !Station coordinates
    hem="N"; if ( s%lat < 0 ) hem="S"
    mer="E"; if ( s%lon < 0 ) mer="W"
    write(clat,"(f11.4,a1)") abs(s%lat),hem
    write(clon,"(f11.4,a1)") abs(s%lon),mer
    write(io,'(a8,3x,a4,3x,2(1x,a12,3x),i7)') s%id(6:11),s%id,clat,clon,int(S%elev)
    ! FORMAT  (A8,2X,A1,A4,A1,2X,2(A1,A16,A1,2X),I7) stnid,q,cname,q,q,clat,q,q,clon,q,ielevm
end subroutine

subroutine write_sounding(io, s)
   implicit none
   integer            ,intent(in)    :: io
   type(igra_sounding),intent(in)    :: s
   integer :: i
   type(datetime) ::this_date,next_date
   integer :: iyr,imo,idy,ihr,isec
   integer :: iyr1,imo1,idy1,ihr1,isec1
   !print*, "> write_sounding"
   character(len=18) :: a,b

   this_date=s%date 
   next_date=this_date + timedelta(hours=6)
   a=this_date%strftime("%Y %m %d %H %S")
   b=next_date%strftime("%Y %m %d %H %S")
   read(a,"(i4,i3,i3,i3,i3)"),iyr ,imo ,idy ,ihr ,isec
   read(b,"(i4,i3,i3,i3,i3)"),iyr1,imo1,idy1,ihr1,isec1
   !c-- CHECKs and repair to implement before writing sounding:
   !=== Non-Reparable problems:      
   !--[x] First level should be at the ground: 
   !--[ ] check pressure monotonical decreasing with levels
   !--[ ] check height   monotonical increaesing with levels
   !--[ ] check variables valid ranges  
   !=== Reparable problems:      
   !--[x] check short sounding                -> discart sounding if valid levels < 4.       
   !--[ ] check missing data at top  
   !--[ ] check missing data at bottom  
   !--[x] check missing value indicator       -> already implemented while reading
   !--[ ] check variables valid ranges  
   !      -   0  < WD < 360  deg
   !      -   0  < WS < ?    m/s
   !      - 175  < T  < 322  K
   !      -   0  < P  < 1040 mbar
   if ( s%istop > 4 ) then
      ! --- Explicit beg/ending times with seconds (format 2.1) (FRR 041123)
      write(io,'(3x,"6201",2x,A8,2(4x,i4,i4,i3,i3,i5),i5,3x,i5)'),s%id(6:11),iyr,imo,idy,ihr,isec,iyr1,imo1,idy1,ihr1,isec1,s%nlev,s%istop-1
      !FORMAT(   3X,'6201',2X,A8,2(4x,i4,i4,i3,i3,i5),i5,3x,i5)
      ! --- Write a comma-delimited file
      write(io,"(4(3x,f6.1,',',f5.0,',',f5.1,',',i3,',',f5.1,','))") (s%r(i)%p,s%r(i)%z,s%r(i)%t,s%r(i)%wd,s%r(i)%ws, I=1,s%istop-1)
      !FORMAT   (4(3X,F6.1,',',F5.0,',',F5.1,',',I3,',',f5.1,a1))

   else
      print '("Warning: Sounding descarted, too few valid levels (date-time: ",A,")")',s%date % isoformat()
   end if
end subroutine

end module
