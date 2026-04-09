module prep_surf
   !
   ! Module to parse ISH files to SURF.dat & PRECIP.dat files for CALMET
   !
   use datetime_module, only: datetime, timedelta, strptime
   use utils_module

   implicit none
   private 
   public ish2surf

   type observation
      type(datetime) :: date             !datetime of observation
      real    :: ws=0,wd=0,tmp=0,pres=0,prate=0,vis=0 !windspeed [m/s], winddir [deg], temp[k], pressure [mb], precipitation [mm]
      integer :: iceil=0,icc=0,irh=0,ipcode=0    !ceiling-height [m], opaque sky cov. [tenths], rel. humidity [%], precip code
      logical :: valid=.false.
   end type

   type station
      character(len=5)               :: name                 !station id
      integer                        :: id                   !station id
      integer                        :: time_zone=0
      real                           :: lat, lon, alt        !lat,lon,altitude
      real(8)                        :: x,y,z                !lat,lon,altitude
      type(observation), allocatable :: O(:)                 !observation array
      integer :: valid_records                               !number of valid records
      logical :: valid=.false.                               !is a valid station?
   end type
   
contains

subroutine ish2surf(sdate, edate, file_list)!, nsta)
    !
    !PURPOSE: Read ISH file and produce surf.dat and precip.dat files
    !
    implicit none
    character(*)  ,intent(in)   :: file_list(:)          !Input Surfac Met. (ISH) files list.

    character(19) ,intent(in)   :: sdate, edate
    integer                     :: NSTA                   !Number of stations
    type(station) ,allocatable  :: S(:)                   !record buffer
    logical                     :: file_exists=.false.    !flag file existence
    character(200)              :: iFile     
    integer                     :: io8=8,io               !file unit numbers for IO 
    integer                     :: ios                    !I/O status indice
    character(500)              :: row

    type(datetime)   :: current_date,next_date,start_date,end_date
    type(timedelta)  :: dt                                        !delta time
    integer          :: N_Hours,hours                             !Run length [hours]
    integer          :: i,j,k,t

    type(observation) :: O

    print '("Parsing ISH file to surf.dat..    "/)'

    start_date = strptime(sdate,'%Y-%m-%d %H:%M:%S') 
    end_date   = strptime(edate,'%Y-%m-%d %H:%M:%S')

    NSTA = size(file_list)                                          !Get Number of Stations to use

    !Calculate Run length [hours]
    dt         = (end_date - start_date)
    N_Hours    = int(dt%total_seconds()/(60*60)) + 1  

    !debug:
    Print*,"  start: ",  start_date%strftime("day: %Y %j. hour: %H %S")     
    Print*,"    end: ",    end_date%strftime("day: %Y %j. hour: %H %S")     
    Print*,"    #NT: ",   N_hours, "(hours)."

    allocate(S(NSTA))                                              !Allocate station to have NSTA
    
    do i=1,NSTA                                                    !For each station:
       iFile=file_list(i)

       if ( .not. allocated(S(i)%O)) allocate(S(i)%O(N_Hours))

       inquire(file=trim(iFile), exist=file_exists)                 !Check if input file exists.
       if ( file_exists ) then

          print '("   Reading file: ",(A))', trim(iFile)
          open (io8, file=trim(iFile), action="read", status='old') !Open Input  file (IGRA file)

          read(io8,'(A)',iostat=ios)row

          !Read Station Header:                                     !Read ISH first row to get staton data.
          write(S(i)%name,'("SS",i2.2)'),i                          !station-id
          S(i)%id =atoi(row(11:15))                                 !station-id
          S(i)%lat=real(atoi(row(29:34))/1000.)                     !latitude
          S(i)%lon=real(atoi(row(35:41))/1000.)                     !longitude
          S(i)%alt=real(atoi(row(47:51))/1.)                        !base altitude

          print '("   Station ID: ",a7," (Lat,Lon,Alt:",f7.3,2x,f8.4,2x,f5.1,")"/)', S(i)%id, S(i)%lat, S(i)%lon, S(i)%alt

          current_date = strptime(row(16:27),"%Y%m%d%H%M")                  !date  YYYYMMDDHHmm
          next_date    = start_date
          do while ( current_date <= end_date)                              !Loop over each record till end_date is reached.

             current_date = strptime(row(16:27),"%Y%m%d%H%M") 
             !print*,"current: ",current_date%strftime("day: %Y %j. hour: %H %S")     

             if ( current_date >= start_date ) then                        !Proceed if start_date has been reached.

                call read_ISH_record(row, O)

                dt=current_date - start_date
                t=int(dt%total_seconds()/3600) + 1
                !print*,"t=",t, O%date%strftime("%Y %m %d %H")

                if ( t > N_Hours) exit 
                if ( O%valid ) then
                   S(i)%O(t)=O
                   S(i)%valid_records=S(i)%valid_records + 1

                   if ( current_date > next_date ) then                      !Check if previous date is what expected (or previous)
                      print '("Warning: No records for this hour: ",A,".")', next_date % isoformat();
                       dt=current_date - next_date                              
                       hours=int(dt%total_seconds()/3600)                  !calculo cuantas horas me salie
                       S(i)%O(t-hours:t)=O                                 !Relleno repitiendo current observ.
                   endif

                end if

                next_date = current_date + timedelta(hours=1)            !Define next expected date
                !print*,"next:   ",next_date%strftime("day: %Y %j. hour: %H %S")
             end if

             read(io8,'(A)',iostat=ios)row
             if (ios /= 0) exit

          end do
          close(io8)                                                    !Close Input  file

          if ( S(i)%valid_records .le. 0.3*N_Hours ) then
             print '("ERROR. No enough valid records for the dates requested found in provided file.")'
             S(i)%valid=.false.
             NSTA=NSTA-1
          end if
       else                                                             !If file doesn't exist print error.
          print '("ERROR. File not found: ",A)',trim(iFile); 
          S(i)%valid=.false.
          NSTA=NSTA-1
       end if 

    enddo ! each station "i"

    !Write output:
    call Write_Surf_dat  ('surf.dat'  , S, start_date, end_date)

    print '(/"surf.dat created succesfully?      ")'
    print '("-----------------------------------"/)'

end subroutine


subroutine read_ish_record(row, o)
     implicit none
     !integer          , intent(in)    :: io
     character(len=*), intent(in) :: row
     type(observation), intent(inout) :: o
     character(105)               :: fix 
     character(400)               :: extra 
     integer :: j
     ! --- Read ISHD record
     !Variable Position Type   Description
     !----------------------------------------
     !RECSIZE   1- 4    int    total record length - 105 = variable record characters
     !USAF_ID   5-10    int    
     !NCEI_ID  11-15    int    
     !DATE     16-23    int    YYYYMMDD
     !TIME     24-27    int    HHMM
     !LAT      29-34    int    Latitude  coordinate [degrees]*1000
     !LON      35-41    int    Longitude coordinate [degrees]*1000
     !REP_TYP  42-46    ascii  Report type code     (!) ONLY ALOW FM-12
     !ELEV     47-51    int    Altitude             [m]*1
     !----------------MANDATORY  DATA -------------------------------------------------
     !WDIR     61-63    int    Wind drection         [deg] * 1
     !WSPD     66-69    int    Wind speed            [m/s] * 10          
     !ICEIL    71-75    int    Ceilingh Height       [m]   * 1
     !VIS      79-84    int    Visibility            [m]   * 1
     !TEMP     88-92    int    Temperature           [ºC]  * 10
     !DP       94-98    int    Dew-point             [ºC]  * 10
     !SLP      22-22    int    Sea level pressure    [mb]  * 10
     !----------------PRECIPIT.  DATA -------------------------------------------------
     !          1- 3    /AA1/
     !PPT_time  4- 5    int    precipitation record duration [hours]
     !PPT_mm    6- 9    int    precipitation quantity [mm]               
     !----------------CLOUD COV. DATA -------------------------------------------------
     !          1- 3    /GF1/
     !ICC       4- 5    int    cloud coverage        [tenths]
     !----------------PRESSURE   DATA -------------------------------------------------
     !          1- 3    /MA1/
     !ALT       4- 8    int    altimeter setting     [hPa*10]
     !PRES     10-14    int    station pressure      [hPa*10]
     !---------------------------------------------------------------------------------
     !Example: 0154720379638822024010100244+36855-084856FM-16+029499999V0203101N005710051819N016093199+00401+00201999999ADD...
     integer :: n_valid  !number of valid values on record

     integer :: nchar
     character(12) :: date !yyyymmddhhmin
     real   :: wd,ws,vis,tmp,dp,slp
     character :: awdqc,awsqc,aceilqc,avisqc,atempqc,adpqc,aslpqc,wdtype !quality checks
     integer:: iceil,irh

     real :: tmpf, dpf    !vars used to compute rh
     real :: pp_mm,pp_t   !vars used to compute prate & pcode

     !Fixed position values:
     fix=row(1:105)
     !read(fix,"(i4,11x,a12,33x,f3.0,2x,f4.1,1x,i5,3x,f6.3,3x,f5.1,1x,f5.1,1x,f5.1,1x)"),nchar,date,wd,ws,iceil,vis,tmp,dp,slp
     !debug: 
     !print*,date,wd,ws,iceil,vis,tmp,dp,slp
     read(fix,'(i4,11x,a12,33x,f3.0,a1,a1,f4.1,a1,i5,a1,2x,f6.3,a1,2x,f5.1,a1,f5.1,a1,f5.1,a1)') nchar,date,wd,awdqc,wdtype,ws,awsqc,iceil,aceilqc,vis,avisqc,tmp,atempqc,dp,adpqc,slp,aslpqc

     !check how many valid fields the record has:
     n_valid=0
     if ( wd    <  999   ) n_valid=n_valid+1
     if ( ws    <  999.9 ) n_valid=n_valid+1
     if ( iceil < 99999  ) n_valid=n_valid+1
     if ( vis   <  999.9 ) n_valid=n_valid+1
     if ( tmp   <  999.9 ) n_valid=n_valid+1
     if (  dp   <  999.9 ) n_valid=n_valid+1
     if ( slp   < 9999.9 ) n_valid=n_valid+1

     !from SMERGE: ---------------------------------------------------
     !convert ceiling height from m to 100s of ft
     if(iceil.gt.99998.or.ishqc(aceilqc).eq.9) then
       iceil=9999
     else
       if(iceil.eq.22000) then
         iceil=999
       else
         iceil=nint(float(iceil)/30.48)
       endif
     endif
     !convert visibility from km to miles
     if(vis.gt.999.98.or.ishqc(avisqc).eq.9) then
       vis=99999.
     else
       vis=vis/1.6093
     endif
     !check temperature and calculate humidity if possible
      if(tmp.gt.999.8.or.ishqc(atempqc).eq.9) then
        tmp=9999.
        irh=9999
      else
        tmp=tmp+273.15
        if(dp.gt.999.8.or.ishqc(adpqc).eq.9) then
          irh=9999
        else
          !compute relative humidity
          tmpf=(tmp-273.15)*1.8+32.
          dpf=dp*1.8+32.
          irh=nint(100.*(((173.-0.1*tmpf+dpf)/(173.+0.9*tmpf))**8))
          irh=max(min(irh,100),1)
        endif
      endif

     !check sea-level pressure
     if(slp.gt.9999.8.or.ishqc(aslpqc).eq.9) slp=9999.

     !end from SMERGE----------------------------------------------

     !decide if valid
     if ( n_valid > 5 ) then
        o%valid=.true.

        o%date=strptime(date,"%Y%m%d%H%M") !record date YYYYMMDDHHmm
        o%wd   =wd /1.0              !deg.
        o%ws   =ws                   !m/s*10 -> m/s
        o%iceil=iceil                !!iceil/1000    !  m    -> km
        o%vis  =vis                  !!vis/1000.     !  m    -> km            
        o%tmp  =tmp                  !!tmp  + 273.15 !  ºC  -> ºK
        o%irh  =irh                  !
        dp     = dp                  !! dp  + 273.15 !  ºC  -> ºK
        o%pres =slp                  !!slp/10.       !hPa*10 -> hPa

        !Additional flag-associated variables:
        extra=row(106:106+nchar)

        !AA1 (precip)
        j=INDEX(extra,'AA1')
        if (j > 0) then
           pp_t  =atoi(extra(j+3:j+4))
           pp_mm =atoi(extra(j+5:j+8))/10.0
           o%prate=real(pp_mm/pp_t)         !prec rate [mm/hr] 
           !ipcode: (precipitation code)
           if (pp_t > 0 .and. o%prate < 2.5                   ) o%ipcode=1 ;
           if (pp_t > 0 .and. o%prate > 2.5 .and. o%prate< 7.6) o%ipcode=2 ;
           if (pp_t > 0 .and. o%prate > 7.6                   ) o%ipcode=3 ;
           if (pp_t < 0 .and. o%prate < 2.5                   ) o%ipcode=19;
           if (pp_t < 0 .and. o%prate > 2.5 .and. o%prate< 7.6) o%ipcode=20;
           if (pp_t < 0 .and. o%prate > 7.6                   ) o%ipcode=21;
        else
           o%prate=9999.
           o%ipcode=99
        endif
        
        !GF1 (cloud cover)
        j=INDEX(extra,'GF1')
        if ( j > 0) then
           o%icc=atoi(extra(j+3:j+4))         !sky cover [octas]
        else  
           o%icc=99
        end if

        !MA1 (pressure)
        if ( slp == 9999.  ) then
           j=INDEX(extra,'MA1')
           if ( j > 0) then
              o%pres=real(atoi(extra(j+9:j+13)))/10.
           end if
        endif

     else
        o%valid=.false.
     end if

end subroutine

integer function ishqc(aqc)
      !ish variables "quality check" (From SMERGE)
      implicit none
      character(len=*):: aqc
      logical         :: lunknown=.true.
      if(ICHAR(aqc).GE.48 .AND. ICHAR(aqc).LE.57) then
         lunknown=.FALSE.
      else
         !Test for known characters
         if(aqc.EQ.'A') lunknown=.FALSE.
         if(aqc.EQ.'C') lunknown=.FALSE.
         if(aqc.EQ.'I') lunknown=.FALSE.
         if(aqc.EQ.'M') lunknown=.FALSE.
         if(aqc.EQ.'P') lunknown=.FALSE.
         if(aqc.EQ.'R') lunknown=.FALSE.
         if(aqc.EQ.'U') lunknown=.FALSE.
      endif
      if(lunknown) then
         write(*,*)'ERROR in subroutine ISHQC'
         write(*,*)'Data quality flag is not known: ',aqc
         write(*,*)'Current flags are 0-9 and A,C,I,M,P,R,U'
         stop
      endif
      !Current character codes indicate values may be used, so test for restricted integer values
      if(aqc.EQ.'3' .OR. aqc.EQ.'7') then
        ishqc=9
      else
        ishqc=1
      endif
      return
end function


subroutine write_Surf_Dat(oFile, S,sdate,edate)
   implicit none
   character(*)  ,intent(in)  :: oFile
   type(station)  ,intent(in) :: S(:)
   type(datetime),intent(in)  :: sdate,edate
   integer                    :: io=7
   type(observation)          :: O
   integer :: i,t
   integer :: NSTA, N_Hours, Time_zone

   character(16) ::  dataset,dataver
   character(64) ::  datamod
   character(16) :: clat,clon
   character(1)  :: hem,mer        !hemispherio & meridiano
   ! --- Configure output variables
   data dataset/'SURF.DAT'/, dataver/'2.1'/
   data datamod/'Hour Start and End Times with Seconds'/

   NSTA   =size(S(:))
   N_Hours=size(S(1)%O(:))

   open (io, file=trim(oFile), action='write', status='replace')   !Open Output file (UP.DAT)

     !Header:
     write(io,'(2a16,a64)') dataset,dataver,datamod
     write(io,*),'1'           !n-comments   
     write(io,*),'ESPG:4326'
     write(io,'("LL",/,"WGS-84   01-01-2001"/,"DEG")')
     !Start / End time
     write(io,'(a3,sp,i3.2,ss,i2.2)') "UTC",0,0
     write(io,*) sdate%strftime("%Y  %j   %H   %S"), edate%strftime(" %Y  %j   %H   %S"), NSTA 
     !Station coordinates
     !write(io,*) NSTA
     do i=1,NSTA
        hem="N"; if ( S(i)%lat < 0 ) hem="S"
        mer="E"; if ( S(i)%lon < 0 ) mer="W"     
        write(clat,"(f12.4,a1)"),abs(S(i)%lat),hem
        write(clon,"(f12.4,a1)"),abs(S(i)%lon),mer
        write(io,'(i12,a12,2a16,f12.2)')S(i)%id,S(i)%name, clat,clon, S(i)%alt
        !write(io,'(a12,f12.3,f12.3,f12.2)') adjustl(S(i)%id), S(i)%x*1E-03, S(i)%y*1E-03, S(i)%z
     enddo
     !Body:
     do t=1,N_Hours
        !Current date
        write(io,*)S(1)%O(t)%date%strftime("%Y %j %-H %-S"),S(1)%O(t+1)%date%strftime(" %Y %j %-H %-S")
        !format(2(3i4,i5,3x))

        !Station record (for each station)
        do i=1, NSTA
           o=S(i)%O(t)
           !write(io,*) o%ws,o%wd,o%iceil,o%icc,o%tmp,o%irh,o%pres,o%ipcode
           write(io,'(1x,f8.3,1x,f8.3,1x,i4,1x,i4,1x,f8.3,1x,i4,1x,f8.3,1x,i4,1x,f8.3)') o%ws,o%wd,o%iceil,o%icc,o%tmp,o%irh,o%pres,o%ipcode!,o%prate
           !format(1x,f8.3,1x,f8.3,1x,i4,1x,i4,1x,f8.3,1x,i4,1x,f8.3,1x,i4,1x,f8.3), (ws(n),wd(n),iceil(n),icc(n),tempk(n),irh(n),pres(n),ipcode(n),xprate(n),n=1,nssta)
        end do
     end do

   close(io)
end subroutine


end module

