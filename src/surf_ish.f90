module prep_surf_ish
   !
   ! Module to parse ISH files to SURF.dat & PRECIP.dat files for CALMET
   !
   use datetime_module, only: datetime, timedelta, strptime
   use utils_module

   implicit none
   private
   public ish2surf

   integer, parameter :: row_len = 20000

   real,    parameter :: rmiss=9999.,vmiss=99999.
   integer, parameter :: imiss=9999

   type observation
      type(datetime) :: date             !datetime of observation
      real    :: ws=rmiss,wd=rmiss,tmp=rmiss,pres=rmiss,prate=0.0,vis=vmiss
      integer :: iceil=imiss,icc=imiss,irh=imiss,ipcode=0
      logical :: valid=.false.
   end type

   type station
      character(len=5)               :: name
      integer                        :: id=0
      integer                        :: time_zone=0
      real                           :: lat=rmiss, lon=rmiss, alt=rmiss
      real(8)                        :: x=0.d0,y=0.d0,z=0.d0
      type(observation), allocatable :: O(:)                 !observation array
      integer :: valid_records=0                             !number of hours with records
      logical :: valid=.false.                               !is a valid station?
   end type

contains

subroutine ish2surf(sdate, edate, time_zone_float, file_list)
    !
    ! PURPOSE: Read ISH file and produce surf.dat and precip.dat files
    !
    implicit none
    character(*)  ,intent(in)   :: file_list(:)           !Input Surfac Met. (ISH) files list.
    character(19) ,intent(in)   :: sdate, edate
    real          ,intent(in)   :: time_zone_float        !Time zone in Hours

    integer                     :: NSTA                   !Number of stations
    type(station) ,allocatable  :: S(:)                   !record buffer
    logical                     :: file_exists=.false.    !flag file existence
    character(200)              :: iFile     
    integer                     :: io8=8                  !file unit number for IO 
    integer                     :: ios                    !I/O status indice
    character(row_len)          :: row

    type(datetime)   :: current_date_utc, current_date, start_date, end_date
    type(timedelta)  :: dt
    integer          :: N_Hours                           !Run length [hours]
    integer          :: i, t
    integer          :: tz_hours, tz_minutes !tz_sign, 
    type(timedelta)  :: time_zone

    type(observation):: O

    Print '("Parsing ISH file to surf.dat..   "/)'

    start_date = strptime(sdate,'%Y-%m-%d %H:%M:%S')
    end_date   = strptime(edate,'%Y-%m-%d %H:%M:%S')

    ! Calculate time zone
    tz_hours   = int(time_zone_float)
    tz_minutes = sign(int(time_zone_float - tz_hours)*60, tz_hours)
    time_zone  = timedelta(hours=tz_hours, minutes=tz_minutes)

    NSTA = size(file_list)

    !Calculate run length [hours]
    dt      = (end_date - start_date)
    N_Hours = int(dt%total_seconds()/(60*60)) + 1

    Print*,"  start: ", start_date%strftime("day: %Y %j. hour: %H %S")
    Print*,"    end: ",   end_date%strftime("day: %Y %j. hour: %H %S")
    Print*,"    #NT: ",  N_hours, "(hours)."

    allocate(S(NSTA))

    station_loop: do i=1,NSTA
       iFile=file_list(i)

       !Allocate one extra record because Write_Surf_Dat writes hour start/end pairs as O(t) and O(t+1).
       if ( .not. allocated(S(i)%O)) allocate(S(i)%O(N_Hours+1))
       call initialize_observations(S(i)%O, start_date)

       inquire(file=trim(iFile), exist=file_exists)
       if ( file_exists ) then

          Print '("   Reading file: ",(A))', trim(iFile)
          open (io8, file=trim(iFile), action="read", status='old') 

          !Read header (first record)
          read(io8,'(A)',iostat=ios)row
          if (ios /= 0) then
             Print '("ERROR. Empty or unreadable file: ",A)', trim(iFile)
             close(io8)
             S(i)%valid=.false.
             cycle
          endif

          !Parse station header                                     !Read ISH first row to get staton data.
          write(S(i)%name,'("SS",i2.2)') i                          !station-id
          S(i)%id =atoi(row(5:10))                                  !USAF/WMO station-id
          S(i)%lat=real(atoi(row(29:34))/1000.)                     !latitude
          S(i)%lon=real(atoi(row(35:41))/1000.)                     !longitude
          S(i)%alt=real(atoi(row(47:51))/1.)                        !base altitude
          S(i)%valid=.true.
          Print '("   Station ID: ",i7," (Lat,Lon,Alt:",f7.3,2x,f8.4,2x,f5.1,")")', S(i)%id, S(i)%lat, S(i)%lon, S(i)%alt

          records_loop: do while (ios == 0)

             call read_ISH_record(row, O, S(i)%alt)
             current_date_utc = O%date
             current_date     = current_date_utc + time_zone
             O%date           = current_date

             if ( current_date >  end_date ) exit
             if ( current_date >= start_date .and. O%valid ) then

                !Print*,"      current: ",current_date_utc%strftime("%Y-%j %H:%S:00 UTC"), "   ",current_date%strftime("%Y-%j %H:%S:00 LST")

                dt = current_date - start_date
                t  = int(dt%total_seconds()/3600) !+ 1

                if ( t >= 1 .and. t <= N_Hours ) then
                   if (.not. S(i)%O(t)%valid) S(i)%valid_records=S(i)%valid_records + 1
                   S(i)%O(t) = O
                endif
             endif

             read(io8,'(A)',iostat=ios) row
             if (ios /= 0) exit 
          end do records_loop
          close(io8)

          if ( S(i)%valid_records .le. 0.3*N_Hours ) then
             print '("     WARNING! No enough valid records for the dates requested found in provided file.")'
             S(i)%valid=.false.
          else
             S(i)%valid=.true.
          endif

       else
          print '("ERROR. File not found: ",A)',trim(iFile)
          S(i)%valid=.false.

       endif

    enddo station_loop 

    !Write output:
    call Write_Surf_dat('surf.dat', S, start_date, end_date, tz_hours, tz_minutes)

    print '(/A," created succesfully?      ")', 'surf.dat'
    print '("-----------------------------------"/)'

end subroutine



subroutine read_ish_record(row, o, elev)
     !                                                                         
     ! Convert one ISH row record into the common SURF observation type.
     !
     implicit none
     character(len=*), intent(in)     :: row
     type(observation), intent(inout) :: o
     real, intent(in) :: elev
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

     integer :: nchar, last_extra
     character(len=*), parameter :: fixfmt = '(i4,11x,a12,33x,f3.0,a1,a1,f4.1,a1,i5,a1,2x,f6.3,a1,2x,f5.1,a1,f5.1,a1,f5.1,a1)'
     character(12) :: date !yyyymmddhhmin
     real   :: wd,ws,vis,tmp,dp,slp,altim,pres
     character :: awdqc,awsqc,aceilqc,avisqc,atempqc,adpqc,aslpqc,wdtype
     character :: alpqc,askcqc,aaltqc,apresqc,apwqc
     integer:: iceil,irh, lpdur, itskc, ioskc, ipw

     real :: tmpf, dpf    !vars used to compute rh
     real :: pp_mm        !vars used to compute prate & pcode

     !Set all observation values to null or missing
     call reset_observation(o)

     !Fixed position values:
     fix=row(1:105)
     !read(fix,"(i4,11x,a12,33x,f3.0,2x,f4.1,1x,i5,3x,f6.3,3x,f5.1,1x,f5.1,1x,f5.1,1x)"),nchar,date,wd,ws,iceil,vis,tmp,dp,slp
     read(fix,fixfmt) nchar,date,wd,awdqc,wdtype,ws,awsqc,iceil,aceilqc,vis,avisqc,tmp,atempqc,dp,adpqc,slp,aslpqc

     !check how many valid fields the record has:
     n_valid=0
     if ( wd    <  999.  .and. qcheck(awdqc)   /= 9 ) n_valid=n_valid+1
     if ( ws    < 9999.  .and. qcheck(awsqc)   /= 9 ) n_valid=n_valid+1
     if ( iceil < 999999 .and. qcheck(aceilqc) /= 9 ) n_valid=n_valid+1
     if ( vis   < 9999.  .and. qcheck(avisqc)  /= 9 ) n_valid=n_valid+1
     if ( tmp   < 9999.  .and. qcheck(atempqc) /= 9 ) n_valid=n_valid+1
     if ( dp    < 9999.  .and. qcheck(adpqc)   /= 9 ) n_valid=n_valid+1
     if ( slp   < 9999.  .and. qcheck(aslpqc)  /= 9 ) n_valid=n_valid+1

     !from SMERGE----------------------------------------------
     !wind hierarchy: calm, variable/missing, calm speed, then QC.
     if ( wdtype.eq.'C' .or. wdtype.eq.'c'     ) then
        ws=0.;    wd=0.
     elseif ( wdtype.eq.'V' .or. wdtype.eq.'v' ) then
        ws=rmiss; wd=rmiss
     elseif ( ws.lt.0.1 .and. qcheck(awsqc).ne.9) then
        ws=0.;    wd=0.
     elseif ( wd.gt.998..or.qcheck(awdqc) .eq. 9) then
        ws=rmiss; wd=rmiss
     elseif ( ws.gt.998..or.qcheck(awsqc) .eq. 9) then
        ws=rmiss; wd=rmiss
     elseif ( wdtype.eq.'9') then
        ws=rmiss; wd=rmiss
     endif

     !convert ceiling height from m to 100s of ft.
     if(iceil>=22000.or.qcheck(aceilqc).eq.9) then
       iceil=imiss
     else
         iceil=nint(float(iceil)/30.48)
     endif

     !convert visibility from km to miles
     if (vis>999.98.or.qcheck(avisqc).eq.9) then
       vis=vmiss
     else
       vis=vis/1.6093
     endif

     !check temperature and compute humidity if possible
     if(tmp>999.8.or.qcheck(atempqc).eq.9) then
       tmp=rmiss
       irh=imiss
     else
       tmp=tmp+273.15
       if(dp>999.8.or.qcheck(adpqc).eq.9) then
         irh=imiss
       else
         !compute relative humidity
         tmpf=(tmp-273.15)*1.8+32.
         dpf=dp*1.8+32.
         irh=nint(100.*(((173.-0.1*tmpf+dpf)/(173.+0.9*tmpf))**8))
         irh=max(min(irh,100),1)
       endif
     endif

     !check sea-level pressure
     if (slp >= 9999.8.or.qcheck(aslpqc).eq.9) slp=rmiss

     !end from SMERGE----------------------------------------------

     o%date=strptime(date(1:10)//"00","%Y%m%d%H%M")
     if (atoi(date(11:12)) > 0) o%date=o%date+timedelta(hours=1)

     !decide if valid
     if ( n_valid > 0 ) then
        o%valid=.true.

        o%wd   =wd /1.0              !deg.
        o%ws   =ws                   !m/s
        o%iceil=iceil                !!iceil/1000    !  m    -> km
        o%vis  =vis                  !!vis/1000.     !  m    -> km            
        o%tmp  =tmp                  !!tmp  + 273.15 !  ºC  -> ºK
        o%irh  =irh                  !
        o%pres =slp

        !Additional flag-associated variables:
        extra=' '
        last_extra=min(len_trim(row),105+nchar)
        if (last_extra > 105) extra(1:last_extra-105)=row(106:last_extra)

        !Precipitation (AA1)
        j=INDEX(extra,'AA1')
        if (j > 0) then
           lpdur = atoi(extra(j+3:j+4))
           pp_mm = atoi(extra(j+5:j+8))/10.0
           alpqc = extra(j+10:j+10)
           if (lpdur == 1 .and. pp_mm < 999.8 .and. qcheck(alpqc) /= 9) then
              o%prate=pp_mm
              if (o%prate < 0.01) then
                 o%ipcode=0
              elseif (o%tmp < 9998. .and. o%tmp >= 273.15) then
                 o%ipcode=1
              elseif (o%tmp < 9998.) then
                 o%ipcode=20
              endif
           endif
        else
           o%prate  = 0.0 !rmiss
           o%ipcode = 0   !imiss
        endif

        !AW1/MW1 present weather can override precipitation type.
        j=INDEX(extra,'AW1')
        if (j <= 0) j=INDEX(extra,'MW1')
        if (j > 0) then
           ipw=atoi(extra(j+3:j+4))
           apwqc=extra(j+5:j+5)
           if (qcheck(apwqc) /= 9) call set_ipcode_from_weather(ipw, o%tmp, o%ipcode)
        endif
        
        !Cloud cover (GF1)
        j=INDEX(extra,'GF1')
        if ( j > 0) then
           itskc=atoi(extra(j+3:j+4))
           ioskc=atoi(extra(j+5:j+6))
           askcqc=extra(j+7:j+7)
           if (ioskc == 99) ioskc=itskc
           if (ioskc == 99 .or. qcheck(askcqc) == 9) then
              o%icc=imiss
           else
              o%icc=sky_cover_to_tenths(ioskc)
           endif
        else  
           o%icc=imiss
        end if
      
        !Pressure (MA1)
        j=INDEX(extra,'MA1')
        if ( j > 0) then
           altim = real(atoi(extra(j+3:j+7)))/10.
           aaltqc = extra(j+8:j+8)
           pres  = real(atoi(extra(j+9:j+13)))/10.
           apresqc = extra(j+14:j+14)
           if (pres < 9999.8 .and. qcheck(apresqc) /= 9) o%pres=pres
           if (o%pres > 9998. .and. altim < 9999.8 .and. qcheck(aaltqc) /= 9) then
              o%pres=pressure_from_altimeter(altim, elev)
           endif
        endif
        if (o%pres > 9998. .and. slp < 9998.) o%pres=pressure_from_slp(slp, elev, o%tmp)

     else
        o%valid = .false.
     endif

end subroutine

subroutine initialize_observations(O, start_date)
   implicit none
   type(observation), intent(inout) :: O(:)
   type(datetime), intent(in)       :: start_date
   integer                          :: t

   do t=1,size(O)
      call reset_observation(O(t))
      O(t)%date = start_date + timedelta(hours=t-1)
   enddo
end subroutine

subroutine reset_observation(O)
   implicit none
   type(observation), intent(inout) :: O

   O%ws    = rmiss
   O%wd    = rmiss
   O%tmp   = rmiss
   O%pres  = rmiss
   O%prate = 0.0   !rmiss
   O%vis   = imiss
   O%iceil = imiss
   O%icc   = imiss
   O%irh   = imiss
   O%ipcode= 0 !imiss
   O%valid =.false.
end subroutine

integer function qcheck(aqc)
      !"quality check" (From SMERGE)
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
        qcheck=9
      else
        qcheck=1
      endif
      return
end function

integer function sky_cover_to_tenths(iskc)
      implicit none
      integer, intent(in) :: iskc

      select case (iskc)
      case (0)
         sky_cover_to_tenths=0
      case (1)
         sky_cover_to_tenths=1
      case (2)
         sky_cover_to_tenths=3
      case (3)
         sky_cover_to_tenths=4
      case (4)
         sky_cover_to_tenths=5
      case (5)
         sky_cover_to_tenths=6
      case (6)
         sky_cover_to_tenths=8
      case (7)
         sky_cover_to_tenths=9
      case (8)
         sky_cover_to_tenths=10
      case default
         sky_cover_to_tenths=imiss
      end select
end function

subroutine set_ipcode_from_weather(ipw, tempk, ipcode)
      implicit none
      integer, intent(in)    :: ipw
      real,    intent(in)    :: tempk
      integer, intent(inout) :: ipcode

      if (((ipw <= 19) .or. (ipw >= 30 .and. ipw <= 35) .or. &
           (ipw >= 40 .and. ipw <= 49)) .and. ipcode == imiss) then
         ipcode=0
      elseif (ipw == 20 .or. ipw == 21 .or. (ipw >= 23 .and. ipw <= 25) .or. &
              (ipw >= 50 .and. ipw <= 69) .or. (ipw >= 80 .and. ipw <= 84) .or. &
              ipw == 91 .or. ipw == 92) then
         ipcode=1
      elseif (ipw == 22 .or. (ipw >= 36 .and. ipw <= 39) .or. &
              (ipw >= 70 .and. ipw <= 79) .or. ipw == 85 .or. ipw == 86 .or. &
              ipw == 96 .or. ipw == 99) then
         if (ipcode /= 1) ipcode=20
      elseif (tempk < 9998.) then
         if (tempk >= 273.15) then
            ipcode=1
         else
            if (ipcode /= 1) ipcode=20
         endif
      endif
end subroutine

real function pressure_from_altimeter(altim, elev)
      implicit none
      real, intent(in) :: altim, elev
      real :: z

      if (altim > 9998.) then
         pressure_from_altimeter=rmiss
      else
         z=max(elev,0.)
         !pressure_from_altimeter=altim*(1. - 2.25577e-5*z)**5.25588
         pressure_from_altimeter=((altim**0.19026)-8.41717e-5*z)**5.25593

      endif
end function

real function pressure_from_slp(slp, elev, tempk)
      implicit none
      real, intent(in) :: slp, elev, tempk
      real :: tk, z

      if (slp > 9998.) then
         pressure_from_slp=rmiss
      else
         tk=tempk
         if (tk > 9998.) tk=288.15
         z=max(elev,0.)
         pressure_from_slp=slp*exp(-z/(29.3*tk))
      endif
end function

subroutine write_Surf_Dat(oFile, S, sdate, edate, tz_hour, tz_min)
   implicit none
   character(*)  ,intent(in)  :: oFile
   type(station) ,intent(in)  :: S(:)
   type(datetime),intent(in)  :: sdate,edate
   integer       ,intent(in)  :: tz_hour,tz_min
   integer                    :: io=7
   type(observation)          :: O
   integer :: i,t
   integer :: NSTA, N_Hours
   type(datetime) :: bdate,edate1

   character(16) ::  dataset,dataver
   character(64) ::  datamod
   character(16) :: clat,clon
   character(1)  :: hem,mer
   character(len=*), parameter :: recfmt='(1x,f8.3,1x,f8.3,1x,i4,1x,i4,1x,f8.3,1x,i4,1x,f8.3,1x,i4,1x,f8.3)'
   character(len=*), parameter :: recfmt_old='(1x,f8.3,1x,f8.3,1x,i4,1x,i4,1x,f8.3,1x,i4,1x,f8.3,1x,i4)'

   data dataset/'SURF.DAT'/, dataver/'2.1'/
   data datamod/'Hour Start and End Times with Seconds'/

   NSTA=0
   do i=1,size(S)
      if (S(i)%valid) NSTA=NSTA+1
   enddo
   N_Hours=size(S(1)%O(:)) - 1

   open (io, file=trim(oFile), action='write', status='replace')

     !Header
     write(io,'(2a16,a64)') dataset,dataver,datamod
     write(io,*)'1'
     write(io,*)'SURF.DAT file created by CALPREP v1.0' !ESPG:4326'
     write(io,'("LL",/,"WGS-84   01-01-2001"/,"DEG")')
     write(io,'(a3,sp,i3.2,ss,i2.2)') "UTC",tz_hour,abs(tz_min)
     write(io,*) sdate%strftime("%Y  %j   %H   %S"), edate%strftime(" %Y  %j   %H   %S"), NSTA

     !Station coordinates
     do i=1,size(S)
        !if (.not. S(i)%valid) cycle
        hem="N"; if ( S(i)%lat < 0 ) hem="S"
        mer="E"; if ( S(i)%lon < 0 ) mer="W"
        write(clat,"(f12.4,a1)") abs(S(i)%lat),hem
        write(clon,"(f12.4,a1)") abs(S(i)%lon),mer
        write(io,'(i12,a12,2a16,f12.2)') S(i)%id,S(i)%name, clat,clon, S(i)%alt
     enddo

     !Body
     do t=1,N_Hours
        !Current date
        bdate =sdate + timedelta(hours=t-1)
        edate1=bdate + timedelta(hours=1)
        write(io,*)bdate%strftime("%Y %j %-H %-S"),edate1%strftime(" %Y %j %-H %-S")

        !Station record (for each station)
        do i=1, size(S)
           !if (.not. S(i)%valid) cycle
           o=S(i)%O(t)
           write(io,recfmt_old) o%ws,o%wd,o%iceil,o%icc,o%tmp,o%irh,o%pres,o%ipcode
           !write(io,recfmt) o%ws,o%wd,o%iceil,o%icc,o%tmp,o%irh,o%pres,o%ipcode,o%prate
        enddo
     enddo

   close(io)
end subroutine


end module
