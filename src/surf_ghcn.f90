module prep_surf_ghcn
   !
   ! Module to parse GHCN hourly pipe-separated files to surf.dat for CALMET.
   !
   use datetime_module, only: datetime, timedelta, strptime
   use utils_module

   implicit none
   private
   public ghcn2surf

   integer, parameter :: row_len = 20000
   integer, parameter :: max_fields = 500
   integer, parameter :: field_len = 256

   type observation
      type(datetime) :: date
      real    :: ws=0,wd=0,tmp=0,pres=0,prate=0,vis=0
      integer :: iceil=0,icc=0,irh=0,ipcode=0
      logical :: valid=.false.
   end type

   type station
      character(len=5)               :: name
      integer                        :: id
      integer                        :: time_zone=0
      real                           :: lat, lon, alt
      real(8)                        :: x,y,z
      type(observation), allocatable :: O(:)
      integer :: valid_records=0
      logical :: valid=.false.
   end type

   type ghcn_header
      integer :: station=0, station_name=0, date=0
      integer :: latitude=0, longitude=0, elevation=0
      integer :: temperature=0, temperature_qc=0
      integer :: dew_point=0, dew_point_qc=0
      integer :: station_pressure=0, station_pressure_qc=0
      integer :: sea_pressure=0, sea_pressure_qc=0
      integer :: wind_direction=0, wind_direction_qc=0
      integer :: wind_speed=0, wind_speed_qc=0
      integer :: precipitation=0, precipitation_qc=0
      integer :: relative_humidity=0, relative_humidity_qc=0
      integer :: visibility=0, visibility_qc=0
      integer :: sky_condition=0, sky_condition_qc=0
      integer :: ceiling_height=0, ceiling_height_qc=0
      integer :: sky_cover_1=0, sky_cover_sum_1=0
      integer :: precipitation_3_hour=0, precipitation_3_hour_qc=0
      integer :: precipitation_6_hour=0, precipitation_6_hour_qc=0
      integer :: precipitation_12_hour=0, precipitation_12_hour_qc=0
      integer :: precipitation_24_hour=0, precipitation_24_hour_qc=0
   end type

contains

subroutine ghcn2surf(sdate, edate, time_zone_float, file_list)
    !
    ! PURPOSE: Read GHCN hourly PSV files and produce surf.dat.
    !
    implicit none
    character(*)  ,intent(in)   :: file_list(:)
    character(19) ,intent(in)   :: sdate, edate
    real          ,intent(in)   :: time_zone_float

    integer                     :: NSTA
    type(station) ,allocatable  :: S(:)
    logical                     :: file_exists=.false.
    character(200)              :: iFile
    integer                     :: io8=8
    integer                     :: ios
    character(row_len)          :: row

    type(datetime)   :: current_date_utc, current_date, next_date, start_date, end_date
    type(timedelta)  :: dt
    integer          :: N_Hours, hours
    integer          :: i, t
    integer          :: tz_sign, tz_hours, tz_minutes
    type(timedelta)  :: time_zone
    type(observation):: O
    type(ghcn_header):: H

    print '("Parsing GHCN file to surf.dat..   "/)'

    start_date = strptime(sdate,'%Y-%m-%d %H:%M:%S')
    end_date   = strptime(edate,'%Y-%m-%d %H:%M:%S')

    ! Calculate time zone in the same way as the ISH reader.
    tz_sign    = sign(1,int(time_zone_float))
    tz_hours   = int(time_zone_float)
    tz_minutes = int(time_zone_float - tz_hours)*60
    time_zone  = timedelta(hours=tz_hours, minutes=tz_sign*tz_minutes)

    NSTA = size(file_list)

    ! Requested run length, inclusive of start and end hours.
    dt      = (end_date - start_date)
    N_Hours = int(dt%total_seconds()/(60*60)) + 1

    print*,"  start: ", start_date%strftime("day: %Y %j. hour: %H %S")
    print*,"    end: ", end_date%strftime("day: %Y %j. hour: %H %S")
    print*,"    #NT: ", N_Hours, "(hours)."

    allocate(S(NSTA))

    do i=1,NSTA
       iFile=file_list(i)

       ! Allocate one extra record because Write_Surf_Dat writes hour start/end
       ! pairs as O(t) and O(t+1).
       if ( .not. allocated(S(i)%O)) allocate(S(i)%O(N_Hours+1))
       call initialize_observations(S(i)%O, start_date)

       inquire(file=trim(iFile), exist=file_exists)
       if ( file_exists ) then

          print '("   Reading file: ",(A))', trim(iFile)
          open (io8, file=trim(iFile), action="read", status='old')

          read(io8,'(A)',iostat=ios) row
          if (ios /= 0) then
             print '("ERROR. Empty GHCN file: ",A)', trim(iFile)
             S(i)%valid=.false.
             NSTA=NSTA-1
             close(io8)
             cycle
          endif

          call parse_ghcn_header(row, H)
          if (.not. ghcn_header_is_valid(H)) then
             print '("ERROR. Required GHCN columns were not found in: ",A)', trim(iFile)
             S(i)%valid=.false.
             NSTA=NSTA-1
             close(io8)
             cycle
          endif

          read(io8,'(A)',iostat=ios) row
          if (ios /= 0) then
             print '("ERROR. No GHCN records found in: ",A)', trim(iFile)
             S(i)%valid=.false.
             NSTA=NSTA-1
             close(io8)
             cycle
          endif

          ! Read station metadata from the first data row.
          call read_ghcn_station(row, H, i, S(i))

          print '("   Station ID: ",i7," (Lat,Lon,Alt:",f7.3,2x,f8.4,2x,f5.1,")"/)', &
                 S(i)%id, S(i)%lat, S(i)%lon, S(i)%alt

          current_date = ghcn_record_date(row, H) + time_zone
          next_date    = start_date
          do while ( current_date <= end_date )

             current_date_utc = ghcn_record_date(row, H)
             current_date     = current_date_utc + time_zone
             if ( current_date > end_date ) exit
             if ( current_date >= start_date ) then

                print*,"current: ",current_date_utc%strftime("%Y-%j %H:%S:00 UTC"),"   ", &
                       current_date%strftime("%Y-%j %H:%S:00 LST")

                call read_ghcn_record(row, H, O)
                O%date = O%date + time_zone

                dt = current_date - start_date
                t  = int(dt%total_seconds()/3600) + 1

                if ( t > N_Hours ) exit
                if ( O%valid ) then
                   S(i)%O(t) = O
                   S(i)%valid_records = S(i)%valid_records + 1

                   if ( current_date > next_date ) then
                      print '("Warning: No records for this hour: ",A,".")', next_date%isoformat()
                      dt    = current_date - next_date
                      hours = int(dt%total_seconds()/3600)
                      S(i)%O(t-hours:t) = O
                   endif
                endif

                next_date = current_date + timedelta(hours=1)
             endif

             read(io8,'(A)',iostat=ios) row
             if (ios /= 0) exit
          enddo
          close(io8)

          if ( S(i)%valid_records .le. 0.3*N_Hours ) then
             print '("ERROR. No enough valid records for the dates requested found in provided file.")'
             S(i)%valid=.false.
             NSTA=NSTA-1
          else
             S(i)%valid=.true.
          endif
       else
          print '("ERROR. File not found: ",A)', trim(iFile)
          S(i)%valid=.false.
          NSTA=NSTA-1
       endif

    enddo

    call Write_Surf_dat('surf.dat', S, start_date, end_date, tz_hours, tz_minutes)

    print '(/"surf.dat created succesfully?      ")'
    print '("-----------------------------------"/)'

end subroutine


subroutine read_ghcn_record(row, H, o)
     !
     ! Convert one GHCN hourly PSV row into the common SURF observation type.
     ! GHCN values in the sample are already in meteorological units:
     ! temperature/dew point [C], pressure [hPa], wind speed [m/s],
     ! visibility [km], cloud base/ceiling [m], and precipitation [mm].
     !
     implicit none
     character(len=*), intent(in)    :: row
     type(ghcn_header), intent(in)   :: H
     type(observation), intent(inout):: o

     character(field_len), allocatable :: F(:)
     integer              :: nfields, n_valid
     real                 :: wd, ws, vis, tmp, dp, slp, stp, pp_mm, pp_t
     integer              :: iceil, irh, icc

     allocate(F(max_fields))
     call split_psv(row, F, nfields)
     call reset_observation(o)

     wd    = read_real_field(F, nfields, H%wind_direction, 9999.)
     ws    = read_real_field(F, nfields, H%wind_speed, 9999.)
     vis   = read_real_field(F, nfields, H%visibility, 9999.)
     tmp   = read_real_field(F, nfields, H%temperature, 9999.)
     dp    = read_real_field(F, nfields, H%dew_point, 9999.)
     stp   = read_real_field(F, nfields, H%station_pressure, 9999.)
     slp   = read_real_field(F, nfields, H%sea_pressure, 9999.)
     iceil = read_int_field (F, nfields, H%ceiling_height, 999999)
     irh   = read_int_field (F, nfields, H%relative_humidity, 9999)
     icc   = read_cloud_cover(F, nfields, H)

     ! Calm winds are commonly stored as wind_direction=999 with wind_speed=0.
     if (ws == 0.) wd = 0.

     ! Apply simple quality filtering. Blank quality flags are accepted because
     ! many derived GHCN fields use them; restricted ISH-like flags are rejected.
     if (.not. ghcnqc(field_value(F,nfields,H%wind_direction_qc))) wd = 9999.
     if (.not. ghcnqc(field_value(F,nfields,H%wind_speed_qc)))     ws = 9999.
     if (.not. ghcnqc(field_value(F,nfields,H%visibility_qc)))     vis = 9999.
     if (.not. ghcnqc(field_value(F,nfields,H%temperature_qc)))    tmp = 9999.
     if (.not. ghcnqc(field_value(F,nfields,H%dew_point_qc)))      dp = 9999.
     if (.not. ghcnqc(field_value(F,nfields,H%sea_pressure_qc)))   slp = 9999.
     if (.not. ghcnqc(field_value(F,nfields,H%station_pressure_qc))) stp = 9999.
     if (.not. ghcnqc(field_value(F,nfields,H%ceiling_height_qc))) iceil = 999999
     if (.not. ghcnqc(field_value(F,nfields,H%relative_humidity_qc))) irh = 9999

     n_valid=0
     if ( wd    <  999.  ) n_valid=n_valid+1
     if ( ws    < 9999.  ) n_valid=n_valid+1
     if ( iceil < 999999 ) n_valid=n_valid+1
     if ( vis   < 9999.  ) n_valid=n_valid+1
     if ( tmp   < 9999.  ) n_valid=n_valid+1
     if ( dp    < 9999.  ) n_valid=n_valid+1
     if ( slp   < 9999. .or. stp < 9999. ) n_valid=n_valid+1

     ! Convert ceiling height from meters to hundreds of feet. GHCN/ISH use
     ! 22000 m for unlimited ceiling.
     if (iceil >= 999999) then
        iceil = 9999
     elseif (iceil == 22000) then
        iceil = 999
     else
        iceil = nint(real(iceil)/30.48)
     endif

     ! Convert visibility from km to miles, matching the ISH reader.
     if (vis >= 9999.) then
        vis = 99999.
     else
        vis = vis/1.6093
     endif

     ! Temperature is written to SURF.dat as Kelvin. Relative humidity is read
     ! directly when available, otherwise estimated from temperature/dew point.
     if (tmp >= 9999.) then
        tmp = 9999.
        irh = 9999
     else
        if (irh >= 9999 .and. dp < 9999.) irh = calc_rh(tmp, dp)
        tmp = tmp + 273.15
     endif

     if (slp >= 9999.) slp = stp
     if (slp >= 9999.) slp = 9999.

     if (n_valid > 5) then
        o%valid = .true.
        o%date  = ghcn_record_date(row, H)
        o%wd    = wd
        o%ws    = ws
        o%iceil = iceil
        o%vis   = vis
        o%tmp   = tmp
        o%irh   = max(min(irh,100),1)
        if (irh >= 9999) o%irh = 9999
        o%pres  = slp
        o%icc   = icc

        call read_precip(F, nfields, H, pp_mm, pp_t)
        if (pp_t > 0.) then
           o%prate = pp_mm/pp_t
           o%ipcode = precip_code(o%prate)
        else
           o%prate = 9999.
           o%ipcode = 99
        endif
     else
        o%valid = .false.
     endif

end subroutine


subroutine parse_ghcn_header(row, H)
   implicit none
   character(len=*), intent(in)  :: row
   type(ghcn_header), intent(out):: H
   character(field_len), allocatable :: F(:)
   integer                       :: nfields

   allocate(F(max_fields))
   call split_psv(row, F, nfields)

   H%station       = find_field(F, nfields, 'STATION')
   H%station_name  = find_field(F, nfields, 'Station_name')
   H%date          = find_field(F, nfields, 'DATE')
   H%latitude      = find_field(F, nfields, 'LATITUDE')
   H%longitude     = find_field(F, nfields, 'LONGITUDE')
   H%elevation     = find_field(F, nfields, 'ELEVATION')

   H%temperature      = find_field(F, nfields, 'temperature')
   H%temperature_qc   = find_field(F, nfields, 'temperature_Quality_Code')
   H%dew_point        = find_field(F, nfields, 'dew_point_temperature')
   H%dew_point_qc     = find_field(F, nfields, 'dew_point_temperature_Quality_Code')
   H%station_pressure = find_field(F, nfields, 'station_level_pressure')
   H%station_pressure_qc = find_field(F, nfields, 'station_level_pressure_Quality_Code')
   H%sea_pressure     = find_field(F, nfields, 'sea_level_pressure')
   H%sea_pressure_qc  = find_field(F, nfields, 'sea_level_pressure_Quality_Code')
   H%wind_direction   = find_field(F, nfields, 'wind_direction')
   H%wind_direction_qc= find_field(F, nfields, 'wind_direction_Quality_Code')
   H%wind_speed       = find_field(F, nfields, 'wind_speed')
   H%wind_speed_qc    = find_field(F, nfields, 'wind_speed_Quality_Code')
   H%precipitation    = find_field(F, nfields, 'precipitation')
   H%precipitation_qc = find_field(F, nfields, 'precipitation_Quality_Code')
   H%relative_humidity= find_field(F, nfields, 'relative_humidity')
   H%relative_humidity_qc = find_field(F, nfields, 'relative_humidity_Quality_Code')
   H%visibility       = find_field(F, nfields, 'visibility')
   H%visibility_qc    = find_field(F, nfields, 'visibility_Quality_Code')
   H%sky_condition    = find_field(F, nfields, 'sky_condition')
   H%sky_condition_qc = find_field(F, nfields, 'sky_condition_Quality_Code')
   H%ceiling_height   = find_field(F, nfields, 'ceiling_height')
   H%ceiling_height_qc= find_field(F, nfields, 'ceiling_height_Quality_Code')
   H%sky_cover_1      = find_field(F, nfields, 'sky_cover_layer_1')
   H%sky_cover_sum_1  = find_field(F, nfields, 'sky_cover_summation_1')
   H%precipitation_3_hour  = find_field(F, nfields, 'precipitation_3_hour')
   H%precipitation_3_hour_qc = find_field(F, nfields, 'precipitation_3_hour_Quality_Code')
   H%precipitation_6_hour  = find_field(F, nfields, 'precipitation_6_hour')
   H%precipitation_6_hour_qc = find_field(F, nfields, 'precipitation_6_hour_Quality_Code')
   H%precipitation_12_hour = find_field(F, nfields, 'precipitation_12_hour')
   H%precipitation_12_hour_qc = find_field(F, nfields, 'precipitation_12_hour_Quality_Code')
   H%precipitation_24_hour = find_field(F, nfields, 'precipitation_24_hour')
   H%precipitation_24_hour_qc = find_field(F, nfields, 'precipitation_24_hour_Quality_Code')
end subroutine


logical function ghcn_header_is_valid(H)
   implicit none
   type(ghcn_header), intent(in) :: H
   ghcn_header_is_valid = H%station > 0 .and. H%date > 0 .and. &
                          H%latitude > 0 .and. H%longitude > 0 .and. H%elevation > 0
end function


subroutine read_ghcn_station(row, H, idx, S)
   implicit none
   character(len=*), intent(in)  :: row
   type(ghcn_header), intent(in) :: H
   integer, intent(in)           :: idx
   type(station), intent(inout)  :: S
   character(field_len), allocatable :: F(:)
   integer                       :: nfields

   allocate(F(max_fields))
   call split_psv(row, F, nfields)

   write(S%name,'("SS",i2.2)') idx
   S%id  = station_numeric_id(field_value(F, nfields, H%station))
   S%lat = read_real_field(F, nfields, H%latitude, 9999.)
   S%lon = read_real_field(F, nfields, H%longitude, 9999.)
   S%alt = read_real_field(F, nfields, H%elevation, 9999.)
end subroutine


type(datetime) function ghcn_record_date(row, H)
   implicit none
   character(len=*), intent(in)  :: row
   type(ghcn_header), intent(in) :: H
   character(field_len), allocatable :: F(:)
   integer                       :: nfields

   allocate(F(max_fields))
   call split_psv(row, F, nfields)
   ghcn_record_date = strptime(field_value(F, nfields, H%date), "%Y-%m-%dT%H:%M:%S")
end function


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

   O%ws=0.
   O%wd=0.
   O%tmp=9999.
   O%pres=9999.
   O%prate=9999.
   O%vis=99999.
   O%iceil=9999
   O%icc=99
   O%irh=9999
   O%ipcode=99
   O%valid=.false.
end subroutine


subroutine split_psv(row, fields, nfields)
   implicit none
   character(len=*), intent(in)  :: row
   character(len=*), intent(out) :: fields(:)
   integer, intent(out)          :: nfields
   integer                       :: i, start, finish, row_end

   fields = ''
   nfields = 0
   start = 1
   row_end = len_trim(row)

   do i=1,row_end+1
      if (i > row_end .or. row(i:i) == '|') then
         nfields = nfields + 1
         if (nfields <= size(fields)) then
            finish = i - 1
            if (finish >= start) fields(nfields) = adjustl(row(start:finish))
         endif
         start = i + 1
      endif
   enddo
end subroutine


integer function find_field(fields, nfields, name)
   implicit none
   character(len=*), intent(in) :: fields(:)
   integer, intent(in)          :: nfields
   character(len=*), intent(in) :: name
   integer                      :: i

   find_field = 0
   do i=1,min(nfields,size(fields))
      if (trim(fields(i)) == trim(name)) then
         find_field = i
         return
      endif
   enddo
end function


character(field_len) function field_value(fields, nfields, idx)
   implicit none
   character(len=*), intent(in) :: fields(:)
   integer, intent(in)          :: nfields, idx

   field_value = ''
   if (idx > 0 .and. idx <= nfields .and. idx <= size(fields)) field_value = trim(fields(idx))
end function


real function read_real_field(fields, nfields, idx, missing)
   implicit none
   character(len=*), intent(in) :: fields(:)
   integer, intent(in)          :: nfields, idx
   real, intent(in)             :: missing
   integer                      :: ios
   character(field_len)         :: value

   read_real_field = missing
   value = field_value(fields,nfields,idx)
   if (len_trim(value) <= 0) return
   read(value,*,iostat=ios) read_real_field
   if (ios /= 0) read_real_field = missing
end function


integer function read_int_field(fields, nfields, idx, missing)
   implicit none
   character(len=*), intent(in) :: fields(:)
   integer, intent(in)          :: nfields, idx, missing
   integer                      :: ios
   character(field_len)         :: value

   read_int_field = missing
   value = field_value(fields,nfields,idx)
   if (len_trim(value) <= 0) return
   read(value,*,iostat=ios) read_int_field
   if (ios /= 0) read_int_field = missing
end function


logical function ghcnqc(aqc)
   !
   ! Keep the ISH convention that quality codes 3 and 7 are restricted.
   ! Blank flags are accepted because several GHCN derived columns leave their
   ! quality flag empty even when the value is usable.
   !
   implicit none
   character(len=*), intent(in) :: aqc
   character(field_len)         :: q

   q = trim(aqc)
   ghcnqc = .true.
   if (len(q) <= 0) return
   if (q == '3' .or. q == '7') ghcnqc = .false.
end function


integer function calc_rh(tmp, dp)
   implicit none
   real, intent(in) :: tmp, dp
   real             :: tmpf, dpf

   tmpf = tmp*1.8 + 32.
   dpf  = dp*1.8 + 32.
   calc_rh = nint(100.*(((173.-0.1*tmpf+dpf)/(173.+0.9*tmpf))**8))
   calc_rh = max(min(calc_rh,100),1)
end function


integer function read_cloud_cover(fields, nfields, H)
   implicit none
   character(len=*), intent(in) :: fields(:)
   integer, intent(in)          :: nfields
   type(ghcn_header), intent(in) :: H
   character(field_len)         :: sky
   integer                      :: pos, ios

   read_cloud_cover = 99

   sky = field_value(fields, nfields, H%sky_condition)
   if (len_trim(sky) >= 2) then
      read(sky(1:2),*,iostat=ios) read_cloud_cover
      if (ios == 0) return
   endif

   sky = field_value(fields, nfields, H%sky_cover_sum_1)
   pos = index(sky, ':')
   if (pos > 0 .and. len_trim(sky) >= pos+2) then
      read(sky(pos+1:pos+2),*,iostat=ios) read_cloud_cover
      if (ios == 0) return
   endif

   sky = field_value(fields, nfields, H%sky_cover_1)
   pos = index(sky, ':')
   if (pos > 0 .and. len_trim(sky) >= pos+2) then
      read(sky(pos+1:pos+2),*,iostat=ios) read_cloud_cover
   endif
end function


subroutine read_precip(fields, nfields, H, pp_mm, pp_t)
   implicit none
   character(len=*), intent(in) :: fields(:)
   integer, intent(in)          :: nfields
   type(ghcn_header), intent(in):: H
   real, intent(out)            :: pp_mm, pp_t

   pp_mm = 9999.
   pp_t  = 0.

   if (ghcnqc(field_value(fields,nfields,H%precipitation_qc))) then
      pp_mm = read_real_field(fields, nfields, H%precipitation, 9999.)
      if (pp_mm < 9999.) then
         pp_t = 1.
         return
      endif
   endif

   call read_precip_period(fields, nfields, H%precipitation_3_hour,  H%precipitation_3_hour_qc,  3.,  pp_mm, pp_t)
   if (pp_t > 0.) return
   call read_precip_period(fields, nfields, H%precipitation_6_hour,  H%precipitation_6_hour_qc,  6.,  pp_mm, pp_t)
   if (pp_t > 0.) return
   call read_precip_period(fields, nfields, H%precipitation_12_hour, H%precipitation_12_hour_qc, 12., pp_mm, pp_t)
   if (pp_t > 0.) return
   call read_precip_period(fields, nfields, H%precipitation_24_hour, H%precipitation_24_hour_qc, 24., pp_mm, pp_t)
end subroutine


subroutine read_precip_period(fields, nfields, value_idx, qc_idx, period, pp_mm, pp_t)
   implicit none
   character(len=*), intent(in) :: fields(:)
   integer, intent(in)          :: nfields, value_idx, qc_idx
   real, intent(in)             :: period
   real, intent(out)            :: pp_mm, pp_t

   pp_mm = 9999.
   pp_t  = 0.
   if (.not. ghcnqc(field_value(fields,nfields,qc_idx))) return
   pp_mm = read_real_field(fields, nfields, value_idx, 9999.)
   if (pp_mm < 9999.) pp_t = period
end subroutine


integer function precip_code(prate)
   implicit none
   real, intent(in) :: prate

   precip_code = 99
   if (prate < 2.5) then
      precip_code = 1
   elseif (prate < 7.6) then
      precip_code = 2
   else
      precip_code = 3
   endif
end function


integer function station_numeric_id(station_name)
   implicit none
   character(len=*), intent(in) :: station_name
   integer                      :: i, start, ios

   station_numeric_id = 0
   start = 0
   do i=len_trim(station_name),1,-1
      if (station_name(i:i) < '0' .or. station_name(i:i) > '9') exit
      start = i
   enddo

   if (start > 0) then
      read(station_name(start:len_trim(station_name)),*,iostat=ios) station_numeric_id
      if (ios /= 0) station_numeric_id = 0
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

   character(16) ::  dataset,dataver
   character(64) ::  datamod
   character(16) :: clat,clon
   character(1)  :: hem,mer

   data dataset/'SURF.DAT'/, dataver/'2.1'/
   data datamod/'Hour Start and End Times with Seconds'/

   NSTA    = count(S(:)%valid)
   N_Hours = size(S(1)%O(:)) - 1

   if (NSTA <= 0) then
      print '("ERROR. No valid GHCN stations available for surf.dat.")'
      return
   endif

   open (io, file=trim(oFile), action='write', status='replace')

     write(io,'(2a16,a64)') dataset,dataver,datamod
     write(io,*) '1'
     write(io,*) 'ESPG:4326'
     write(io,'("LL",/,"WGS-84   01-01-2001"/,"DEG")')
     write(io,'(a3,sp,i3.2,ss,i2.2)') "UTC",tz_hour,tz_min
     write(io,*) sdate%strftime("%Y  %j   %H   %S"), edate%strftime(" %Y  %j   %H   %S"), NSTA

     do i=1,size(S)
        if (.not. S(i)%valid) cycle
        hem="N"; if ( S(i)%lat < 0 ) hem="S"
        mer="E"; if ( S(i)%lon < 0 ) mer="W"
        write(clat,"(f12.4,a1)") abs(S(i)%lat),hem
        write(clon,"(f12.4,a1)") abs(S(i)%lon),mer
        write(io,'(i12,a12,2a16,f12.2)') S(i)%id,S(i)%name, clat,clon, S(i)%alt
     enddo

     do t=1,N_Hours
        write(io,*) S(1)%O(t)%date%strftime("%Y %j %-H %-S"), &
                    S(1)%O(t+1)%date%strftime(" %Y %j %-H %-S")

        do i=1,size(S)
           if (.not. S(i)%valid) cycle
           o=S(i)%O(t)
           write(io,'(1x,f8.3,1x,f8.3,1x,i4,1x,i4,1x,f8.3,1x,i4,1x,f8.3,1x,i4,1x,f8.3)') &
                    o%ws,o%wd,o%iceil,o%icc,o%tmp,o%irh,o%pres,o%ipcode
        enddo
     enddo

   close(io)
end subroutine


end module
