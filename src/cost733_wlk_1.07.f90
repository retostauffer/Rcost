subroutine wlk733(dates,ydays,nt,lat,ni,lon,nj,u700,v700,z500,z925,pw,pwclim,nsector,weight,maindirthreshold,ZAMG,ires,rres,verbose)
  ! determine WLK(DWD)-like weather types using ERA40 data
  ! andreas.philipp@geo.uni-augsburg.de
  implicit none

  integer, intent(IN) :: ni, nj, nt, nsector, ZAMG, verbose
  real(kind=8), intent(IN) :: maindirthreshold
  real(kind=8), dimension(nj), intent(IN) :: lon
  real(kind=8), dimension(ni), intent(IN) :: lat
  real(kind=8), dimension(ni,nj,nt), intent(IN) :: u700, v700, z500, z925, pw
  real(kind=8), dimension(366), intent(IN) :: pwclim
  real(kind=8), dimension(ni,nj), intent(IN) :: weight

  integer, dimension(nt,6), intent(INOUT) :: ires
  real(kind=8), dimension(nt,3), intent(INOUT) :: rres

  integer, dimension(nt), intent(IN) :: dates ! Integer dates (YYYYmmdd, for output only)
  integer, dimension(nt), intent(IN) :: ydays ! Day of the year 0,..,365

  integer :: t,x,y
  real(kind=8) :: diflon, diflat
  real(kind=8) :: minlon, maxlon, minlat, maxlat 
  real(kind=8),allocatable :: cyc925(:),cyc500(:)


  integer,allocatable :: mainsector(:)
  real(kind=8) :: ltm(31,12),div(31,12)
  real(kind=8),allocatable :: pwi(:)
  integer, allocatable :: humid(:)

  integer :: nclass
  integer, allocatable :: intclass(:)
  character(len=5), allocatable :: charclass(:)
  character(len=5) :: classstring
  integer :: classnumber,i,c925,c500,hum

  integer, allocatable :: clsize(:)
   
  if ( verbose .ge. 4 ) then
    do i = 1,ni
      write(*,"(999f2.0)")weight(i,:)
    enddo
  endif

  ! Compute grid extent
  ! compute grid and grid increments
  minlon = MINVAL(lon)
  maxlon = MAXVAL(lon)
  minlat = MINVAL(lat)
  maxlat = MAXVAL(lat)
  diflon=(maxlon-minlon) / real(nj - 1)
  diflat=(maxlat-minlat) / real(ni - 1)

  ! Some user readable information
  if ( verbose .gt. 0 ) then
    write(*,"(5xa)") "Domain definition:"
    write(*,"(10X,1A,I8)")              "Timesteps ", nt
    write(*,"(10X,1A,I8,1A,I8,1A)")     "Gridsize  ", ni, "(rows;lat) x  ", nj, "(cols;lon)"
    write(*,"(10X,1A,F8.3,1A,F8.3)")    "Longitude ", minlon, " to ", maxlon
    write(*,"(10X,1A,F8.3,1A,F8.3)")    "Latitude  ", minlat, " to ", maxlat
    write(*,"(10X,1A,F8.3,1A,F8.3,1A)") "Increment ", diflon, " (lon), ", diflat, " (lat)"
  endif

  ! __________________________________________________________________________
  ! DATE INFORMATION FOR RESULT ARRAY
  do t=1,nt
    ires(t,1) = dates(t)
    ires(t,2) = ydays(t)
  enddo

  ! __________________________________________________________________________
  ! MAIN WIND SECTOR
  ! Using 700hPa wind direction based on meridional and zonal wind speeds (u,v)
  allocate(mainsector(nt))
  if ( ZAMG .eq. 1 ) then
     do t=1,nt
        call windsectorZAMG(u700(:,:,t),v700(:,:,t),ni,nj,weight,maindirthreshold,mainsector(t),verbose)
        ires(t,3) = mainsector(t)
     enddo
  else
     do t=1,nt
        call windsector(u700(:,:,t),v700(:,:,t),ni,nj,nsector,weight,maindirthreshold,mainsector(t),verbose)
        ires(t,3) = mainsector(t)
     enddo
  endif

  ! __________________________________________________________________________
  ! CYCLONICITY
  ! FOR Z925
  allocate(cyc925(nt))
  do t=1,nt
     call cycindex(z925(:,:,t),lat,ni,lon,nj,weight, cyc925(t))
     rres(t,1) = cyc925(t)
  enddo

  ! FOR Z500
  allocate(cyc500(nt))
  do t=1,nt
     call cycindex(z500(:,:,t),lat,ni,lon,nj,weight, cyc500(t))
     rres(t,2) = cyc500(t)
  enddo

  ! __________________________________________________________________________
  ! PW-INDEX
  ! AREA MEAN IS PW-INDEX pwi
  allocate(pwi(nt))
  do t=1,nt
     pwi(t)    = SUM(pw(1:ni,1:nj,t) * weight(1:ni,1:nj))/SUM(weight)
     rres(t,3) = pwi(t)
  enddo
  allocate(humid(nt))
  ! DECIDE WET OR DRY
  humid=0 ! Default: dry weather situation
  do t=1,nt
     if ( pwclim(ydays(t)+1) == -999. ) then
        humid = -90 ! pwclim unspecified
     elseif ( pwi(t) > pwclim(ydays(t)+1) ) then
        humid(t)=1 ! Humid weather situation
     endif
  enddo

  ! __________________________________________________________________________
  ! ASSIGN TO TYPES
  allocate(intclass(nt))
  intclass(1:nt)=0
  allocate(charclass(nt))
  charclass(1:nt)="     "

  nclass=(nsector+1)*2*2*2
  
  allocate(clsize(nclass))
  clsize=0

  do t=1,nt

     ! cyclonicity 925 hPa
     charclass(t)(3:3)="A"
     ires(t,4) = 0
     if(cyc925(t)>0.D0)then
        charclass(t)(3:3)="C"
        intclass(t)=intclass(t)+(nclass/2)
        ires(t,4) = 1
     endif

     ! cyclonicity 500 hPa
     charclass(t)(4:4)="A"
     ires(t,5) = 0
     if(cyc500(t)>0.D0)then
        charclass(t)(4:4)="C"
        intclass(t)=intclass(t)+(nclass/4)
        ires(t,5) = 1
     endif
     
     ! humidity
     charclass(t)(5:5)="D"    ! Default try D (dry)
     ires(t,6) = 0
     if ( humid(t) .lt. 0 ) then
        charclass(t)(5:5)="X" ! if humid(t) undefined: X (missing)
        ires(t,6) = -99
     elseif (humid(t)==1)then
        charclass(t)(5:5)="W" ! if humid==1 W (wet)
        intclass(t)=intclass(t)+(nclass/8)
        ires(t,6) = 1
     endif

     ! main wind sector
     !write(charclass(t)(1:2),"(1i2.2)")mainsector(t)
     intclass(t)=intclass(t)+mainsector(t)+1

     ! OUTPUT OF TYPE
     write(*,"(2i4,3X,a)") t, intclass(t), charclass(t)
     clsize(intclass(t))=clsize(intclass(t))+1
  enddo

!  ! __________________________________________________________________________
!  ! GENERATE KEY FOR TYPE NUMBERS AND STRINGS
!  write(*,"(a,i6)")"number of classes =",nclass
!  write(*,"(a)")"First two letters denote main wind sector number counting clockwise:"
!  write(*,"(a)")"sector 01=NE, 02=SE, 03=SW, 04=NW, 00=undefined (varying directions)"
!  write(*,"(a)")"Third letter denotes >A<nticyclonicity or >C<yclonicity at z925. "
!  write(*,"(a)")"4th letter denotes >A<nticyclonicity or >C<yclonicity at z500. "
!  write(*,"(a)")"5th letter denotes >D<ry or >W<et."
!
!  write(*,*)"nclass =",nclass
!  classnumber=0
!
!  do c925=0,1
!     if(c925==0)then
!        classstring(3:3)="A"
!     else
!        classstring(3:3)="C"
!     endif
!     do c500=0,1
!        if(c500==0)then
!           classstring(4:4)="A"
!        else
!           classstring(4:4)="C"
!        endif
!        do hum=0,1
!           if(hum==0)then
!              classstring(5:5)="D"
!           else
!              classstring(5:5)="W"
!           endif
!           do i=0,nsector
!              write(classstring(1:2),"(1i2.2)")i
!              classnumber=classnumber+1
!              write(*,"(1i6,1a10,1i10)")classnumber,classstring,clsize(classnumber)
!           enddo
!        enddo
!     enddo
!  enddo


end subroutine wlk733

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cycindex(phi,lat,ni,lon,nj,weight, cycindx)
!!!inputs!!! call cycindex(z925(:,:,t),lat,ni,lon,nj,weight, cyc925(t))
  ! calculate mean curvature-index for a field
  ! andreas.philipp@geo.uni-augsburg.de
  implicit none
  integer :: ni, nj, i, j
  real(kind=8) :: phi(ni,nj),lat(ni),lon(nj)
  real(kind=8) :: weight(ni,nj)
  real(kind=8) :: distance
  real(kind=8) :: nabla2(ni,nj),deltax,deltay
  real(kind=8) :: cycindx

  ! cyclonicity index calculated with individual distances
  do i=2,ni-1
     do j=2,nj-1
        deltay = distance(lon(j),lat(i-1),lon(j),lat(i+1))
        deltax = distance(lon(j-1),lat(i),lon(j+1),lat(i))
        nabla2(i,j) = ( phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1)-4.D0*phi(i,j) ) / ( deltax*deltay )
     enddo
  enddo

  ! western and eastern edge
  do i=2,ni-1
     nabla2(i,1)  = nabla2(i,2)
     nabla2(i,nj) = nabla2(i,nj-1)
  enddo

  ! southern and northern edge
  do j=2,nj-1
     nabla2(1,j)  = nabla2(2,j)
     nabla2(ni,j) = nabla2(ni-1,j)
  enddo

  ! corners
  nabla2(1,1)   = nabla2(2,2)
  nabla2(1,nj)  = nabla2(2,nj-1)
  nabla2(ni,1)  = nabla2(ni-1,2)
  nabla2(ni,nj) = nabla2(ni-1,nj-1)

  ! weighted mean
  cycindx=0.D0
  do i=1,ni
     do j=1,nj
        cycindx=cycindx+nabla2(i,j)*weight(i,j)
     enddo
  enddo
  !write(*,*) ">>> ", cycindx, SUM(weight)
  cycindx=cycindx / sum(weight) * 10000.

end subroutine cycindex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=8) function distance(lon1,lat1,lon2,lat2)
  ! return distance between two points on sphere
  ! andreas.philipp@geo.uni-augsburg.de
  implicit none
  real(kind=8) :: lon1,lat1,lon2,lat2,rlon1,rlat1,rlon2,rlat2
  real(kind=8),parameter :: pi=3.141592653589793,radius=6371.0087714
  real(kind=8) :: rad,deg,angle
  rad=pi/180.D0
  deg=180.D0/pi
  rlon1=lon1*rad
  rlat1=lat1*rad
  rlon2=lon2*rad
  rlat2=lat2*rad
  angle=sin(rlat1)*sin(rlat2) + cos(rlat1)*cos(rlat2)*cos(rlon1-rlon2)
  angle=acos(angle)*deg
  distance=(radius * angle * pi) / 180.D0
end function distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine windsectorZAMG(u,v,ni,nj,weight,maindirthreshold,mainsector,verbose)

   ! calculate main wind direction for a grid of u and v
   ! andreas.philipp@geo.uni-augsburg.de, modified version Reto Stauffer
   implicit none
   integer :: ni, nj, i, j, k
   real(kind=8) :: maindirthreshold
   real(kind=8) :: u(ni,nj),v(ni,nj)
   real(kind=8) :: weight(ni,nj)
   real(kind=8) :: winddir
   real(kind=8) :: speed(ni,nj),direction(ni,nj)
   integer      :: sector(ni,nj),mainsector,idx,verbose
   real(kind=8) :: sectorwidth
   real(kind=8) :: maxfreq
 
   ! alternative
   real(kind=8) :: a1,a2
   integer s,mainsec90
   real(kind=8), dimension(8) :: freq ! count frequencies.
                                      ! freq(1) = N  (class 1) 345-015,
                                      ! freq(2) = NE (class 2) 015-075,
                                      ! freq(3) = E  (class 3) 075-105,
                                      ! freq(4) = SE (class 4) 105-165,
                                      ! freq(5) = S  (class 5) 165-195,
                                      ! freq(6) = SW (class 6) 195-255,
                                      ! freq(7) = W  (class 7) 255-285,
                                      ! freq(8) = NW (class 8) 285-345,
   real(kind=8) :: maindir90
   real(kind=8) :: largestangle
   logical :: sectorshift
 
   ! Compute wind direction in degrees
   do i=1,ni
      do j=1,nj
         direction(i,j)=winddir(u(i,j),v(i,j))
      enddo
   enddo

   ! As suggested by Thomas Kennert (Klien, Endbericht Trendanalyse)
   ! for the ZAMG Alpine classification: cont wind directions slightly
   ! different.
   freq = 0.D0
   do i = 1,ni
      do j = 1,nj
         ! 345 - 015 degrees: class 1
         if      ( direction(i,j) .lt.  15.D0 ) then
            idx = 1
         ! 015 - 075 degrees: class 2
         else if ( direction(i,j) .lt.  75.D0 ) then
            idx = 2
         ! 075 - 105 degrees: class 3
         else if ( direction(i,j) .lt. 105.D0 ) then
            idx = 3
         ! 105 - 165 degrees: class 4
         else if ( direction(i,j) .lt. 165.D0 ) then
            idx = 4
         ! 165 - 195 degrees, class 5
         else if ( direction(i,j) .lt. 195.D0 ) then
            idx = 5
         ! 195 - 255 degrees, class 6
         else if ( direction(i,j) .lt. 255.D0 ) then
            idx = 6
         ! 255 - 285 degrees, class 7
         else if ( direction(i,j) .lt. 285.D0 ) then
            idx = 7
         ! 285 - 345 degrees, class 8
         else if ( direction(i,j) .lt. 345.D0 ) then
            idx = 8
         ! 345 - 015 degrees: class 1 again
         else if ( direction(i,j) .le. 360.D0 ) then
            idx = 1
         ! If not classified, stop
         else 
            write(*,*) "Problems with classification. Wind direction unclassified!"
            write(*,*) idx, direction(i,j)
            stop 9
         end if

         freq(idx) = freq(idx) + weight(i,j)
         !write(*,"(I5,2X,I5,2X,F10.5,2X,F10.5,2X,I5)") y, x, weight(x,y), direction(x,y), idx
            
      end do
   end do

   ! Verbose output if required
   if ( verbose .ge. 2 ) then
     write(*,"(3X,A)") "Wind sector classification"
     do k=1,8
        write(*,"(6X,A,I2,2X,F12.1,5X,A,F5.1,A,2X,F5.1,A)") "Sector", k, freq(k), &
              "(sector angle:", sectorwidth*(k-1), "--", sectorwidth*k, ")"
     enddo
     write(*,"(6X,A,2X,F8.1)") "Sum ", SUM(freq)
     write(*,"(6X,A,2X,F8.1)") "Max ", maxfreq
     write(*,"(6X,A,2X,I8)")   "Sec ", mainsector
     write(*,"(6X,A,2X,F8.6)") "Rate", maxfreq / sum(weight)
     write(*,"(6X,A,F5.3)")    "maindirthrehsold: ", maindirthreshold
   endif

   if ( SUM(freq) > SUM(weight) ) then
      write(*,"(AF10.2AF10.2A)") "ERROR: sum of freq (",SUM(freq),") is bigger than sum of weights (",SUM(weight),")"
      stop 3
   endif

   mainsector = -99 ! Default, undefined value
   if ( MAXVAL(freq)/SUM(weight) .lt. maindirthreshold ) then
      mainsector = 0 ! variable,  no distinct main wind direction
   else
      do i=1,8
         if ( freq(i) == MAXVAL(freq) ) mainsector = i
      enddo
   endif
   !write(*,*) mainsector
 

end subroutine windsectorZAMG


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine windsector(u,v,ni,nj,nsector,weight,maindirthreshold,mainsector,verbose)
  ! calculate main wind direction for a grid of u and v
  ! andreas.philipp@geo.uni-augsburg.de
  implicit none
  integer      :: ni, nj, i, j, k, s
  real(kind=8) :: maindirthreshold
  real(kind=8) :: u(ni,nj),v(ni,nj)
  real(kind=8) :: weight(ni,nj)
  real(kind=8) :: winddir
  real(kind=8) :: speed(ni,nj),direction(ni,nj)
  integer      :: nsector,sector(ni,nj),mainsector,verbose
  real(kind=8) :: sectorwidth
  real(kind=8) :: freq(nsector),maxfreq

  sectorwidth=360.D0/nsector

  do i=1,ni
     do j=1,nj
        direction(i,j)=winddir(u(i,j),v(i,j))
     enddo
  enddo

  ! JUST COUNT FOR nsector SECTORS
  ! determine sectors
  do i=1,ni
     do j=1,nj
        do s=1,nsector+1
           !write(*,"(i4,4f10.4)")i,direction(x,y),sectorwidth*(i-1)
           sector(i,j)=s-1
           if(direction(i,j) < sectorwidth*(s-1))exit
        enddo
        !write(*,*)sector(x,y)
     enddo
  enddo

  ! determine frequency for sectors
  freq=0.D0
  do i=1,ni
     do j=1,nj
        freq(sector(i,j))=freq(sector(i,j))+weight(i,j)
     enddo
  enddo
  maxfreq=0.D0
  do s=1,nsector
     if ( freq(s) .gt. maxfreq ) then
        maxfreq    = freq(s)
        mainsector = s
     endif
  enddo

  ! Verbose output if required
  if ( verbose .ge. 2 ) then
    write(*,"(3X,A)") "Wind sector classification"
    write(*,"(3X,A,I3)") "Sectors used:",nsector
    do k=1,8
       write(*,"(6X,A,I2,2X,F12.1,5X,A,F5.1,A,2X,F5.1,A)") "Sector", k, freq(k), &
             "(sector angle:", sectorwidth*(k-1), "--", sectorwidth*k, ")"
    enddo
    write(*,"(6X,A,2X,F8.1)") "Sum ", SUM(freq)
    write(*,"(6X,A,2X,F8.1)") "Max ", maxfreq
    write(*,"(6X,A,2X,I8)")   "Sec ", mainsector
    write(*,"(6X,A,2X,F8.6)") "Rate", maxfreq / sum(weight)
    write(*,"(6X,A,F5.3)")    "maindirthrehsold: ", maindirthreshold
  endif

  ! if maxfreq < 66.7(maindirthreshold) then main winddirection is undefined
  if( maxfreq/sum(weight) .lt. maindirthreshold ) then
     if ( verbose .gt. 0 ) write(*,"(3X,A)") "Relative frequency lower than maindirthreshold: set wind direction to 0"
     mainsector=0
     return
  endif

end subroutine windsector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=8) function winddir(u,v)
  ! returns direction of wind origin in grad
  ! andreas.philipp@geo.uni-augsburg.de
  implicit none
  real(kind=8) :: u,v
  real(kind=8),parameter :: pi=3.14159265358979323846
  ! zero zonal wind
  if(u==0.D0)then
     if(v==0.D0)then
        winddir=-999.D0
        return
     endif
     if(v>0.D0)then
        winddir=180.D0
        return
     else
        winddir=0.D0
        return
     endif
  endif
  ! zero meridional wind
  if(v==0.D0)then
     if(u>0.D0)then
        winddir=270.D0
        return
     else
        winddir=90.D0
        return
     endif
  endif
  ! rest
  if(u>0.D0)then
     if(v>0.D0)then
        winddir=270.D0-atan(abs(v/u))*180.D0/pi
        return
     else
        winddir=360.D0-atan(abs(u/v))*180.D0/pi
        return
     endif
  else
     if(v>0.D0)then
        winddir=180.D0-atan(abs(u/v))*180.D0/pi
        return
     else
        winddir=90.D0-atan(abs(v/u))*180.D0/pi
        return
     endif
  endif
end function winddir

