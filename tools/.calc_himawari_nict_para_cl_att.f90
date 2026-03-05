program test
! ひまわり 8 号標準データを read_hsd.rb で NetCDF に変換した
! データ処理を示唆補正するプログラム.
! 水平はピクセル座標で格納されており, 緯度経度も 2 次元でデータに含まれる.
! IR データから STPK のルーチンを用いて, 視差補正を行う.
! 補正時に必要となる雲頂輝度温度と高度の対応関係は,
! 任意のサウンディングデータ (テキスト) から鉛直方向線形内挿で作成する.
! サウンディングデータのリストは時刻情報と合わせて 1 つのテキストファイルに
! 収める. 
! calc_himawari_para_cl と同じネームリストファイルで実行する.
! 堀之内研 hsv 視差補正用プログラム (hsv でのデータ変換ではこれを使用). 

! サウンディングデータリストファイルのフォーマット
! 1 行目: サウンディングデータ
! 2 行目: 時刻 (YYYYMMDDHHNNSS)

! 各サウンディングデータの高度座標は全時刻で同じ高度座標値を
! 持っていると仮定している. 
! 再解析を用いる場合は指定気圧面でデータを作ればよい.
! 各サウンディングデータのファイルフォーマットは sndcord で指定.
! 't' = 温度データの入っているカラム.
! 'z' = 高度データの入っているカラム.
! 'x' = 読み飛ばしカラム.

  use gtool_history
  use Math_Const
  use file_operate
  use Basis
  use statistics
  use map_function
  use max_min
  use Typhoon_Analy


  implicit none

  character(100), parameter :: colname='column    '
  character(100), parameter :: linname='line      '
  character(100), parameter :: lonname='longitude '
  character(100), parameter :: latname='latitude  '
  character(100), parameter :: tsname='start_time '
  character(100), parameter :: tename='end_time   '
  character(100), parameter :: tname='obs_time   '
  character(100), parameter :: lonlname='longitude             '
  character(100), parameter :: latlname='latitude              '
  character(100), parameter :: tslname='observation start time '
  character(100), parameter :: telname='observation end time   '
  character(100), parameter :: tlname='half time of start-end '
  character(100), parameter :: londname='degrees_east               '
  character(100), parameter :: latdname='degrees_north              '
  character(100), parameter :: tsdname='days since 1858-11-17 0:0:0'
  character(100), parameter :: tedname='days since 1858-11-17 0:0:0'
  character(100), parameter :: tdname='days since 1858-11-17 0:0:0'
  character(100), parameter :: undefname='_FillValue'

  character(100), dimension(12), parameter :: addattr=  &
  &               (/'satellite_name                 ',  &
  &                 'column_number                  ',  &
  &                 'line_number                    ',  &
  &                 'equator_radius                 ',  &
  &                 'polar_radius                   ',  &
  &                 'satellite_longitude            ',  &
  &                 'satellite_latitude             ',  &
  &                 'satellite_earth_center_distance',  &
  &                 'nadir_longitude                ',  &
  &                 'nadir_latitude                 ',  &
  &                 'wavelength                     ',  &
  &                 'band_number                    '/)

  integer, parameter :: axnum=2
  real, parameter :: outundef=-999.0

  !-- namelist variables
  integer :: sndskip
  real :: rundef, dlon, dlat, dlona, dlata
  double precision :: dundef, dlth
  character(20) :: sndcundef
  character(100) :: sndcord
  character(1000) :: fname, vname, addvname, snd_data, addfname
  character(100), dimension(axnum) :: axname, addaxname
  character(100), dimension(axnum) :: axintername, addaxintername
  character(100), dimension(axnum) :: lonlatname, addlonlatname  ! (/lonname2d, latname2d/)
  logical :: flag_para

  integer :: stat, ncid, varid, aerr, itmp, iattr
  type(dtime) :: stime, etime, rtime, sndtime, sndetime
  integer :: i, j, l, id, itime, nsnd, tempord, zphord, nanstat
  integer :: nt, nts, nfs, nl, nlad, hour, sec
  integer, dimension(axnum) :: dimn, dimon, diman, dimoan
  real :: t1, t2
  real, allocatable, dimension(:) :: lonout, latout, lonouta, latouta
  real, allocatable, dimension(:,:) :: val, tbb, zph
  real, allocatable, dimension(:,:,:) :: tbbc, zphc, tbbca, zphca
  real, allocatable, dimension(:,:) :: vala, tbba, zpha
  real, pointer, dimension(:) :: a1
  double precision :: dattr
  double precision :: tmptime, obstime
  double precision :: londmin, londmax, latdmin, latdmax
  double precision :: londamin, londamax, latdamin, latdamax
  double precision, allocatable, dimension(:) :: timee, times, timesnd
  double precision, allocatable, dimension(:) :: timeea, timesa
  double precision, allocatable, dimension(:) :: temp1d, zph1d
  double precision, allocatable, dimension(:,:) :: sndtemp, sndzph
  double precision, allocatable, dimension(:) :: xinter, yinter, xintera, yintera
  double precision, allocatable, dimension(:) :: lonoutd, latoutd, lonoutda, latoutda
  double precision, allocatable, dimension(:,:) :: lond, latd, londa, latda
  double precision, allocatable, dimension(:,:) :: londr, latdr, londcr, latdcr
  double precision, allocatable, dimension(:,:) :: londra, latdra, londcra, latdcra
  double precision, allocatable, dimension(:,:) :: londc, latdc, londca, latdca
  double precision, allocatable, dimension(:,:) :: tbbd, zphd, tbbcd, zphcd
  double precision, allocatable, dimension(:,:) :: tbbda, zphda, tbbcda, zphcda
  character(1000) :: tmpc, oname, cattr
  character(100) :: vlname, vdname, addvlname, addvdname
  character(100) :: ftitle_sub, fvname_sub, ffoot_sub
  character(100), allocatable, dimension(:,:) :: sndcval
  character(1000), allocatable, dimension(:) :: cval, caval
  character(1000), allocatable, dimension(:,:) :: ccval
  logical :: add_flag

  integer, external :: NF_OPEN, NF_INQ_VARID, NF_GET_VAR1_DOUBLE, NF_CLOSE
  integer, external :: NF_INQ_DIMID, NF_INQ_DIMLEN
  type(GT_HISTORY) :: vhst, vhsta

  namelist /basic /fname, vname, axname, lonlatname, dlon, dlat,  &
  &                addfname, addvname, addaxname, addlonlatname, dlona, dlata,  &
  &                axintername, addaxintername
  namelist /option /rundef, dundef, snd_data, sndcord, sndskip, sndcundef, dlth, flag_para
  read(5,nml=basic)
  read(5,nml=option)

  dlth=dlth*pi_dp/180.0d0

  if(trim(adjustl(addfname))/='')then
     add_flag=.true.
  else
     add_flag=.false.
  end if

  if(flag_para.eqv..true.)then
     ftitle_sub="Parallax correction"
     fvname_sub="Parallax corrected"
     ffoot_sub=".para"
  else
     write(*,*) "*** MESSAGE (main) ***: This mode is not performing parallax collection."
     ftitle_sub="Pixcel to lat-lon"
     fvname_sub="Pix to LL interpolated"
     ffoot_sub=".cl2ll"
  end if

  vlname='brightness_temperature'
  vdname='K'

  stime%year_d=1858
  stime%month_d=11
  stime%day_d=17
  stime%hour_d=0
  stime%min_d=0
  stime%sec_d=0

  nsnd=len_trim(adjustl(sndcord))
  zphord=0
  tempord=0
  do i=1,nsnd
     select case (sndcord(i:i))
     case ('z')
        zphord=i
     case ('t')
        tempord=i
     end select
  end do

  nl=line_number_counter( trim(adjustl(fname)) )
  allocate(cval(nl))
  call read_file_text( trim(adjustl(fname)), 1, nl, cval )

  if(add_flag.eqv..true.)then
     nlad=line_number_counter( trim(adjustl(fname)) )
     if(nl/=nlad)then
        write(*,*) "*** ERROR (main) ***: nl is not identical to nlad."
        stop
     end if
     allocate(caval(nl))
     call read_file_text( trim(adjustl(addfname)), 1, nl, caval )
  else
     nlad=0
  end if

  nts=line_number_counter( trim(adjustl(snd_data)) )
  allocate(ccval(2,nts))
  allocate(timesnd(nts))
  call read_file_text( trim(adjustl(snd_data)), 2, nts, ccval )

  tmpc=trim(adjustl(ccval(2,1)))
  sndtime%year_d=c2i_convert( tmpc(1:4) )
  sndtime%month_d=c2i_convert( tmpc(5:6) )
  sndtime%day_d=c2i_convert( tmpc(7:8) )
  sndtime%hour_d=c2i_convert( tmpc(9:10) )
  sndtime%min_d=c2i_convert( tmpc(11:12) )
  sndtime%sec_d=c2i_convert( tmpc(13:14) )
  timesnd(1)=0.0d0

  nfs=line_number_counter( trim(adjustl(ccval(1,1))) )-sndskip
  allocate(sndcval(nsnd,nfs))
  allocate(sndtemp(nfs,nts))
  allocate(sndzph(nfs,nts))
  allocate(temp1d(nfs))
  allocate(zph1d(nfs))
  sndtemp=dundef
  sndzph=dundef
  call read_file_text( trim(adjustl(ccval(1,1))), nsnd, nfs, sndcval, skip=sndskip )
  do i=1,nfs
     if((trim(adjustl(sndcval(tempord,i)))/=trim(adjustl(sndcundef))).and.  &
  &     (trim(adjustl(sndcval(zphord,i)))/=trim(adjustl(sndcundef))))then
        sndtemp(i,1)=dble(c2r_convert( trim(adjustl(sndcval(tempord,i))) ))
        sndzph(i,1)=dble(c2r_convert( trim(adjustl(sndcval(zphord,i))) ))
     end if
  end do

  do i=2,nts
     tmpc=trim(adjustl(ccval(2,i)))
     sndetime%year_d=c2i_convert( tmpc(1:4) )
     sndetime%month_d=c2i_convert( tmpc(5:6) )
     sndetime%day_d=c2i_convert( tmpc(7:8) )
     sndetime%hour_d=c2i_convert( tmpc(9:10) )
     sndetime%min_d=c2i_convert( tmpc(11:12) )
     sndetime%sec_d=c2i_convert( tmpc(13:14) )
     timesnd(i)=dble( counter_sec( sndtime, sndetime ) )
     call read_file_text( trim(adjustl(ccval(1,i))), nsnd, nfs, sndcval, skip=sndskip )
     do j=1,nfs
        if(trim(adjustl(sndcval(tempord,j)))/=trim(adjustl(sndcundef)).and.  &
  &        trim(adjustl(sndcval(zphord,j)))/=trim(adjustl(sndcundef)))then
           sndtemp(j,i)=dble(c2r_convert( trim(adjustl(sndcval(tempord,j))) ))
           sndzph(j,i)=dble(c2r_convert( trim(adjustl(sndcval(zphord,j))) ))
        end if
     end do
  end do

  do l=1,nl

     do i=1,axnum
        nullify(a1)
        call HistoryGetPointer( trim(adjustl(cval(l))),  &
  &                             trim(adjustl(axname(i))), a1 )
        dimn(i)=size(a1)
        deallocate(a1)
     end do

     if(add_flag.eqv..true.)then
        do i=1,axnum
           nullify(a1)
           call HistoryGetPointer( trim(adjustl(caval(l))),  &
  &                                trim(adjustl(addaxname(i))), a1 )
           diman(i)=size(a1)
           deallocate(a1)
        end do
     end if

!    変数の要素数も後ろの NF_INQ_DIMID, NF_INQ_DIMLEN で見るようにした.
!     nullify(a1)
!     call HistoryGetPointer( trim(cval(l)), trim(tename), a1 )
!     nt=size(a1)
!     deallocate(a1)

     allocate(xinter(dimn(1)),stat=aerr)
     allocate(yinter(dimn(2)),stat=aerr)
     allocate(lond(dimn(1),dimn(2)),stat=aerr)
     allocate(latd(dimn(1),dimn(2)),stat=aerr)
     allocate(londr(dimn(1),dimn(2)),stat=aerr)
     allocate(latdr(dimn(1),dimn(2)),stat=aerr)
     allocate(londc(dimn(1),dimn(2)),stat=aerr)
     allocate(latdc(dimn(1),dimn(2)),stat=aerr)
     allocate(londcr(dimn(1),dimn(2)),stat=aerr)
     allocate(latdcr(dimn(1),dimn(2)),stat=aerr)
     allocate(tbbd(dimn(1),dimn(2)),stat=aerr)
     allocate(zphd(dimn(1),dimn(2)),stat=aerr)
     allocate(tbb(dimn(1),dimn(2)),stat=aerr)
     allocate(zph(dimn(1),dimn(2)),stat=aerr)
     allocate(val(dimn(1),dimn(2)),stat=aerr)

     if(add_flag.eqv..true.)then
        allocate(xintera(diman(1)),stat=aerr)
        allocate(yintera(diman(2)),stat=aerr)
        allocate(londa(diman(1),diman(2)),stat=aerr)
        allocate(latda(diman(1),diman(2)),stat=aerr)
        allocate(londra(diman(1),diman(2)),stat=aerr)
        allocate(latdra(diman(1),diman(2)),stat=aerr)
        allocate(londca(diman(1),diman(2)),stat=aerr)
        allocate(latdca(diman(1),diman(2)),stat=aerr)
        allocate(londcra(diman(1),diman(2)),stat=aerr)
        allocate(latdcra(diman(1),diman(2)),stat=aerr)
        allocate(tbbda(diman(1),diman(2)),stat=aerr)
        allocate(zphda(diman(1),diman(2)),stat=aerr)
        allocate(tbba(diman(1),diman(2)),stat=aerr)
        allocate(zpha(diman(1),diman(2)),stat=aerr)
        allocate(vala(diman(1),diman(2)),stat=aerr)
     end if
!     if(aerr/=0)then
!        write(*,*) "Detect allocating error; stop."
!        stop
!     end if

     call HistoryGet( trim(cval(l)), trim(axintername(1)), xinter )
     call HistoryGet( trim(cval(l)), trim(axintername(2)), yinter )
     call HistoryGet( trim(cval(l)), trim(lonlatname(1)), lond )
     call HistoryGet( trim(cval(l)), trim(lonlatname(2)), latd )
     nanstat=check_NaN( lond, latd )
     !nanstat=check_NaN( lond, latd, rep_val=dundef )

     if(nanstat==0)then
        call min_val_2d( lond, itmp, itmp, londmin, undef=dundef )
        call max_val_2d( lond, itmp, itmp, londmax, undef=dundef )
        call min_val_2d( latd, itmp, itmp, latdmin, undef=dundef )
        call max_val_2d( latd, itmp, itmp, latdmax, undef=dundef )

        londmin=dint(londmin)  ! 小数点以下切り捨て
        latdmin=dint(latdmin)  ! 小数点以下切り捨て
        londmax=dint(londmax)+1.0d0  ! 小数点以下切り捨て
        latdmax=dint(latdmax)+1.0d0  ! 小数点以下切り捨て

        dimon(1)=idint((londmax-londmin)/dble(dlon))+1
        dimon(2)=idint((latdmax-latdmin)/dble(dlat))+1

        allocate(lonout(dimon(1)),stat=aerr)
        allocate(latout(dimon(2)),stat=aerr)
        allocate(lonoutd(dimon(1)),stat=aerr)
        allocate(latoutd(dimon(2)),stat=aerr)
        allocate(tbbcd(dimon(1),dimon(2)),stat=aerr)
        allocate(zphcd(dimon(1),dimon(2)),stat=aerr)
        allocate(tbbc(dimon(1),dimon(2),1),stat=aerr)
        allocate(zphc(dimon(1),dimon(2),1),stat=aerr)

        lonout=(/((real(londmin)+real(i-1)*dlon),i=1,dimon(1))/)
        latout=(/((real(latdmin)+real(i-1)*dlat),i=1,dimon(2))/)
        lonoutd=(/((londmin+dble(i-1)*dble(dlon)),i=1,dimon(1))/)
        latoutd=(/((latdmin+dble(i-1)*dble(dlat)),i=1,dimon(2))/)

        if(add_flag.eqv..true.)then
           call HistoryGet( trim(caval(l)), trim(addaxintername(1)), xintera )
           call HistoryGet( trim(caval(l)), trim(addaxintername(2)), yintera )
           call HistoryGet( trim(caval(l)), trim(addlonlatname(1)), londa )
           call HistoryGet( trim(caval(l)), trim(addlonlatname(2)), latda )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             'long_name', addvlname )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             'units', addvdname )

           nanstat=check_NaN( londa, latda )
           !nanstat=check_NaN( londa, latda, rep_val=dundef )  ! nanstat is dummy
           call min_val_2d( londa, itmp, itmp, londamin, undef=dundef )
           call max_val_2d( londa, itmp, itmp, londamax, undef=dundef )
           call min_val_2d( latda, itmp, itmp, latdamin, undef=dundef )
           call max_val_2d( latda, itmp, itmp, latdamax, undef=dundef )

           londamin=dint(londamin)  ! 小数点以下切り捨て
           latdamin=dint(latdamin)  ! 小数点以下切り捨て
           londamax=dint(londamax)+1.0d0  ! 小数点以下切り捨て
           latdamax=dint(latdamax)+1.0d0  ! 小数点以下切り捨て

           dimoan(1)=idint((londamax-londamin)/dble(dlona))+1
           dimoan(2)=idint((latdamax-latdamin)/dble(dlata))+1

           allocate(lonouta(dimoan(1)),stat=aerr)
           allocate(latouta(dimoan(2)),stat=aerr)
           allocate(lonoutda(dimoan(1)),stat=aerr)
           allocate(latoutda(dimoan(2)),stat=aerr)
           allocate(tbbcda(dimoan(1),dimoan(2)),stat=aerr)
           allocate(zphcda(dimoan(1),dimoan(2)),stat=aerr)
           allocate(tbbca(dimoan(1),dimoan(2),1),stat=aerr)
           allocate(zphca(dimoan(1),dimoan(2),1),stat=aerr)

           lonouta=(/((real(londamin)+real(i-1)*dlona),i=1,dimoan(1))/)
           latouta=(/((real(latdamin)+real(i-1)*dlata),i=1,dimoan(2))/)
           lonoutda=(/((londamin+dble(i-1)*dble(dlona)),i=1,dimoan(1))/)
           latoutda=(/((latdamin+dble(i-1)*dble(dlata)),i=1,dimoan(2))/)
        end if
!        call HistoryGet( trim(cval(l)), trim(tname), time(1) )
!--------- NetCDF original library ---------------
! start_time 変数のみ, gtool5 では正常に読み込めないので,
! netcdf のプリミティブ関数を用いて読み込む.
        stat=NF_OPEN( trim(adjustl(cval(l))), 0, ncid )
        stat=NF_INQ_DIMID( ncid, trim(tsname), varid )
        stat=NF_INQ_DIMLEN( ncid, varid, nt )
write(*,*) "checkid", nt
        stat=NF_INQ_VARID( ncid, trim(tsname), varid )
        allocate(times(nt),stat=aerr)
        stat=NF_GET_VAR1_DOUBLE( ncid, varid, 1, times(1:nt) )
        stat=NF_CLOSE( ncid )

        stat=NF_OPEN( trim(adjustl(cval(l))), 0, ncid )
        stat=NF_INQ_DIMID( ncid, trim(tename), varid )
        stat=NF_INQ_DIMLEN( ncid, varid, nt )
write(*,*) "checkid", nt
        stat=NF_INQ_VARID( ncid, trim(tename), varid )
        allocate(timee(nt),stat=aerr)
        stat=NF_GET_VAR1_DOUBLE( ncid, varid, 1, timee(1:nt) )
        stat=NF_CLOSE( ncid )
        if(add_flag.eqv..true.)then
           stat=NF_OPEN( trim(adjustl(caval(l))), 0, ncid )
           stat=NF_INQ_DIMID( ncid, trim(tsname), varid )
           stat=NF_INQ_DIMLEN( ncid, varid, nt )
write(*,*) "checkid", nt
           stat=NF_INQ_VARID( ncid, trim(tsname), varid )
           allocate(timesa(nt),stat=aerr)
           stat=NF_GET_VAR1_DOUBLE( ncid, varid, 1, timesa(1:nt) )
           stat=NF_CLOSE( ncid )

           stat=NF_OPEN( trim(adjustl(caval(l))), 0, ncid )
           stat=NF_INQ_DIMID( ncid, trim(tename), varid )
           stat=NF_INQ_DIMLEN( ncid, varid, nt )
write(*,*) "checkid", nt
           stat=NF_INQ_VARID( ncid, trim(tename), varid )
           allocate(timeea(nt),stat=aerr)
           stat=NF_GET_VAR1_DOUBLE( ncid, varid, 1, timeea(1:nt) )
           stat=NF_CLOSE( ncid )
        end if
!--------- NetCDF original library ---------------
        call HistoryGet( trim(cval(l)), trim(vname), val )
        if(add_flag.eqv..true.)then
           call HistoryGet( trim(caval(l)), trim(addvname), vala )
        end if

        do id=1,nt
           obstime = 0.5d0*(times(id)+timee(id))
           hour=idint(obstime*24.0d0)
write(*,*) "time checkkkk", times(id), timee(id), id, l, idint(obstime*24.0d0),  &
        int(real(obstime*24.0d0)), obstime*24.0d0
           call time_zone_convert( hour, stime, etime )
           sec=idint((obstime*24.0d0-dble(hour))*3600.0d0)
           call sec_convert( sec, etime, rtime )
write(*,*) "time check", times(id), timee(id), obstime, hour, sec, stime, etime, rtime
        end do

        tmptime=dble( counter_sec( sndtime, rtime ) )

        temp1d=dundef
        zph1d=dundef
        call interpo_search_1d( timesnd, tmptime, itime, stdopt=.true. )
        do i=1,nfs
           if(sndtemp(i,itime)/=dundef.and.sndtemp(i,itime+1)/=dundef.and.  &
  &           sndzph(i,itime)/=dundef.and.sndzph(i,itime+1)/=dundef)then
              call interpolation_1d( timesnd(itime:itime+1), sndtemp(i,itime:itime+1),  &
  &                                  tmptime, temp1d(i) )
              call interpolation_1d( timesnd(itime:itime+1), sndzph(i,itime:itime+1),  &
  &                                  tmptime, zph1d(i) )
           end if
        end do

        call conv_unit_2d( lond, latd, londr, latdr, pi_dp/180.0d0,  &
  &                        dundef )
        if(add_flag.eqv..true.)then
           call conv_unit_2d( londa, latda, londra, latdra, pi_dp/180.0d0,  &
  &                           dundef )
        end if

        call conv_r4_r8( val, tbbd, rundef=rundef, dundef=dundef )

        call convert_Tbb2Zph( tbbd, zphd, temp1d, zph1d, undef=dundef )

        if(flag_para.eqv..true.)then
           call Parallax_Himawari( londr, latdr, zphd, londcr, latdcr, undef=dundef )
        else
           londcr=londr
           latdcr=latdr
        end if

        call conv_unit_2d( londcr, latdcr, londc, latdc, 180.0d0/pi_dp, dundef )

!        call cpu_time(t1)
        call tri_interpolation_2d( londc, latdc, zphd, tbbd,  &
  &                                lonoutd(1:dimon(1)), latoutd(1:dimon(2)),  &
  &                                zphcd, tbbcd, undef=dundef, jflag='u' )
!        call tri_interpolation_2d( londcr, latdcr, tbbd,  &
!  &                                londr(1:dimn(1),1), latdr(1,1:dimn(2)),  &
!  &                                tbbcd, undef=dundef, jflag='l' )
!        call cpu_time(t2)
        if(add_flag.eqv..true.)then
           call conv_r4_r8( vala, tbbda, rundef=rundef, dundef=dundef )

           call auto_interpolation_2d( xinter(1:dimn(1)),  &
  &                                    yinter(1:dimn(2)),  &
  &                                    xintera(1:diman(1)),  &
  &                                    yintera(1:diman(2)),  &
  &                                    zphd, zphda, undef=dundef, undefr=dundef )

           if(flag_para.eqv..true.)then
              call Parallax_Himawari( londra, latdra, zphda, londcra, latdcra,  &
  &                                   undef=dundef )
           else
              londcra=londra
              latdcra=latdra
           end if

           call conv_unit_2d( londcra, latdcra, londca, latdca, 180.0d0/pi_dp,  &
  &                           dundef )

           call tri_interpolation_2d( londca, latdca, zphda, tbbda,  &
  &                                   lonoutda(1:dimoan(1)), latoutda(1:dimoan(2)),  &
  &                                   zphcda, tbbcda, undef=dundef, jflag='u' )
!           call tri_interpolation_2d( londcra, latdcra, tbbda,  &
!  &                                   londra(1:diman(1),1),  &
!  &                                   latdra(1,1:diman(2)),  &
!  &                                   tbbcda, undef=dundef )
        end if
!        write(*,*) "check 1", t2-t1
!        call nearest_neighbor_2d( londcr, latdcr, zphd,  &
!  &                               londr(1:dimn(1),1), latdr(1,1:dimn(2)), zphcd, dundef, dlth )
!        call nearest_neighbor_2d( londcr, latdcr, tbbd,  &
!  &                               londr(1:dimn(1),1), latdr(1,1:dimn(2)), tbbcd, dundef, dlth )
!        call cpu_time(t1)
!        write(*,*) "check 2", t1-t2
!        stop

        call conv_r8_r4( tbbcd, tbbc(1:dimon(1),1:dimon(2),1),  &
  &                      rundef=outundef, dundef=dundef )
        call conv_r8_r4( zphcd, zphc(1:dimon(1),1:dimon(2),1),  &
  &                      rundef=outundef, dundef=dundef )
        if(add_flag.eqv..true.)then
           call conv_r8_r4( tbbcda, tbbca(1:dimoan(1),1:dimoan(2),1),  &
  &                         rundef=outundef, dundef=dundef )
           call conv_r8_r4( zphcda, zphca(1:dimoan(1),1:dimoan(2),1),  &
  &                         rundef=outundef, dundef=dundef )
        end if

  !-- NetCDF ファイルの定義
        oname=trim(adjustl(cval(l)))//trim(adjustl(ffoot_sub))//'.nc'
        call HistoryCreate( &                        ! ヒストリー作成
  &          file=trim(adjustl(oname)),  &
  &          title=trim(adjustl(ftitle_sub)), &
  &          source='test',   &
  &          institution='',       &
  &          dims=(/lonname, latname, tname/),  &
  &          dimsizes=(/dimon(1),dimon(2),1/),  &
  &          longnames=(/lonlname, latlname, tlname/),  &
  &          units=(/londname, latdname, tdname/), history=vhst,  &
  &          xtypes=(/'float ', 'float ', 'double'/) )

        call HistoryAddVariable( &                   ! 変数定義
  &          varname=trim(adjustl(vname)),  &
  &          dims=(/lonname, latname, tname/),  &
  &          longname=vlname,  &
  &          units=vdname, xtype='float', history=vhst )

        call HistoryAddVariable( &                   ! 変数定義
  &          varname='gph',  &
  &          dims=(/lonname, latname, tname/),  &
  &          longname=trim(adjustl(fvname_sub))//' geopotential height',  &
  &          units='m', xtype='float', history=vhst )

        call HistoryPut( lonname, lonout, history=vhst )
        call HistoryPut( latname, latout, history=vhst )
        call HistoryPut( tname, 0.5d0*(times(1)+timee(1)), history=vhst )

        if(add_flag.eqv..true.)then
        !-- NetCDF ファイルの定義
           oname=trim(adjustl(caval(l)))//trim(adjustl(ffoot_sub))//'.nc'
           call HistoryCreate( &                        ! ヒストリー作成
  &             file=trim(adjustl(oname)),  &
  &             title=trim(adjustl(ftitle_sub)), &
  &             source='test',   &
  &             institution='',       &
  &             dims=(/lonname, latname, tname/),  &
  &             dimsizes=(/dimoan(1),dimoan(2),1/),  &
  &             longnames=(/lonlname, latlname, tlname/),  &
  &             units=(/londname, latdname, tdname/), history=vhsta,  &
  &             xtypes=(/'float ', 'float ', 'double'/) )

           call HistoryAddVariable( &                   ! 変数定義
  &             varname=trim(adjustl(addvname)),  &
  &             dims=(/lonname, latname, tname/),  &
  &             longname=addvlname,  &
  &             units=addvdname, xtype='float', history=vhsta )

           call HistoryAddVariable( &                   ! 変数定義
  &             varname='gph',  &
  &             dims=(/lonname, latname, tname/),  &
  &             longname=trim(adjustl(fvname_sub))//' geopotential height',  &
  &             units='m', xtype='float', history=vhsta )

           call HistoryPut( lonname, lonouta, history=vhsta )
           call HistoryPut( latname, latouta, history=vhsta )
           call HistoryPut( tname, 0.5d0*(timesa(1)+timeea(1)), history=vhsta )
        end if

        call HistoryPut( trim(adjustl(vname)), tbbc, history=vhst )

        call HistoryPut( 'gph', zphc, history=vhst )

        call HistoryAddAttr( trim(adjustl(vname)), trim(adjustl(undefname)),  &
  &                          outundef, history=vhst )
        call HistoryAddAttr( 'gph', trim(adjustl(undefname)),  &
  &                          outundef, history=vhst )

        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(1))), cattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(1))),  &
  &                          trim(adjustl(cattr)), history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(2))), iattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(2))),  &
  &                          iattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(3))), iattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(3))),  &
  &                          iattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(4))), dattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(4))),  &
  &                          dattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(5))), dattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(5))),  &
  &                          dattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(6))), dattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(6))),  &
  &                          dattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(7))), dattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(7))),  &
  &                          dattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(8))), dattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(8))),  &
  &                          dattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(9))), dattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(9))),  &
  &                          dattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(10))), dattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(10))),  &
  &                          dattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(11))), dattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(11))),  &
  &                          dattr, history=vhst )
        call HistoryGetAttr( trim(adjustl(cval(l))), trim(adjustl(vname)),  &
  &                          trim(adjustl(addattr(12))), iattr )
        call HistoryAddAttr( trim(adjustl(vname)), '+'//trim(adjustl(addattr(12))),  &
  &                          iattr, history=vhst )

        call HistoryClose( history=vhst )

        if(add_flag.eqv..true.)then
           call HistoryPut( trim(adjustl(addvname)), tbbca, history=vhsta )

           call HistoryPut( 'gph', zphca, history=vhsta )

           call HistoryAddAttr( trim(adjustl(addvname)), trim(adjustl(undefname)),  &
  &                             outundef, history=vhsta )
           call HistoryAddAttr( 'gph', trim(adjustl(undefname)),  &
  &                             outundef, history=vhsta )

           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(1))), cattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(1))),  &
  &                             trim(adjustl(cattr)), history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(2))), iattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(2))),  &
  &                             iattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(3))), iattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(3))),  &
  &                             iattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(4))), dattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(4))),  &
  &                             dattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(5))), dattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(5))),  &
  &                             dattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(6))), dattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(6))),  &
  &                             dattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(7))), dattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(7))),  &
  &                             dattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(8))), dattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(8))),  &
  &                             dattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(9))), dattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(9))),  &
  &                             dattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(10))), dattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(10))),  &
  &                             dattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(11))), dattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(11))),  &
  &                             dattr, history=vhsta )
           call HistoryGetAttr( trim(adjustl(caval(l))), trim(adjustl(addvname)),  &
  &                             trim(adjustl(addattr(12))), iattr )
           call HistoryAddAttr( trim(adjustl(addvname)), '+'//trim(adjustl(addattr(12))),  &
  &                             iattr, history=vhsta )

           call HistoryClose( history=vhsta )
        end if

        deallocate(lonout)
        deallocate(latout)
        deallocate(lonoutd)
        deallocate(latoutd)
        deallocate(tbbcd)
        deallocate(zphcd)
        deallocate(tbbc)
        deallocate(zphc)
        deallocate(times)
        deallocate(timee)

        if(add_flag.eqv..true.)then
           deallocate(lonouta)
           deallocate(latouta)
           deallocate(lonoutda)
           deallocate(latoutda)
           deallocate(tbbcda)
           deallocate(zphcda)
           deallocate(tbbca)
           deallocate(zphca)
           deallocate(timesa)
           deallocate(timeea)
        end if

     else

        write(*,*) "*** MESSAGE (main) ***: Data is NaN. Skip to make file."

     end if

     deallocate(lond)
     deallocate(latd)
     deallocate(londr)
     deallocate(latdr)
     deallocate(londcr)
     deallocate(latdcr)
     deallocate(tbbd)
     deallocate(zphd)
     deallocate(tbb)
     deallocate(zph)
     deallocate(val)

     if(add_flag.eqv..true.)then
        deallocate(londa)
        deallocate(latda)
        deallocate(londra)
        deallocate(latdra)
        deallocate(londcra)
        deallocate(latdcra)
        deallocate(tbbda)
        deallocate(zphda)
        deallocate(tbba)
        deallocate(zpha)
        deallocate(vala)
     end if

  end do

contains

integer function check_NaN( val1, val2, rep_val )
  implicit none
  double precision, intent(inout) :: val1(:,:)
  double precision, intent(inout) :: val2(size(val1,1),size(val1,2))
  double precision, intent(in), optional :: rep_val
  integer ii, jj, ix, jy, nanst

  ix=size(val1,1)
  jy=size(val1,2)

  nanst=0

  if(present(rep_val))then
     write(*,*) "*** MESSAGE (check_NaN) ***: NaN values are replaced with", rep_val
     do jj=1,jy
        do ii=1,ix
!ORG           if(val1(ii,jj)/=val1(ii,jj))then
           if(val1(ii,jj)*0.0/=0.0)then  ! NaN は任意の演算結果が NaN になるので, 0 を掛けた結果がゼロでなければ NaN.
              val1(ii,jj)=rep_val
           end if
!ORG           if(val2(ii,jj)/=val2(ii,jj))then
           if(val2(ii,jj)*0.0/=0.0)then
              val2(ii,jj)=rep_val
           end if
        end do
     end do
  else
     do jj=1,jy
        do ii=1,ix
!ORG           if(val1(ii,jj)/=val1(ii,jj))then
           if(val1(ii,jj)*0.0/=0.0)then
              nanst=1
              exit
           end if
!ORG           if(val2(ii,jj)/=val2(ii,jj))then
           if(val2(ii,jj)*0.0/=0.0)then
              nanst=1
              exit
           end if
        end do
     end do
  end if

  check_NaN=nanst

  return

end function check_NaN

subroutine conv_unit_2d( ivlon, ivlat, ovlon, ovlat, coe_uni, undef )
  implicit none
  double precision, intent(in) :: ivlon(:,:)
  double precision, intent(in) :: ivlat(size(ivlon,1),size(ivlon,2))
  double precision, intent(inout) :: ovlon(size(ivlon,1),size(ivlon,2))
  double precision, intent(inout) :: ovlat(size(ivlon,1),size(ivlon,2))
  double precision :: coe_uni
  double precision :: undef
  integer :: ii, jj, ix, jy

  ix=size(ivlon,1)
  jy=size(ivlon,2)

  ovlon=undef
  ovlat=undef

  do jj=1,jy
     do ii=1,ix
        if(ivlon(ii,jj)/=undef.and.ivlat(ii,jj)/=undef)then
           ovlon(ii,jj)=ivlon(ii,jj)*coe_uni
           ovlat(ii,jj)=ivlat(ii,jj)*coe_uni
        end if
     end do
  end do

end subroutine conv_unit_2d

subroutine conv_r4_r8( rval, dval, rundef, dundef )
  implicit none
  real, intent(in) :: rval(:,:)
  double precision, intent(inout) :: dval(size(rval,1),size(rval,2))
  real, intent(in) :: rundef
  double precision, intent(in) :: dundef
  integer :: ii, jj, ix, jy

  ix=size(rval,1)
  jy=size(rval,2)

  dval=dundef

  do jj=1,jy
     do ii=1,ix
        if(rval(ii,jj)/=rundef)then
           dval(ii,jj)=dble(rval(ii,jj))
        end if
     end do
  end do

end subroutine conv_r4_r8

subroutine conv_r8_r4( dval, rval, rundef, dundef )
  implicit none
  double precision, intent(in) :: dval(:,:)
  real, intent(inout) :: rval(size(dval,1),size(dval,2))
  real, intent(in) :: rundef
  double precision, intent(in) :: dundef
  integer :: ii, jj, ix, jy

  ix=size(dval,1)
  jy=size(dval,2)

  rval=rundef

  do jj=1,jy
     do ii=1,ix
        if(dval(ii,jj)/=dundef)then
           rval(ii,jj)=real(dval(ii,jj))
        end if
     end do
  end do

end subroutine conv_r8_r4

subroutine convert_Tbb2Zph( tval, zval, t1d, z1d, undef )
! 与えられた t1d, z1d という鉛直 1 次元の気温, 高度プロファイルから
! 各 2 次元水平面内における tval に対応する高度 zval を内挿処理して返す.
  implicit none
  double precision, intent(in) :: tval(:,:)
  double precision, intent(inout) :: zval(size(tval,1),size(tval,2))
  double precision, intent(in) :: t1d(:)
  double precision, intent(in) :: z1d(size(t1d))
  double precision, intent(in) :: undef
  integer :: ii, jj, kk, ix, jy, kz, itmin
  double precision :: tmin, dzdt

  ix=size(tval,1)
  jy=size(tval,2)
  kz=size(t1d)

  zval=undef

  ! 気温の最低値をとり, これより下の格子でのみ内挿処理を行う.
  call min_val_1d( t1d, itmin, tmin, undef=undef )

  do jj=1,jy
     do ii=1,ix
        if(tval(ii,jj)/=undef)then
           if(tval(ii,jj)>=tmin)then
              do kk=1,itmin-1
                 if((t1d(kk)-tval(ii,jj))*(t1d(kk+1)-tval(ii,jj))<0.0d0)then
                    dzdt=(z1d(kk+1)-z1d(kk))/(t1d(kk+1)-t1d(kk))
                    zval(ii,jj)=z1d(kk)+dzdt*(tval(ii,jj)-t1d(kk))
                    exit
                 end if
              end do
           else  ! replacing tval to tmin
              zval(ii,jj)=z1d(itmin)
           end if
        end if
     end do
  end do

end subroutine convert_Tbb2Zph

subroutine nearest_neighbor_2d( ilon, ilat, iv, olon, olat, ov, undef,  &
  &                             dl_thres )
  ! lon, lat は全て単位 [rad], dl_thres も [rad]. 
  use Math_Const

  implicit none

  double precision, dimension(:,:), intent(in) :: ilon
  double precision, dimension(size(ilon,1),size(ilon,2)), intent(in) :: ilat
  double precision, dimension(size(ilon,1),size(ilon,2)), intent(in) :: iv
  double precision, dimension(:), intent(in) :: olon
  double precision, dimension(:), intent(in) :: olat
  double precision, dimension(size(olon),size(olat)), intent(inout) :: ov
  double precision, intent(in) :: undef
  double precision, intent(in) :: dl_thres

  integer :: k, l, ix, iy, icounter
  integer :: nsi, nti, nxo, nyo
  integer :: ileft, iright, ibot, itop
  integer, dimension(size(olat)) :: boundxs, boundxe
  double precision :: tmprad, fleft, fright, fbot, ftop
  double precision :: drad(size(olon),size(olat))

  nsi=size(ilon,1)
  nti=size(ilon,2)
  nxo=size(olon)
  nyo=size(olat)

  ov=undef
  drad=undef

  do k=1,nti
     do l=1,nsi
        call nearest_neighbor_search_1d( olon, ilon(l,k), ix )
        call nearest_neighbor_search_1d( olat, ilat(l,k), iy )
        if(abs(ilon(l,k)-olon(ix))<=dl_thres.and.  &
  &        abs(ilat(l,k)-olat(iy))<=dl_thres)then
           tmprad=ll2radi( dble(ilon(l,k)), dble(ilat(l,k)),  &
  &                        dble(olon(ix)), dble(olat(iy)), forcef=.true. )
           if(drad(ix,iy)==dble(undef))then
              drad(ix,iy)=tmprad
              ov(ix,iy)=iv(l,k)
           else
              if(tmprad<drad(ix,iy))then
                 drad(ix,iy)=tmprad
                 ov(ix,iy)=iv(l,k)
              else if((tmprad==drad(ix,iy)).and.(ov(ix,iy)==undef).and.  &
  &                   (iv(l,k)/=undef))then
                 ov(ix,iy)=iv(l,k)
              end if
           end if
        end if
     end do
  end do

!-- データ領域のみ探索するフラグ
  boundxs=0
  boundxe=0
  do k=4,nyo-3
     do l=4,nxo-3
        if(ov(l,k)/=undef)then  ! とある k で左端から最初に undef でなくなる点
           boundxs(k)=l
           exit
        end if
     end do
     do l=nxo-3,4,-1
        if(ov(l,k)/=undef)then  ! とある k で右端から最初に undef でなくなる点
           boundxe(k)=l
           exit
        end if
     end do
  end do

  icounter=1

  do while (icounter>0)
     icounter=0
     do k=4,nyo-3
        if(boundxs(k)/=0)then
           do l=boundxs(k)+3,boundxe(k)-3
              if(ov(l,k)==undef)then
                 ileft=0
                 iright=0
                 ibot=0
                 itop=0

                 if(ov(l-1,k)/=undef)then
                    ileft=1
                 else
                    if(ov(l-2,k)/=undef)then
                       ileft=2
                    else
                       if(ov(l-3,k)/=undef)then
                          ileft=3
                       end if
                    end if
                 end if

                 if(ov(l+1,k)/=undef)then
                    iright=1
                 else
                    if(ov(l+2,k)/=undef)then
                       iright=2
                    else
                       if(ov(l+3,k)/=undef)then
                          iright=3
                       end if
                    end if
                 end if

                 if(ov(l,k-1)/=undef)then
                    ibot=1
                 else
                    if(ov(l,k-2)/=undef)then
                       ibot=2
                    else
                       if(ov(l,k-3)/=undef)then
                          ibot=3
                       end if
                    end if
                 end if

                 if(ov(l,k+1)/=undef)then
                    itop=1
                 else
                    if(ov(l,k+2)/=undef)then
                       itop=2
                    else
                       if(ov(l,k+3)/=undef)then
                          itop=3
                       end if
                    end if
                 end if

                 if(ileft/=0.and.iright/=0.and.ibot/=0.and.itop/=0)then
                    fleft=dble(iright)/dble(ileft+iright)
                    fright=dble(ileft)/dble(ileft+iright)
                    fbot=dble(itop)/dble(ibot+itop)
                    ftop=dble(ibot)/dble(ibot+itop)
                    ov(l,k)=0.5d0*(fleft*ov(l-ileft,k)+fright*ov(l+iright,k)  &
  &                               +fbot*ov(l,k-ibot)+ftop*ov(l,k+itop))
                    icounter=icounter+1
                 else
                    if(ileft/=0.and.iright/=0)then
                       fleft=dble(iright)/dble(ileft+iright)
                       fright=dble(ileft)/dble(ileft+iright)
                       ov(l,k)=fleft*ov(l-ileft,k)+fright*ov(l+iright,k)
                       icounter=icounter+1
                    else if(ibot/=0.and.itop/=0)then
                       fbot=dble(itop)/dble(ibot+itop)
                       ftop=dble(ibot)/dble(ibot+itop)
                       ov(l,k)=fbot*ov(l,k-ibot)+ftop*ov(l,k+itop)
                       icounter=icounter+1
                    end if
                 end if

              end if
           end do
        end if
     end do
  end do

end subroutine nearest_neighbor_2d

!-------------------------------
!-------------------------------

end program
