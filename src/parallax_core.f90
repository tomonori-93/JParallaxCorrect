!-----------------------------------------------------------------------
!     Copyright (C) 2020-2026 Satoki Tsujino. All rights reserved.
!-----------------------------------------------------------------------

module parallax_core
!! Fortran library for parallax correction of geostationary meteorological satellites

  implicit none

  double precision, parameter :: pi_dp=3.14159265358979d0

contains

subroutine Parallax_Correct( lon_cld, lat_cld, h_cld, lon_cor, lat_cor,  &
  &                          re, rp, hsat, psat, lsat, missing_value )
!! Core library for parallax correction
  double precision, intent(in) :: lon_cld(:,:)  !! Longitude for each pixel [rad]
  double precision, intent(in) :: lat_cld(size(lon_cld,1),size(lon_cld,2))
                                                   !! Latitude for each pixel [rad]
  double precision, intent(in) :: h_cld(size(lon_cld,1),size(lon_cld,2))
                                                   !! Height for each pixel [m]
  double precision, intent(out) :: lon_cor(size(lon_cld,1),size(lon_cld,2))
                                                 !! Parallax-corrected longitude [rad]
  double precision, intent(out) :: lat_cor(size(lon_cld,1),size(lon_cld,2))
                                                 !! Parallax-corrected latitude [rad]
  !-- Intrinsic parameters for each satellite and constants
  double precision, intent(in) :: re  !! equator_radius =6378.1370d3
  double precision, intent(in) :: rp  !! polar_radius =6356.7523d3
  double precision, intent(in) :: hsat  !! satellite_earth_center_distance [m]
  double precision, intent(in) :: psat  !! satellite_latitude [rad]
  double precision, intent(in) :: lsat  !! satellite_longitude [rad]
  double precision, intent(in) :: missing_value  !! missing value

  double precision :: x0, y0, z0, xs, ys, zs
  double precision :: rrat  !! re/rp
  double precision :: h, lcld, pcld, lcldr, pcldr, lcldc, pcldc
  double precision :: xx0, yy0, zz0, xxs, yys, zzs, xys0, x20, y20, z20, onemr
  double precision :: xyz0, gamcoe, gamcoe2, tparm
  double precision :: ztmp, zeps, parfunc, pardfunc, ztmpa, ztmpb, ztmpc
  double precision :: xa, ya, za, lata, lona
  integer :: i, j, nx, ny

  rrat=re/rp
  onemr=1.0d0-rrat**2

  nx=size(lon_cld,1)
  ny=size(lon_cld,2)

!  if(present(missing_value))then
!     rmissing_value=missing_value
!  else
!     rmissing_value=-100.0d0
!  end if

  do j=1,ny
     do i=1,nx

        pcldr=datan2(dtan(lat_cld(i,j)),rrat**2)  ! geogra -> geocen
        lcldr=lon_cld(i,j)
        h=h_cld(i,j)

        x0=rrat*dsin(0.5d0*pi_dp-pcldr)*dcos(lcldr)
        y0=rrat*dsin(0.5d0*pi_dp-pcldr)*dsin(lcldr)
        z0=dcos(0.5d0*pi_dp-pcldr)
        xs=(hsat/rp)*dcos(lsat)*dcos(psat)
        ys=(hsat/rp)*dsin(lsat)*dcos(psat)
        zs=(hsat/rp)*dsin(psat)
        h=h/rp

        xx0=xs-x0
        yy0=ys-y0
        zz0=zs-z0
        xxs=x0*zs-z0*xs
        yys=y0*zs-z0*ys
        x20=xx0**2
        y20=yy0**2
        z20=zz0**2
        xys0=xxs*xx0+yys*yy0
        xyz0=x20+y20+z20

        if(h/=missing_value.and.z0/=0.0d0)then

           ztmp=z0
           ztmpa=ztmp
           ztmpb=0.8d0*ztmp
           ztmpc=0.5d0*(ztmpa+ztmpb)
           parfunc=(((xxs+ztmpa*xx0)**2+(yys+ztmpa*yy0)**2  &
  &                 +(zz0*rrat*ztmpa)**2-(zz0*rrat)**2)  &
  &                *(1.0d0-onemr*(ztmpa**2))  &
  &                +((ztmpa*h)**2)*((rrat**2)*(x20+y20)+z20)-(h*zz0)**2)**2  &
  &               -4.0d0*((xys0*ztmpa+xyz0*(ztmpa**2)-z20)**2)  &
  &                *(1.0d0-onemr*(ztmpa**2))*((h*rrat)**2)
           pardfunc=(((xxs+ztmpc*xx0)**2+(yys+ztmpc*yy0)**2  &
  &                  +(zz0*rrat*ztmpc)**2-(zz0*rrat)**2)  &
  &                 *(1.0d0-onemr*(ztmpc**2))  &
  &                 +((ztmpc*h)**2)*((rrat**2)*(x20+y20)+z20)-(h*zz0)**2)**2  &
  &                -4.0d0*((xys0*ztmpc+xyz0*(ztmpc**2)-z20)**2)  &
  &                 *(1.0d0-onemr*(ztmpc**2))*((h*rrat)**2)
           do while (abs(ztmpa-ztmpb)>1.0d-12)
              if(parfunc*pardfunc<0.0d0)then
                 ztmpb=ztmpc
                 ztmpc=0.5d0*(ztmpa+ztmpb)
                 pardfunc=(((xxs+ztmpc*xx0)**2+(yys+ztmpc*yy0)**2  &
  &                        +(zz0*rrat*ztmpc)**2-(zz0*rrat)**2)  &
  &                       *(1.0d0-onemr*(ztmpc**2))  &
  &                       +((ztmpc*h)**2)*((rrat**2)*(x20+y20)+z20)-(h*zz0)**2)**2  &
  &                      -4.0d0*((xys0*ztmpc+xyz0*(ztmpc**2)-z20)**2)  &
  &                       *(1.0d0-onemr*(ztmpc**2))*((h*rrat)**2)
              else
                 ztmpa=ztmpc
                 ztmpc=0.5d0*(ztmpa+ztmpb)
                 parfunc=pardfunc
                 pardfunc=(((xxs+ztmpc*xx0)**2+(yys+ztmpc*yy0)**2  &
  &                        +(zz0*rrat*ztmpc)**2-(zz0*rrat)**2)  &
  &                       *(1.0d0-onemr*(ztmpc**2))  &
  &                       +((ztmpc*h)**2)*((rrat**2)*(x20+y20)+z20)-(h*zz0)**2)**2  &
  &                      -4.0d0*((xys0*ztmpc+xyz0*(ztmpc**2)-z20)**2)  &
  &                       *(1.0d0-onemr*(ztmpc**2))*((h*rrat)**2)
              end if
           end do
           ztmp=ztmpc

           za=ztmp
           gamcoe=1.0d0/dsqrt(1.0d0-onemr*(za**2))
           gamcoe2=((rrat+h*gamcoe)/(1.0d0+h*rrat*gamcoe))**2
           tparm=((x0*xx0+y0*yy0+z0*zz0*gamcoe2)/(x20+y20+z20*gamcoe2))  &
  &             *(-1.0d0+dsqrt(1.0d0+((x20+y20+z20*gamcoe2)  &
  &                                  /((x0*xx0+y0*yy0+z0*zz0*gamcoe2)**2))  &
  &                                  *((2.0d0+gamcoe*h*rrat)*gamcoe*h*rrat  &
  &                                   -z0*z0*(gamcoe2-rrat**2))))
           if(tparm<0.0d0.or.tparm>1.0d0)then
              tparm=((x0*xx0+y0*yy0+z0*zz0*gamcoe2)/(x20+y20+z20*gamcoe2))  &
  &                *(-1.0d0-dsqrt(1.0d0+((x20+y20+z20*gamcoe2)  &
  &                                     /((x0*xx0+y0*yy0+z0*zz0*gamcoe2)**2))  &
  &                                     *((2.0d0+gamcoe*h*rrat)*gamcoe*h*rrat  &
  &                                      -z0*z0*(gamcoe2-rrat**2))))
              if(tparm<0.0d0.or.tparm>1.0d0)then
                 tparm=missing_value
                 lon_cor(i,j)=missing_value
                 lat_cor(i,j)=missing_value
                 cycle
!                 write(*,*) "*** ERROR (Parallax_Himawari) ***: tparm is invalid."
              end if
           end if

!        if(zs-z0/=0.0d0)then
!           write(*,*) "tparm check", tparm, (-z0+za+(h*rrat*za)*gamcoe)/(zs-z0)
!           tparm=(-z0+za+(h*rrat*za)*gamcoe)/(zs-z0)
!        else
!           tparm=0.0d0
!        end if
           xa=((1.0d0-tparm)*x0+tparm*xs)/(1.0d0+(h/rrat)*gamcoe)
           ya=((1.0d0-tparm)*y0+tparm*ys)/(1.0d0+(h/rrat)*gamcoe)

           lata=datan2(dtan(0.5d0*pi_dp-dacos(za)),1.0d0/(rrat**2))  ! geocen -> geogra
           lona=datan2(ya,xa)

           lon_cor(i,j)=lona
           lat_cor(i,j)=lata

        else if(h/=missing_value.and.z0==0.0d0)then

           tparm=((xys0-x20-y20)/((xs-x0)**2+(ys-y0)**2))  &
  &             *(-1.0d0+dsqrt(1.0d0+((xs-x0)**2+(ys-y0)**2)  &
  &                                  *(2.0*rrat*h+h**2)/((xys0-x20-y20)**2)))
           if(tparm<0.0d0.or.tparm>1.0d0)then
              tparm=((xys0-x20-y20)/((xs-x0)**2+(ys-y0)**2))  &
  &                *(-1.0d0-dsqrt(1.0d0+((xs-x0)**2+(ys-y0)**2)  &
  &                                     *(2.0*rrat*h+h**2)/((xys0-x20-y20)**2)))
              if(tparm<0.0d0.or.tparm>1.0d0)then
                 tparm=missing_value
                 lon_cor(i,j)=missing_value
                 lat_cor(i,j)=missing_value
                 cycle
!                 write(*,*) "*** ERROR (Parallax_Himawari) ***: tparm is invalid."
              end if
           end if

           xa=((1.0d0-tparm)*x0+tparm*xs)/(1.0d0+(h/rrat))
           ya=((1.0d0-tparm)*y0+tparm*ys)/(1.0d0+(h/rrat))

           lona=datan2(ya,xa)

           lon_cor(i,j)=lona
           lat_cor(i,j)=pcldr

        else

           tparm=missing_value
           lon_cor(i,j)=missing_value
           lat_cor(i,j)=missing_value

        end if

     end do
  end do

end subroutine Parallax_Correct


subroutine tri_interpolation_2d( x_in, y_in, iv, ivad, x_out, y_out, ov,  &
  &                              ovad, missing_value, jflag )
!-- 2 次元不等間隔格子点 (x_in,y_in) 上で定義された物理量 iv を
!-- 同じ次元の等間隔格子 (x_out,y_out) に
!-- 三角点補間 (tri_interpolation) するルーチン (補間値は ov).
!-- [注意]: x_out, y_out は等間隔でなければならない.
  implicit none
  double precision, dimension(:,:), intent(in) :: x_in  ! 補間前のオリジナル座標 x 成分
  double precision, dimension(size(x_in,1),size(x_in,2)), intent(in) :: y_in  ! 補間前のオリジナル座標 y 成分
  double precision, dimension(size(x_in,1),size(x_in,2)), intent(in) :: iv  ! x_in, y_in で定義された物理量
  double precision, dimension(size(x_in,1),size(x_in,2)), intent(in) :: ivad  ! iv 以外に追加物理量
  double precision, dimension(:), intent(in) :: x_out  ! 補間する座標 x 成分 (x_in と同じ単位, 等間隔)
  double precision, dimension(:), intent(in) :: y_out  ! 補間する座標 y 成分 (y_in と同じ単位, 等間隔)
  double precision, dimension(size(x_out),size(y_out)), intent(out) :: ov  ! 補間された値
  double precision, dimension(size(x_out),size(y_out)), intent(out) :: ovad  ! ivad の補間された値
  double precision, intent(in) :: missing_value
  character(1), intent(in) :: jflag  ! 内挿点の値が更新される基準.
                                     ! 'u' : 値が大きいと更新, 'l' : 値が小さいと更新.

  integer :: k, l, m, ix, jy, icounter
  integer :: nsi, nti, nxo, nyo, ixmin, ixmax, jymin, jymax, itmp
  integer, dimension(2) :: isqr
  double precision :: intx_out(size(x_out)), inty_out(size(y_out))
  double precision, dimension(size(x_in,1),size(y_in,2)) :: intx_in, inty_in
  double precision :: dlon, dlat, x_outmin, y_outmin, ov_tmp, ovad_tmp
  double precision, dimension(4) :: sqrlon, sqrlat, sqrval, sqrvalad
  double precision, dimension(3) :: interp_lon, interp_lat, interp_val, interp_valad
  logical :: calc_flag

  nsi=size(x_in,1)
  nti=size(x_in,2)
  nxo=size(x_out)
  nyo=size(y_out)

  dlon=x_out(2)-x_out(1)
  dlat=y_out(2)-y_out(1)
  x_outmin=x_out(1)
  y_outmin=y_out(1)

  ov=missing_value
  ovad=missing_value

!-- x_out, y_out を x_out(1) = 1, y_out(1) = 1 として格子点番号に変換

  intx_out=(/((x_out(ix)-x_outmin)/dlon+1.0d0,ix=1,nxo)/)
  inty_out=(/((y_out(jy)-y_outmin)/dlat+1.0d0,jy=1,nyo)/)

!-- x_in, y_in を x_out, y_out 系での格子点番号 (実数) に変換

  intx_in=missing_value  ! x_in, y_in は missing_value が入っている (特に領域外側).
  inty_in=missing_value
  do l=1,nti
     do k=1,nsi
        if(x_in(k,l)/=missing_value.and.y_in(k,l)/=missing_value)then
           intx_in(k,l)=(x_in(k,l)-x_outmin)/dlon+1.0d0
           inty_in(k,l)=(y_in(k,l)-y_outmin)/dlat+1.0d0
        end if
     end do
  end do

!-- x_in, y_in を左下から順に三角形分割し, その三角形内に
!-- 1. 標的格子 (x_out, y_out) が含まれているかチェック,
!-- 2. 含まれていれば三角形で線形内挿.

  do l=1,nti-1
     do k=1,nsi-1
        sqrlon=(/intx_in(k,l), intx_in(k+1,l), intx_in(k,l+1), intx_in(k+1,l+1)/)
        sqrlat=(/inty_in(k,l), inty_in(k+1,l), inty_in(k,l+1), inty_in(k+1,l+1)/)
        sqrval=(/iv(k,l), iv(k+1,l), iv(k,l+1), iv(k+1,l+1)/)
        sqrvalad=(/ivad(k,l), ivad(k+1,l), ivad(k,l+1), ivad(k+1,l+1)/)

        !-- 対象とする隣接 4 点が全て未定義でないことの確認
        calc_flag=.true.
        do m=1,4
           if(sqrlon(m)==missing_value)then
              calc_flag=.false.
              exit
           end if
           if(sqrlat(m)==missing_value)then
              calc_flag=.false.
              exit
           end if
           if(sqrval(m)==missing_value)then
              calc_flag=.false.
              exit
           end if
        end do
        if(calc_flag.eqv..false.)then  ! 計算しないなら, cycle で次の k,l へ
           cycle
        end if

        !-- 隣接 4 点から三角形に分けるため, 対角線を特定
        !-- (線分の交点の有無で判断)
        !-- selopt で長い方の対角点を取得している. 
        !-- 対角線の短い三角形を使うので, 後の処理のため長い方の対角点を取得.
        call check_square_intersect( sqrlon, sqrlat, isqr, selopt='l' )

        !-- 分けた三角形で三角形内に内挿点候補があるかのチェック
        do m=1,2  ! 三角形は 2 つある.

           icounter=1
           do itmp=1,4
              if(isqr(m)/=itmp)then  ! 対角点以外で m 番目の三角形を構築
                 interp_lon(icounter)=sqrlon(itmp)
                 interp_lat(icounter)=sqrlat(itmp)
                 interp_val(icounter)=sqrval(itmp)
                 interp_valad(icounter)=sqrvalad(itmp)
                 icounter=icounter+1
              end if
           end do

           ixmin=idint(dmin1(interp_lon(1),interp_lon(2),interp_lon(3)))-1
           ixmax=idint(dmax1(interp_lon(1),interp_lon(2),interp_lon(3)))+2
           jymin=idint(dmin1(interp_lat(1),interp_lat(2),interp_lat(3)))-1
           jymax=idint(dmax1(interp_lat(1),interp_lat(2),interp_lat(3)))+2

           if(ixmin<1.or.ixmax>nxo.or.jymin<1.or.jymax>nyo)then
              ! ここの範囲は以下の ix, jy ループの端と連動する.
              !-- 内挿点候補がそもそもない. -> 次のループへ (reported by Tsukada).
              exit
           end if

           if(ixmin<=ixmax.and.jymin<=jymax)then  ! 内挿点候補がある場合.
              do jy=jymin,jymax  ! この jymin, jymax, ixmin, ixmax は ov の点
                 do ix=ixmin,ixmax
                    if(check_triclose( interp_lon, interp_lat,  &
  &                              (/intx_out(ix),inty_out(jy)/) ).eqv..true.)then
                    ! 実際に内挿点がある場合.
                       call tri_interpolation( interp_lon, interp_lat, interp_val,  &
  &                                            (/intx_out(ix), inty_out(jy)/), ov_tmp )
                       call tri_interpolation( interp_lon, interp_lat, interp_valad,  &
  &                                            (/intx_out(ix), inty_out(jy)/), ovad_tmp )
                       if(ov(ix,jy)/=missing_value)then
                          if(jflag(1:1)=='l'.and.(ov(ix,jy)>ov_tmp))then  ! ov_tmp の方が小さいとき更新.
                             ov(ix,jy)=ov_tmp
                             ovad(ix,jy)=ovad_tmp
                          else if(jflag(1:1)=='u'.and.(ov(ix,jy)<ov_tmp))then  ! ov_tmp の方が大きいとき更新.
                             ov(ix,jy)=ov_tmp
                             ovad(ix,jy)=ovad_tmp
                          end if
                       else
                          ov(ix,jy)=ov_tmp
                          ovad(ix,jy)=ovad_tmp
                       end if
                    end if
                 end do
              end do
           end if

        end do

     end do
  end do

end subroutine tri_interpolation_2d


subroutine tri_interpolation( x, y, val, point, oval )
  ! 三角形の要素内の線形内挿ルーチン
  ! 要素内における内挿点 point(1:2) によって分割される領域面積に対して,
  ! oval = \sum^{3}_{i}(S_(i)*val(i)) / S, 
  ! S = \sum^{3}_{i}(S_(i))
  ! 各面積はベクトルの外積が平行四辺形の面積になるという性質より, 
  ! 三角形の各点と要素点の座標から求める.
  implicit none
  double precision, intent(in) :: x(3)    ! 三角形の各頂点 x 座標
  double precision, intent(in) :: y(3)    ! 三角形の各頂点 y 座標
  double precision, intent(in) :: val(3)  ! 三角形の各頂点での値
  double precision, intent(in) :: point(2)! 内挿点 point(1)<=x 座標, point(2)<=y 座標
  double precision, intent(out) :: oval  ! 内挿点での値
  double precision :: Stot, S(3), xp, yp

  xp=point(1)
  yp=point(2)
  Stot=(x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1))
  S(1)=(x(3)-x(2))*(yp-y(3))-(xp-x(3))*(y(3)-y(2))
  S(2)=(x(1)-x(3))*(yp-y(1))-(xp-x(1))*(y(1)-y(3))
  S(3)=(x(2)-x(1))*(yp-y(2))-(xp-x(2))*(y(2)-y(1))

  oval=(val(1)*S(1)+val(2)*S(2)+val(3)*S(3))/Stot

end subroutine tri_interpolation


subroutine check_square_intersect( x, y, inum, selopt )
!-- 四角形の 4 点の座標から短い方の対角点を求める.
!-- selopt = 'l' にすれば長い方が返される.
  implicit none
  double precision, intent(in) :: x(4)  ! 四角形の 4 点 x 座標
  double precision, intent(in) :: y(4)  ! 四角形の 4 点 y 座標
  integer, intent(out) :: inum(2)  ! 最短対角点の座標
  character(1), intent(in) :: selopt  ! 'l' = 対角線の長い方を選ぶ.
                                      ! 's' = 対角線の短い方を選ぶ.
                                      ! デフォルト: 's'
  logical :: flag_long

  flag_long=.false.

  !if(present(selopt))then
     if(selopt(1:1)=='l')then
        flag_long=.true.
     end if
  !end if

  !-- 四角形 (ABCD) の 2 点を結ぶ線のとり方は 3 種類のみ.
  !-- AB-CD, AC-BD, AD-BC

  if(check_intersect( (/x(1),x(2)/), (/y(1),y(2)/),  &
  &                   (/x(3),x(4)/), (/y(3),y(4)/) ).eqv..true.)then  ! AB-CD
     if((x(1)-x(2))**2+(y(1)-y(2))**2<(x(3)-x(4))**2+(y(3)-y(4))**2)then ! AB
        if(flag_long.eqv..true.)then
           inum(1:2)=(/3,4/)
        else
           inum(1:2)=(/1,2/)
        end if
     else  ! CD
        if(flag_long.eqv..true.)then
           inum(1:2)=(/1,2/)
        else
           inum(1:2)=(/3,4/)
        end if
     end if
  else
     if(check_intersect( (/x(1),x(3)/), (/y(1),y(3)/),  &
  &                      (/x(2),x(4)/), (/y(2),y(4)/) ).eqv..true.)then  ! AC-BD
        if((x(1)-x(3))**2+(y(1)-y(3))**2<(x(2)-x(4))**2+(y(2)-y(4))**2)then ! AC
           if(flag_long.eqv..true.)then
              inum(1:2)=(/2,4/)
           else
              inum(1:2)=(/1,3/)
           end if
        else  ! BD
           if(flag_long.eqv..true.)then
              inum(1:2)=(/1,3/)
           else
              inum(1:2)=(/2,4/)
           end if
        end if
     else  ! 上記以外でなければ AD-BC 確定
        if((x(1)-x(4))**2+(y(1)-y(4))**2<(x(2)-x(3))**2+(y(2)-y(3))**2)then ! AC
           if(flag_long.eqv..true.)then
              inum(1:2)=(/2,3/)
           else
              inum(1:2)=(/1,4/)
           end if
        else  ! BD
           if(flag_long.eqv..true.)then
              inum(1:2)=(/1,4/)
           else
              inum(1:2)=(/2,3/)
           end if
        end if
     end if
  end if

end subroutine check_square_intersect


logical function check_intersect( x1, y1, x2, y2 )
! 2 本の線分の交差を判定する. 
! 求め方: 与えられた {x1,y1} と {x2,y2} で定義される 2 本の線分が
! 交差していることを線分の交差判定から求める. 

  implicit none

  double precision, intent(in) :: x1(2)  ! 1 本目の線分 x 方向格子点番号
  double precision, intent(in) :: y1(2)  ! 1 本目の線分 y 方向格子点番号
  double precision, intent(in) :: x2(2)  ! 2 本目の線分 x 方向格子点番号
  double precision, intent(in) :: y2(2)  ! 2 本目の線分 y 方向格子点番号
  double precision :: xa, ya, xb, yb, xc, yc, xd, yd, t1, t2

  check_intersect=.false.

  xa=x1(1)
  xb=x1(2)
  xc=x2(1)
  xd=x2(2)
  ya=y1(1)
  yb=y1(2)
  yc=y2(1)
  yd=y2(2)

  !-- {x1,y1} を通る直線と, 線分 {x2,y2} の交差判定
  t1=(xa-xb)*(yc-ya)+(ya-yb)*(xa-xc)
  t2=(xa-xb)*(yd-ya)+(ya-yb)*(xa-xd)
  if(t1*t2<0.0d0)then
     !-- {x2,y2} を通る直線と, 線分 {x1,y1} の交差判定
     t1=(xc-xd)*(ya-yc)+(yc-yd)*(xc-xa)
     t2=(xc-xd)*(yb-yc)+(yc-yd)*(xc-xb)
     if(t1*t2<0.0d0)then  ! 両方負でなければ両方の線分が交差していない.
        check_intersect=.true.
        return
     end if
  end if

  return

end function check_intersect


logical function check_triclose( xposi, yposi, ival )
! 三角形で囲まれた閉曲線領域内に ival 点があるかどうかをチェックする.
! この判定は与えられる xposi, yposi が時計回り, 反時計回りどちらで
! 定義されていても, 構成する三角形の重心から判断するので問題ない.
! 求め方: ival と重心のなす線分が, 三角形の各辺のいずれとも交差
!         していないことを判断基準にする. 線分の交差判定は 1 線分と
!         1 直線の交差判定を交互に行う.

  implicit none

  double precision, intent(in) :: xposi(3)  ! 三角形頂点の x 方向格子点番号
  double precision, intent(in) :: yposi(3)  ! 三角形頂点の y 方向格子点番号
  double precision, intent(in) :: ival(2)  ! チェックする点の x, y 格子点番号
  integer :: ii
  double precision, dimension(4) :: xt, yt
  double precision :: xg, yg, xi, yi, xa, xb, ya, yb, t1, t2

  xt(1:3)=xposi(1:3)
  yt(1:3)=yposi(1:3)
  xt(4)=xposi(1)
  yt(4)=yposi(1)
  xi=ival(1)
  yi=ival(2)

  check_triclose=.true.

!-- 三角形重心の計算
  xg=(xt(1)+xt(2)+xt(3))/3.0d0
  yg=(yt(1)+yt(2)+yt(3))/3.0d0

!-- 各辺について, 重心と ival の線分との交差を判定.

  do ii=1,3
     xa=xt(ii)
     xb=xt(ii+1)
     ya=yt(ii)
     yb=yt(ii+1)

     !-- 三角形の一辺を通る直線と, 三角形重心と判定点の線分の交差判定
     t1=(xa-xb)*(yg-ya)+(ya-yb)*(xa-xg)
     t2=(xa-xb)*(yi-ya)+(ya-yb)*(xa-xi)
     if(t1*t2<0.0d0)then
        !-- 三角形重心と判定点を通る直線と, 三角形の一辺の交差判定
        t1=(xg-xi)*(ya-yg)+(yg-yi)*(xg-xa)
        t2=(xg-xi)*(yb-yg)+(yg-yi)*(xg-xb)
        if(t1*t2<0.0d0)then  ! 両方負でなければ両方の線分が交差していない.
           check_triclose=.false.
           return
        end if
     end if

  end do

  return

end function check_triclose


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
  tmin=t1d(1)
  itmin=1
  do kk=2,kz
     if(tmin>t1d(kk))then
        tmin=t1d(kk)
        itmin=kk
     end if
  end do

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


end module parallax_core
