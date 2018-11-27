c     mortran 2.0     (version of 7/04/75 mod 7/4/87 (ajc))              
      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)          zipper
      implicit double precision(a-h,o-z)                                    zipper
      data sml0,eps0,big0,mnlam0,rsqmax0,pmin0,exmx0  /1.0d-5,1.0d-6,9.9    zipper
     *d35,5,0.999,1.0d-9,250.0/
      sml=sml0                                                              zipper
      eps=eps0                                                              zipper
      big=big0                                                              zipper
      mnlam=mnlam0                                                          zipper
      rsqmax=rsqmax0                                                        zipper
      pmin=pmin0                                                            zipper
      exmx=exmx0                                                            zipper
      return                                                                zipper
      entry chg_fract_dev(arg)                                              zipper
      sml0=arg                                                              zipper
      return                                                                zipper
      entry chg_dev_max(arg)                                                zipper
      rsqmax0=arg                                                           zipper
      return                                                                zipper
      entry chg_min_flmin(arg)                                              zipper
      eps0=arg                                                              zipper
      return                                                                zipper
      entry chg_big(arg)                                                    zipper
      big0=arg                                                              zipper
      return                                                                zipper
      entry chg_min_lambdas(irg)                                            zipper
      mnlam0=irg                                                            zipper
      return                                                                zipper
      entry chg_min_null_prob(arg)                                          zipper
      pmin0=arg                                                             zipper
      return                                                                zipper
      entry chg_max_exp(arg)                                                zipper
      exmx0=arg                                                             zipper
      return                                                                zipper
      end                                                                   zipper
      subroutine elnet  (ka,parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,u    zipper
     *lam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                    zipper
      double precision x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni)     zipper
      double precision ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)              zipper
      integer jd(*),ia(nx),nin(nlam)                                        zipper
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 10021                                     zipper
      jerr=10000                                                            zipper
      return                                                                zipper
10021 continue                                                              zipper
      allocate(vq(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      vq=max(0d0,vp)                                                        zipper
      vq=vq*ni/sum(vq)                                                      zipper
      if(ka .ne. 1)goto 10041                                               zipper
      call elnetu  (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,    zipper
     *isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 10051                                                            zipper
10041 continue                                                              zipper
      call elnetn (parm,no,ni,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam,thr,i    zipper
     *sd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
10051 continue                                                              zipper
10031 continue                                                              zipper
      deallocate(vq)                                                        zipper
      return                                                                zipper
      end                                                                   zipper
      subroutine elnetu  (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ula    zipper
     *m,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                    zipper
      double precision x(no,ni),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)      zipper
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)             zipper
      integer jd(*),ia(nx),nin(nlam)                                        zipper
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam           
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                           zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(xm(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(xs(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(ju(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(xv(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(vlam(1:nlam),stat=jerr)                                      zipper
      if(jerr.ne.0) return                                                  zipper
      call chkvars(no,ni,x,ju)                                              zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  zipper
      if(maxval(ju) .gt. 0)goto 10071                                       zipper
      jerr=7777                                                             zipper
      return                                                                zipper
10071 continue                                                              zipper
      call standard(no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr)          zipper
      if(jerr.ne.0) return                                                  zipper
      cl=cl/ys                                                              zipper
      if(isd .le. 0)goto 10091                                              zipper
10100 do 10101 j=1,ni                                                       zipper
      cl(:,j)=cl(:,j)*xs(j)                                                 zipper
10101 continue                                                              zipper
10102 continue                                                              zipper
10091 continue                                                              zipper
      if(flmin.ge.1.0) vlam=ulam/ys                                         zipper
      call elnet1(parm,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,vlam,thr,maxi    zipper
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  zipper
10110 do 10111 k=1,lmu                                                      zipper
      alm(k)=ys*alm(k)                                                      zipper
      nk=nin(k)                                                             zipper
10120 do 10121 l=1,nk                                                       zipper
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          zipper
10121 continue                                                              zipper
10122 continue                                                              zipper
      a0(k)=0.0                                                             zipper
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))           zipper
10111 continue                                                              zipper
10112 continue                                                              zipper
      deallocate(xm,xs,g,ju,xv,vlam)                                        zipper
      return                                                                zipper
      end                                                                   zipper
      subroutine standard (no,ni,x,y,w,isd,intr,ju,g,xm,xs,ym,ys,xv,jerr    zipper
     *)
      implicit double precision(a-h,o-z)                                    zipper
      double precision x(no,ni),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)      zipper
      integer ju(ni)                                                        zipper
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                           zipper
      if(jerr.ne.0) return                                                  zipper
      w=w/sum(w)                                                            zipper
      v=sqrt(w)                                                             zipper
      if(intr .ne. 0)goto 10141                                             zipper
      ym=0.0                                                                zipper
      y=v*y                                                                 zipper
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                         zipper
      y=y/ys                                                                zipper
10150 do 10151 j=1,ni                                                       zipper
      if(ju(j).eq.0)goto 10151                                              zipper
      xm(j)=0.0                                                             zipper
      x(:,j)=v*x(:,j)                                                       zipper
      xv(j)=dot_product(x(:,j),x(:,j))                                      zipper
      if(isd .eq. 0)goto 10171                                              zipper
      xbq=dot_product(v,x(:,j))**2                                          zipper
      vc=xv(j)-xbq                                                          zipper
      xs(j)=sqrt(vc)                                                        zipper
      x(:,j)=x(:,j)/xs(j)                                                   zipper
      xv(j)=1.0+xbq/vc                                                      zipper
      goto 10181                                                            zipper
10171 continue                                                              zipper
      xs(j)=1.0                                                             zipper
10181 continue                                                              zipper
10161 continue                                                              zipper
10151 continue                                                              zipper
10152 continue                                                              zipper
      goto 10191                                                            zipper
10141 continue                                                              zipper
10200 do 10201 j=1,ni                                                       zipper
      if(ju(j).eq.0)goto 10201                                              zipper
      xm(j)=dot_product(w,x(:,j))                                           zipper
      x(:,j)=v*(x(:,j)-xm(j))                                               zipper
      xv(j)=dot_product(x(:,j),x(:,j))                                      zipper
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        zipper
10201 continue                                                              zipper
10202 continue                                                              zipper
      if(isd .ne. 0)goto 10221                                              zipper
      xs=1.0                                                                zipper
      goto 10231                                                            zipper
10221 continue                                                              zipper
10240 do 10241 j=1,ni                                                       zipper
      if(ju(j).eq.0)goto 10241                                              zipper
      x(:,j)=x(:,j)/xs(j)                                                   zipper
10241 continue                                                              zipper
10242 continue                                                              zipper
      xv=1.0                                                                zipper
10231 continue                                                              zipper
10211 continue                                                              zipper
      ym=dot_product(w,y)                                                   zipper
      y=v*(y-ym)                                                            zipper
      ys=sqrt(dot_product(y,y))                                             zipper
      y=y/ys                                                                zipper
10191 continue                                                              zipper
10131 continue                                                              zipper
      g=0.0                                                                 zipper
10250 do 10251 j=1,ni                                                       zipper
      if(ju(j).ne.0) g(j)=dot_product(y,x(:,j))                             zipper
10251 continue                                                              zipper
10252 continue                                                              zipper
      deallocate(v)                                                         zipper
      return                                                                zipper
      end                                                                   zipper
      subroutine elnet1 (beta,ni,ju,vp,cl,g,no,ne,nx,x,nlam,flmin,ulam,t    zipper
     *hr,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                    zipper
      double precision vp(ni),g(ni),x(no,ni),ulam(nlam),ao(nx,nlam)         zipper
      double precision rsqo(nlam),almo(nlam),xv(ni)                         zipper
      double precision cl(2,ni)                                             zipper
      integer ju(ni),ia(nx),kin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: a,da                      
      integer, dimension (:), allocatable :: mm                                 
      double precision, dimension (:,:), allocatable :: c                       
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      if(jerr.ne.0) return;                                                     
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)                zipper
      allocate(a(1:ni),stat=jerr)                                           zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(mm(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(da(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      bta=beta                                                              zipper
      omb=1.0-bta                                                           zipper
      if(flmin .ge. 1.0)goto 10271                                          zipper
      eqs=max(eps,flmin)                                                    zipper
      alf=eqs**(1.0/(nlam-1))                                               zipper
10271 continue                                                              zipper
      rsq=0.0                                                               zipper
      a=0.0                                                                 zipper
      mm=0                                                                  zipper
      nlp=0                                                                 zipper
      nin=nlp                                                               zipper
      iz=0                                                                  zipper
      mnl=min(mnlam,nlam)                                                   zipper
      alm=0.0                                                               zipper
10280 do 10281 m=1,nlam                                                     zipper
      if(flmin .lt. 1.0)goto 10301                                          zipper
      alm=ulam(m)                                                           zipper
      goto 10291                                                            zipper
10301 if(m .le. 2)goto 10311                                                zipper
      alm=alm*alf                                                           zipper
      goto 10291                                                            zipper
10311 if(m .ne. 1)goto 10321                                                zipper
      alm=big                                                               zipper
      goto 10331                                                            zipper
10321 continue                                                              zipper
      alm=0.0                                                               zipper
10340 do 10341 j=1,ni                                                       zipper
      if(ju(j).eq.0)goto 10341                                              zipper
      if(vp(j).le.0.0)goto 10341                                            zipper
      alm=max(alm,abs(g(j))/vp(j))                                          zipper
10341 continue                                                              zipper
10342 continue                                                              zipper
      alm=alf*alm/max(bta,1.0d-3)                                           zipper
10331 continue                                                              zipper
10291 continue                                                              zipper
      dem=alm*omb                                                           zipper
      ab=alm*bta                                                            zipper
      rsq0=rsq                                                              zipper
      jz=1                                                                  zipper
10350 continue                                                              zipper
10351 continue                                                              zipper
      if(iz*jz.ne.0) go to 10360                                            zipper
      nlp=nlp+1                                                             zipper
      dlx=0.0                                                               zipper
10370 do 10371 k=1,ni                                                       zipper
      if(ju(k).eq.0)goto 10371                                              zipper
      ak=a(k)                                                               zipper
      u=g(k)+ak*xv(k)                                                       zipper
      v=abs(u)-vp(k)*ab                                                     zipper
      a(k)=0.0                                                              zipper
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    zipper
     *em)))
      if(a(k).eq.ak)goto 10371                                              zipper
      if(mm(k) .ne. 0)goto 10391                                            zipper
      nin=nin+1                                                             zipper
      if(nin.gt.nx)goto 10372                                               zipper
10400 do 10401 j=1,ni                                                       zipper
      if(ju(j).eq.0)goto 10401                                              zipper
      if(mm(j) .eq. 0)goto 10421                                            zipper
      c(j,nin)=c(k,mm(j))                                                   zipper
      goto 10401                                                            zipper
10421 continue                                                              zipper
      if(j .ne. k)goto 10441                                                zipper
      c(j,nin)=xv(j)                                                        zipper
      goto 10401                                                            zipper
10441 continue                                                              zipper
      c(j,nin)=dot_product(x(:,j),x(:,k))                                   zipper
10401 continue                                                              zipper
10402 continue                                                              zipper
      mm(k)=nin                                                             zipper
      ia(nin)=k                                                             zipper
10391 continue                                                              zipper
      del=a(k)-ak                                                           zipper
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      zipper
      dlx=max(xv(k)*del**2,dlx)                                             zipper
10450 do 10451 j=1,ni                                                       zipper
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                               zipper
10451 continue                                                              zipper
10452 continue                                                              zipper
10371 continue                                                              zipper
10372 continue                                                              zipper
      if(dlx.lt.thr)goto 10352                                              zipper
      if(nin.gt.nx)goto 10352                                               zipper
      if(nlp .le. maxit)goto 10471                                          zipper
      jerr=-m                                                               zipper
      return                                                                zipper
10471 continue                                                              zipper
10360 continue                                                              zipper
      iz=1                                                                  zipper
      da(1:nin)=a(ia(1:nin))                                                zipper
10480 continue                                                              zipper
10481 continue                                                              zipper
      nlp=nlp+1                                                             zipper
      dlx=0.0                                                               zipper
10490 do 10491 l=1,nin                                                      zipper
      k=ia(l)                                                               zipper
      ak=a(k)                                                               zipper
      u=g(k)+ak*xv(k)                                                       zipper
      v=abs(u)-vp(k)*ab                                                     zipper
      a(k)=0.0                                                              zipper
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    zipper
     *em)))
      if(a(k).eq.ak)goto 10491                                              zipper
      del=a(k)-ak                                                           zipper
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                      zipper
      dlx=max(xv(k)*del**2,dlx)                                             zipper
10500 do 10501 j=1,nin                                                      zipper
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                  zipper
10501 continue                                                              zipper
10502 continue                                                              zipper
10491 continue                                                              zipper
10492 continue                                                              zipper
      if(dlx.lt.thr)goto 10482                                              zipper
      if(nlp .le. maxit)goto 10521                                          zipper
      jerr=-m                                                               zipper
      return                                                                zipper
10521 continue                                                              zipper
      goto 10481                                                            zipper
10482 continue                                                              zipper
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                      zipper
10530 do 10531 j=1,ni                                                       zipper
      if(mm(j).ne.0)goto 10531                                              zipper
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))            zipper
10531 continue                                                              zipper
10532 continue                                                              zipper
      jz=0                                                                  zipper
      goto 10351                                                            zipper
10352 continue                                                              zipper
      if(nin .le. nx)goto 10551                                             zipper
      jerr=-10000-m                                                         zipper
      goto 10282                                                            zipper
10551 continue                                                              zipper
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 zipper
      kin(m)=nin                                                            zipper
      rsqo(m)=rsq                                                           zipper
      almo(m)=alm                                                           zipper
      lmu=m                                                                 zipper
      if(m.lt.mnl)goto 10281                                                zipper
      if(flmin.ge.1.0)goto 10281                                            zipper
      me=0                                                                  zipper
10560 do 10561 j=1,nin                                                      zipper
      if(ao(j,m).ne.0.0) me=me+1                                            zipper
10561 continue                                                              zipper
10562 continue                                                              zipper
      if(me.gt.ne)goto 10282                                                zipper
      if(rsq-rsq0.lt.sml*rsq)goto 10282                                     zipper
      if(rsq.gt.rsqmax)goto 10282                                           zipper
10281 continue                                                              zipper
10282 continue                                                              zipper
      deallocate(a,mm,c,da)                                                 zipper
      return                                                                zipper
      end                                                                   zipper
      subroutine elnetn (parm,no,ni,x,y,w,jd,vp,cl,ne,nx,nlam,flmin,ulam    zipper
     *,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                    zipper
      double precision vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)      zipper
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)             zipper
      integer jd(*),ia(nx),nin(nlam)                                        zipper
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam             
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(xs(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(ju(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(xv(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                  zipper
      allocate(vlam(1:nlam),stat=jerr)                                      zipper
      if(jerr.ne.0) return                                                  zipper
      call chkvars(no,ni,x,ju)                                              zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  zipper
      if(maxval(ju) .gt. 0)goto 10581                                       zipper
      jerr=7777                                                             zipper
      return                                                                zipper
10581 continue                                                              zipper
      call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)           zipper
      if(jerr.ne.0) return                                                  zipper
      cl=cl/ys                                                              zipper
      if(isd .le. 0)goto 10601                                              zipper
10610 do 10611 j=1,ni                                                       zipper
      cl(:,j)=cl(:,j)*xs(j)                                                 zipper
10611 continue                                                              zipper
10612 continue                                                              zipper
10601 continue                                                             zipper
      if(flmin.ge.1.0) vlam=ulam/ys                                        zipper
      call elnet2(parm,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam,thr,maxi   zipper
     *t,xv,  lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 zipper
10620 do 10621 k=1,lmu                                                     zipper
      alm(k)=ys*alm(k)                                                     zipper
      nk=nin(k)                                                            zipper
10630 do 10631 l=1,nk                                                      zipper
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         zipper
10631 continue                                                             zipper
10632 continue                                                             zipper
      a0(k)=0.0                                                            zipper
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          zipper
10621 continue                                                             zipper
10622 continue                                                             zipper
      deallocate(xm,xs,ju,xv,vlam)                                         zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine standard1 (no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)   zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)           zipper
      integer ju(ni)                                                       zipper
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      w=w/sum(w)                                                           zipper
      v=sqrt(w)                                                            zipper
      if(intr .ne. 0)goto 10651                                            zipper
      ym=0.0                                                               zipper
      y=v*y                                                                zipper
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                        zipper
      y=y/ys                                                               zipper
10660 do 10661 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 10661                                             zipper
      xm(j)=0.0                                                            zipper
      x(:,j)=v*x(:,j)                                                      zipper
      xv(j)=dot_product(x(:,j),x(:,j))                                     zipper
      if(isd .eq. 0)goto 10681                                             zipper
      xbq=dot_product(v,x(:,j))**2                                         zipper
      vc=xv(j)-xbq                                                         zipper
      xs(j)=sqrt(vc)                                                       zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
      xv(j)=1.0+xbq/vc                                                     zipper
      goto 10691                                                           zipper
10681 continue                                                             zipper
      xs(j)=1.0                                                            zipper
10691 continue                                                             zipper
10671 continue                                                             zipper
10661 continue                                                             zipper
10662 continue                                                             zipper
      go to 10700                                                          zipper
10651 continue                                                             zipper
10710 do 10711 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 10711                                             zipper
      xm(j)=dot_product(w,x(:,j))                                          zipper
      x(:,j)=v*(x(:,j)-xm(j))                                              zipper
      xv(j)=dot_product(x(:,j),x(:,j))                                     zipper
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       zipper
10711 continue                                                             zipper
10712 continue                                                             zipper
      if(isd .ne. 0)goto 10731                                             zipper
      xs=1.0                                                               zipper
      goto 10741                                                           zipper
10731 continue                                                             zipper
10750 do 10751 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 10751                                             zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
10751 continue                                                             zipper
10752 continue                                                             zipper
      xv=1.0                                                               zipper
10741 continue                                                             zipper
10721 continue                                                             zipper
      ym=dot_product(w,y)                                                  zipper
      y=v*(y-ym)                                                           zipper
      ys=sqrt(dot_product(y,y))                                            zipper
      y=y/ys                                                               zipper
10700 continue                                                             zipper
      deallocate(v)                                                        zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine elnet2(beta,ni,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,ulam,th   zipper
     *r,maxit,xv,  lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam)        zipper
      double precision rsqo(nlam),almo(nlam),xv(ni)                        zipper
      double precision cl(2,ni)                                            zipper
      integer ju(ni),ia(nx),kin(nlam)                                      zipper
      double precision, dimension (:), allocatable :: a,g                       
      integer, dimension (:), allocatable :: mm,ix                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               zipper
      allocate(a(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(g(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ix(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      bta=beta                                                             zipper
      omb=1.0-bta                                                          zipper
      ix=0                                                                 zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 10771                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
10771 continue                                                             zipper
      rsq=0.0                                                              zipper
      a=0.0                                                                zipper
      mm=0                                                                 zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      iz=0                                                                 zipper
      mnl=min(mnlam,nlam)                                                  zipper
      alm=0.0                                                              zipper
10780 do 10781 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 10781                                             zipper
      g(j)=abs(dot_product(y,x(:,j)))                                      zipper
10781 continue                                                             zipper
10782 continue                                                             zipper
10790 do 10791 m=1,nlam                                                    zipper
      alm0=alm                                                             zipper
      if(flmin .lt. 1.0)goto 10811                                         zipper
      alm=ulam(m)                                                          zipper
      goto 10801                                                           zipper
10811 if(m .le. 2)goto 10821                                               zipper
      alm=alm*alf                                                          zipper
      goto 10801                                                           zipper
10821 if(m .ne. 1)goto 10831                                               zipper
      alm=big                                                              zipper
      goto 10841                                                           zipper
10831 continue                                                             zipper
      alm0=0.0                                                             zipper
10850 do 10851 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 10851                                             zipper
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           zipper
10851 continue                                                             zipper
10852 continue                                                             zipper
      alm0=alm0/max(bta,1.0d-3)                                            zipper
      alm=alf*alm0                                                         zipper
10841 continue                                                             zipper
10801 continue                                                             zipper
      dem=alm*omb                                                          zipper
      ab=alm*bta                                                           zipper
      rsq0=rsq                                                             zipper
      jz=1                                                                 zipper
      tlam=bta*(2.0*alm-alm0)                                              zipper
10860 do 10861 k=1,ni                                                      zipper
      if(ix(k).eq.1)goto 10861                                             zipper
      if(ju(k).eq.0)goto 10861                                             zipper
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       zipper
10861 continue                                                             zipper
10862 continue                                                             zipper
10870 continue                                                             zipper
10871 continue                                                             zipper
      if(iz*jz.ne.0) go to 10360                                           zipper
10880 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
10890 do 10891 k=1,ni                                                      zipper
      if(ix(k).eq.0)goto 10891                                             zipper
      gk=dot_product(y,x(:,k))                                             zipper
      ak=a(k)                                                              zipper
      u=gk+ak*xv(k)                                                        zipper
      v=abs(u)-vp(k)*ab                                                    zipper
      a(k)=0.0                                                             zipper
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   zipper
     *em)))
      if(a(k).eq.ak)goto 10891                                             zipper
      if(mm(k) .ne. 0)goto 10911                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 10892                                              zipper
      mm(k)=nin                                                            zipper
      ia(nin)=k                                                            zipper
10911 continue                                                             zipper
      del=a(k)-ak                                                          zipper
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       zipper
      y=y-del*x(:,k)                                                       zipper
      dlx=max(xv(k)*del**2,dlx)                                            zipper
10891 continue                                                             zipper
10892 continue                                                             zipper
      if(nin.gt.nx)goto 10872                                              zipper
      if(dlx .ge. thr)goto 10931                                           zipper
      ixx=0                                                                zipper
10940 do 10941 k=1,ni                                                      zipper
      if(ix(k).eq.1)goto 10941                                             zipper
      if(ju(k).eq.0)goto 10941                                             zipper
      g(k)=abs(dot_product(y,x(:,k)))                                      zipper
      if(g(k) .le. ab*vp(k))goto 10961                                     zipper
      ix(k)=1                                                              zipper
      ixx=1                                                                zipper
10961 continue                                                             zipper
10941 continue                                                             zipper
10942 continue                                                             zipper
      if(ixx.eq.1) go to 10880                                             zipper
      goto 10872                                                           zipper
10931 continue                                                             zipper
      if(nlp .le. maxit)goto 10981                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
10981 continue                                                             zipper
10360 continue                                                             zipper
      iz=1                                                                 zipper
10990 continue                                                             zipper
10991 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
11000 do 11001 l=1,nin                                                     zipper
      k=ia(l)                                                              zipper
      gk=dot_product(y,x(:,k))                                             zipper
      ak=a(k)                                                              zipper
      u=gk+ak*xv(k)                                                        zipper
      v=abs(u)-vp(k)*ab                                                    zipper
      a(k)=0.0                                                             zipper
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   zipper
     *em)))
      if(a(k).eq.ak)goto 11001                                             zipper
      del=a(k)-ak                                                          zipper
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       zipper
      y=y-del*x(:,k)                                                       zipper
      dlx=max(xv(k)*del**2,dlx)                                            zipper
11001 continue                                                             zipper
11002 continue                                                             zipper
      if(dlx.lt.thr)goto 10992                                             zipper
      if(nlp .le. maxit)goto 11021                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
11021 continue                                                             zipper
      goto 10991                                                           zipper
10992 continue                                                             zipper
      jz=0                                                                 zipper
      goto 10871                                                           zipper
10872 continue                                                             zipper
      if(nin .le. nx)goto 11041                                            zipper
      jerr=-10000-m                                                        zipper
      goto 10792                                                           zipper
11041 continue                                                             zipper
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                zipper
      kin(m)=nin                                                           zipper
      rsqo(m)=rsq                                                          zipper
      almo(m)=alm                                                          zipper
      lmu=m                                                                zipper
      if(m.lt.mnl)goto 10791                                               zipper
      if(flmin.ge.1.0)goto 10791                                           zipper
      me=0                                                                 zipper
11050 do 11051 j=1,nin                                                     zipper
      if(ao(j,m).ne.0.0) me=me+1                                           zipper
11051 continue                                                             zipper
11052 continue                                                             zipper
      if(me.gt.ne)goto 10792                                               zipper
      if(rsq-rsq0.lt.sml*rsq)goto 10792                                    zipper
      if(rsq.gt.rsqmax)goto 10792                                          zipper
10791 continue                                                             zipper
10792 continue                                                             zipper
      deallocate(a,mm,g,ix)                                                zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine chkvars(no,ni,x,ju)                                       zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni)                                            zipper
      integer ju(ni)                                                       zipper
11060 do 11061 j=1,ni                                                      zipper
      ju(j)=0                                                              zipper
      t=x(1,j)                                                             zipper
11070 do 11071 i=2,no                                                      zipper
      if(x(i,j).eq.t)goto 11071                                            zipper
      ju(j)=1                                                              zipper
      goto 11072                                                           zipper
11071 continue                                                             zipper
11072 continue                                                             zipper
11061 continue                                                             zipper
11062 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine uncomp(ni,ca,ia,nin,a)                                    zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision ca(*),a(ni)                                         zipper
      integer ia(*)                                                        zipper
      a=0.0                                                                zipper
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                  zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine modval(a0,ca,ia,nin,n,x,f)                                zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision ca(nin),x(n,*),f(n)                                 zipper
      integer ia(nin)                                                      zipper
      f=a0                                                                 zipper
      if(nin.le.0) return                                                  zipper
11080 do 11081 i=1,n                                                       zipper
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      zipper
11081 continue                                                             zipper
11082 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spelnet  (ka,parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam   zipper
     *,flmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)         zipper
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            zipper
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           zipper
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 11101                                    zipper
      jerr=10000                                                           zipper
      return                                                               zipper
11101 continue                                                             zipper
      allocate(vq(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      vq=max(0d0,vp)                                                       zipper
      vq=vq*ni/sum(vq)                                                     zipper
      if(ka .ne. 1)goto 11121                                              zipper
      call spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,u   zipper
     *lam,thr,isd,  intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      goto 11131                                                           zipper
11121 continue                                                             zipper
      call spelnetn (parm,no,ni,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,flmin,ul   zipper
     *am,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
11131 continue                                                             zipper
11111 continue                                                             zipper
      deallocate(vq)                                                       zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spelnetu  (parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,f   zipper
     *lmin,ulam,thr,isd,intr,  maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),w(no),vp(ni),ulam(nlam),cl(2,ni)         zipper
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            zipper
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           zipper
      double precision, dimension (:), allocatable :: xm,xs,g,xv,vlam           
      integer, dimension (:), allocatable :: ju                                 
      allocate(g(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xv(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(vlam(1:nlam),stat=jerr)                                     zipper
      if(jerr.ne.0) return                                                 zipper
      call spchkvars(no,ni,x,ix,ju)                                        zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 11151                                      zipper
      jerr=7777                                                            zipper
      return                                                               zipper
11151 continue                                                             zipper
      call spstandard(no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys,xv,jer   zipper
     *r)
      if(jerr.ne.0) return                                                 zipper
      cl=cl/ys                                                             zipper
      if(isd .le. 0)goto 11171                                             zipper
11180 do 11181 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
11181 continue                                                             zipper
11182 continue                                                             zipper
11171 continue                                                             zipper
      if(flmin.ge.1.0) vlam=ulam/ys                                        zipper
      call spelnet1(parm,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   zipper
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 zipper
11190 do 11191 k=1,lmu                                                     zipper
      alm(k)=ys*alm(k)                                                     zipper
      nk=nin(k)                                                            zipper
11200 do 11201 l=1,nk                                                      zipper
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         zipper
11201 continue                                                             zipper
11202 continue                                                             zipper
      a0(k)=0.0                                                            zipper
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          zipper
11191 continue                                                             zipper
11192 continue                                                             zipper
      deallocate(xm,xs,g,ju,xv,vlam)                                       zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spstandard (no,ni,x,ix,jx,y,w,ju,isd,intr,g,xm,xs,ym,ys   zipper
     *,xv,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),w(no),g(ni),xm(ni),xs(ni),xv(ni)         zipper
      integer ix(*),jx(*),ju(ni)                                           zipper
      w=w/sum(w)                                                           zipper
      if(intr .ne. 0)goto 11221                                            zipper
      ym=0.0                                                               zipper
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     zipper
      y=y/ys                                                               zipper
11230 do 11231 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11231                                             zipper
      xm(j)=0.0                                                            zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          zipper
      if(isd .eq. 0)goto 11251                                             zipper
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            zipper
      vc=xv(j)-xbq                                                         zipper
      xs(j)=sqrt(vc)                                                       zipper
      xv(j)=1.0+xbq/vc                                                     zipper
      goto 11261                                                           zipper
11251 continue                                                             zipper
      xs(j)=1.0                                                            zipper
11261 continue                                                             zipper
11241 continue                                                             zipper
11231 continue                                                             zipper
11232 continue                                                             zipper
      goto 11271                                                           zipper
11221 continue                                                             zipper
11280 do 11281 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11281                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             zipper
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 zipper
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       zipper
11281 continue                                                             zipper
11282 continue                                                             zipper
      if(isd .ne. 0)goto 11301                                             zipper
      xs=1.0                                                               zipper
      goto 11311                                                           zipper
11301 continue                                                             zipper
      xv=1.0                                                               zipper
11311 continue                                                             zipper
11291 continue                                                             zipper
      ym=dot_product(w,y)                                                  zipper
      y=y-ym                                                               zipper
      ys=sqrt(dot_product(w,y**2))                                         zipper
      y=y/ys                                                               zipper
11271 continue                                                             zipper
11211 continue                                                             zipper
      g=0.0                                                                zipper
11320 do 11321 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11321                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      g(j)=dot_product(w(jx(jb:je))*y(jx(jb:je)),x(jb:je))/xs(j)           zipper
11321 continue                                                             zipper
11322 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spelnet1(beta,ni,g,no,w,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   zipper
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision g(ni),vp(ni),x(*),ulam(nlam),w(no)                  zipper
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam)                   zipper
      double precision xm(ni),xs(ni),xv(ni),cl(2,ni)                       zipper
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          zipper
      double precision, dimension (:), allocatable :: a,da                      
      integer, dimension (:), allocatable :: mm                                 
      double precision, dimension (:,:), allocatable :: c                       
      allocate(c(1:ni,1:nx),stat=jerr)                                          
      if(jerr.ne.0) return;                                                     
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               zipper
      allocate(a(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(da(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      bta=beta                                                             zipper
      omb=1.0-bta                                                          zipper
      alm=0.0                                                              zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 11341                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
11341 continue                                                             zipper
      rsq=0.0                                                              zipper
      a=0.0                                                                zipper
      mm=0                                                                 zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      iz=0                                                                 zipper
      mnl=min(mnlam,nlam)                                                  zipper
11350 do 11351 m=1,nlam                                                    zipper
      if(flmin .lt. 1.0)goto 11371                                         zipper
      alm=ulam(m)                                                          zipper
      goto 11361                                                           zipper
11371 if(m .le. 2)goto 11381                                               zipper
      alm=alm*alf                                                          zipper
      goto 11361                                                           zipper
11381 if(m .ne. 1)goto 11391                                               zipper
      alm=big                                                              zipper
      goto 11401                                                           zipper
11391 continue                                                             zipper
      alm=0.0                                                              zipper
11410 do 11411 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11411                                             zipper
      if(vp(j).le.0.0)goto 11411                                           zipper
      alm=max(alm,abs(g(j))/vp(j))                                         zipper
11411 continue                                                             zipper
11412 continue                                                             zipper
      alm=alf*alm/max(bta,1.0d-3)                                          zipper
11401 continue                                                             zipper
11361 continue                                                             zipper
      dem=alm*omb                                                          zipper
      ab=alm*bta                                                           zipper
      rsq0=rsq                                                             zipper
      jz=1                                                                 zipper
11420 continue                                                             zipper
11421 continue                                                             zipper
      if(iz*jz.ne.0) go to 10360                                           zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
11430 do 11431 k=1,ni                                                      zipper
      if(ju(k).eq.0)goto 11431                                             zipper
      ak=a(k)                                                              zipper
      u=g(k)+ak*xv(k)                                                      zipper
      v=abs(u)-vp(k)*ab                                                    zipper
      a(k)=0.0                                                             zipper
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   zipper
     *em)))
      if(a(k).eq.ak)goto 11431                                             zipper
      if(mm(k) .ne. 0)goto 11451                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 11432                                              zipper
11460 do 11461 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11461                                             zipper
      if(mm(j) .eq. 0)goto 11481                                           zipper
      c(j,nin)=c(k,mm(j))                                                  zipper
      goto 11461                                                           zipper
11481 continue                                                             zipper
      if(j .ne. k)goto 11501                                               zipper
      c(j,nin)=xv(j)                                                       zipper
      goto 11461                                                           zipper
11501 continue                                                             zipper
      c(j,nin)=  (row_prod(j,k,ix,jx,x,w)-xm(j)*xm(k))/(xs(j)*xs(k))       zipper
11461 continue                                                             zipper
11462 continue                                                             zipper
      mm(k)=nin                                                            zipper
      ia(nin)=k                                                            zipper
11451 continue                                                             zipper
      del=a(k)-ak                                                          zipper
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     zipper
      dlx=max(xv(k)*del**2,dlx)                                            zipper
11510 do 11511 j=1,ni                                                      zipper
      if(ju(j).ne.0) g(j)=g(j)-c(j,mm(k))*del                              zipper
11511 continue                                                             zipper
11512 continue                                                             zipper
11431 continue                                                             zipper
11432 continue                                                             zipper
      if(dlx.lt.thr)goto 11422                                             zipper
      if(nin.gt.nx)goto 11422                                              zipper
      if(nlp .le. maxit)goto 11531                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
11531 continue                                                             zipper
10360 continue                                                             zipper
      iz=1                                                                 zipper
      da(1:nin)=a(ia(1:nin))                                               zipper
11540 continue                                                             zipper
11541 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
11550 do 11551 l=1,nin                                                     zipper
      k=ia(l)                                                              zipper
      ak=a(k)                                                              zipper
      u=g(k)+ak*xv(k)                                                      zipper
      v=abs(u)-vp(k)*ab                                                    zipper
      a(k)=0.0                                                             zipper
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   zipper
     *em)))
      if(a(k).eq.ak)goto 11551                                             zipper
      del=a(k)-ak                                                          zipper
      rsq=rsq+del*(2.0*g(k)-del*xv(k))                                     zipper
      dlx=max(xv(k)*del**2,dlx)                                            zipper
11560 do 11561 j=1,nin                                                     zipper
      g(ia(j))=g(ia(j))-c(ia(j),mm(k))*del                                 zipper
11561 continue                                                             zipper
11562 continue                                                             zipper
11551 continue                                                             zipper
11552 continue                                                             zipper
      if(dlx.lt.thr)goto 11542                                             zipper
      if(nlp .le. maxit)goto 11581                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
11581 continue                                                             zipper
      goto 11541                                                           zipper
11542 continue                                                             zipper
      da(1:nin)=a(ia(1:nin))-da(1:nin)                                     zipper
11590 do 11591 j=1,ni                                                      zipper
      if(mm(j).ne.0)goto 11591                                             zipper
      if(ju(j).ne.0) g(j)=g(j)-dot_product(da(1:nin),c(j,1:nin))           zipper
11591 continue                                                             zipper
11592 continue                                                             zipper
      jz=0                                                                 zipper
      goto 11421                                                           zipper
11422 continue                                                             zipper
      if(nin .le. nx)goto 11611                                            zipper
      jerr=-10000-m                                                        zipper
      goto 11352                                                           zipper
11611 continue                                                             zipper
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                zipper
      kin(m)=nin                                                           zipper
      rsqo(m)=rsq                                                          zipper
      almo(m)=alm                                                          zipper
      lmu=m                                                                zipper
      if(m.lt.mnl)goto 11351                                               zipper
      if(flmin.ge.1.0)goto 11351                                           zipper
      me=0                                                                 zipper
11620 do 11621 j=1,nin                                                     zipper
      if(ao(j,m).ne.0.0) me=me+1                                           zipper
11621 continue                                                             zipper
11622 continue                                                             zipper
      if(me.gt.ne)goto 11352                                               zipper
      if(rsq-rsq0.lt.sml*rsq)goto 11352                                    zipper
      if(rsq.gt.rsqmax)goto 11352                                          zipper
11351 continue                                                             zipper
11352 continue                                                             zipper
      deallocate(a,mm,c,da)                                                zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spelnetn(parm,no,ni,x,ix,jx,y,w,jd,vp,cl,ne,nx,nlam,flm   zipper
     *in,ulam,  thr,isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),vp(ni),y(no),w(no),ulam(nlam),cl(2,ni)         zipper
      double precision ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)            zipper
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           zipper
      double precision, dimension (:), allocatable :: xm,xs,xv,vlam             
      integer, dimension (:), allocatable :: ju                                 
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xv(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(vlam(1:nlam),stat=jerr)                                     zipper
      if(jerr.ne.0) return                                                 zipper
      call spchkvars(no,ni,x,ix,ju)                                        zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 11641                                      zipper
      jerr=7777                                                            zipper
      return                                                               zipper
11641 continue                                                             zipper
      call spstandard1(no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,xv,jerr   zipper
     *)
      if(jerr.ne.0) return                                                 zipper
      cl=cl/ys                                                             zipper
      if(isd .le. 0)goto 11661                                             zipper
11670 do 11671 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
11671 continue                                                             zipper
11672 continue                                                             zipper
11661 continue                                                             zipper
      if(flmin.ge.1.0) vlam=ulam/ys                                        zipper
      call spelnet2(parm,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flmin,vla   zipper
     *m,thr,maxit,  xm,xs,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 zipper
11680 do 11681 k=1,lmu                                                     zipper
      alm(k)=ys*alm(k)                                                     zipper
      nk=nin(k)                                                            zipper
11690 do 11691 l=1,nk                                                      zipper
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                         zipper
11691 continue                                                             zipper
11692 continue                                                             zipper
      a0(k)=0.0                                                            zipper
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))          zipper
11681 continue                                                             zipper
11682 continue                                                             zipper
      deallocate(xm,xs,ju,xv,vlam)                                         zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spstandard1 (no,ni,x,ix,jx,y,w,ju,isd,intr,xm,xs,ym,ys,   zipper
     *xv,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),w(no),xm(ni),xs(ni),xv(ni)               zipper
      integer ix(*),jx(*),ju(ni)                                           zipper
      w=w/sum(w)                                                           zipper
      if(intr .ne. 0)goto 11711                                            zipper
      ym=0.0                                                               zipper
      ys=sqrt(dot_product(w,y**2)-dot_product(w,y)**2)                     zipper
      y=y/ys                                                               zipper
11720 do 11721 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11721                                             zipper
      xm(j)=0.0                                                            zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          zipper
      if(isd .eq. 0)goto 11741                                             zipper
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            zipper
      vc=xv(j)-xbq                                                         zipper
      xs(j)=sqrt(vc)                                                       zipper
      xv(j)=1.0+xbq/vc                                                     zipper
      goto 11751                                                           zipper
11741 continue                                                             zipper
      xs(j)=1.0                                                            zipper
11751 continue                                                             zipper
11731 continue                                                             zipper
11721 continue                                                             zipper
11722 continue                                                             zipper
      return                                                               zipper
11711 continue                                                             zipper
11760 do 11761 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11761                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             zipper
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 zipper
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       zipper
11761 continue                                                             zipper
11762 continue                                                             zipper
      if(isd .ne. 0)goto 11781                                             zipper
      xs=1.0                                                               zipper
      goto 11791                                                           zipper
11781 continue                                                             zipper
      xv=1.0                                                               zipper
11791 continue                                                             zipper
11771 continue                                                             zipper
      ym=dot_product(w,y)                                                  zipper
      y=y-ym                                                               zipper
      ys=sqrt(dot_product(w,y**2))                                         zipper
      y=y/ys                                                               zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spelnet2(beta,ni,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,nlam,flm   zipper
     *in,ulam,  thr,maxit,xm,xs,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision y(no),w(no),x(*),vp(ni),ulam(nlam),cl(2,ni)         zipper
      double precision ao(nx,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni),x   zipper
     *v(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          zipper
      double precision, dimension (:), allocatable :: a,g                       
      integer, dimension (:), allocatable :: mm,iy                              
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               zipper
      allocate(a(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(g(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(iy(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      bta=beta                                                             zipper
      omb=1.0-bta                                                          zipper
      alm=0.0                                                              zipper
      iy=0                                                                 zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 11811                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
11811 continue                                                             zipper
      rsq=0.0                                                              zipper
      a=0.0                                                                zipper
      mm=0                                                                 zipper
      o=0.0                                                                zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      iz=0                                                                 zipper
      mnl=min(mnlam,nlam)                                                  zipper
11820 do 11821 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11821                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    zipper
11821 continue                                                             zipper
11822 continue                                                             zipper
11830 do 11831 m=1,nlam                                                    zipper
      alm0=alm                                                             zipper
      if(flmin .lt. 1.0)goto 11851                                         zipper
      alm=ulam(m)                                                          zipper
      goto 11841                                                           zipper
11851 if(m .le. 2)goto 11861                                               zipper
      alm=alm*alf                                                          zipper
      goto 11841                                                           zipper
11861 if(m .ne. 1)goto 11871                                               zipper
      alm=big                                                              zipper
      goto 11881                                                           zipper
11871 continue                                                             zipper
      alm0=0.0                                                             zipper
11890 do 11891 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 11891                                             zipper
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           zipper
11891 continue                                                             zipper
11892 continue                                                             zipper
      alm0=alm0/max(bta,1.0d-3)                                            zipper
      alm=alf*alm0                                                         zipper
11881 continue                                                             zipper
11841 continue                                                             zipper
      dem=alm*omb                                                          zipper
      ab=alm*bta                                                           zipper
      rsq0=rsq                                                             zipper
      jz=1                                                                 zipper
      tlam=bta*(2.0*alm-alm0)                                              zipper
11900 do 11901 k=1,ni                                                      zipper
      if(iy(k).eq.1)goto 11901                                             zipper
      if(ju(k).eq.0)goto 11901                                             zipper
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       zipper
11901 continue                                                             zipper
11902 continue                                                             zipper
11910 continue                                                             zipper
11911 continue                                                             zipper
      if(iz*jz.ne.0) go to 10360                                           zipper
10880 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
11920 do 11921 k=1,ni                                                      zipper
      if(iy(k).eq.0)goto 11921                                             zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           zipper
      ak=a(k)                                                              zipper
      u=gk+ak*xv(k)                                                        zipper
      v=abs(u)-vp(k)*ab                                                    zipper
      a(k)=0.0                                                             zipper
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   zipper
     *em)))
      if(a(k).eq.ak)goto 11921                                             zipper
      if(mm(k) .ne. 0)goto 11941                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 11922                                              zipper
      mm(k)=nin                                                            zipper
      ia(nin)=k                                                            zipper
11941 continue                                                             zipper
      del=a(k)-ak                                                          zipper
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       zipper
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         zipper
      o=o+del*xm(k)/xs(k)                                                  zipper
      dlx=max(xv(k)*del**2,dlx)                                            zipper
11921 continue                                                             zipper
11922 continue                                                             zipper
      if(nin.gt.nx)goto 11912                                              zipper
      if(dlx .ge. thr)goto 11961                                           zipper
      ixx=0                                                                zipper
11970 do 11971 j=1,ni                                                      zipper
      if(iy(j).eq.1)goto 11971                                             zipper
      if(ju(j).eq.0)goto 11971                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      g(j)=abs(dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(j))    zipper
      if(g(j) .le. ab*vp(j))goto 11991                                     zipper
      iy(j)=1                                                              zipper
      ixx=1                                                                zipper
11991 continue                                                             zipper
11971 continue                                                             zipper
11972 continue                                                             zipper
      if(ixx.eq.1) go to 10880                                             zipper
      goto 11912                                                           zipper
11961 continue                                                             zipper
      if(nlp .le. maxit)goto 12011                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
12011 continue                                                             zipper
10360 continue                                                             zipper
      iz=1                                                                 zipper
12020 continue                                                             zipper
12021 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
12030 do 12031 l=1,nin                                                     zipper
      k=ia(l)                                                              zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      gk=dot_product(y(jx(jb:je))+o,w(jx(jb:je))*x(jb:je))/xs(k)           zipper
      ak=a(k)                                                              zipper
      u=gk+ak*xv(k)                                                        zipper
      v=abs(u)-vp(k)*ab                                                    zipper
      a(k)=0.0                                                             zipper
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d   zipper
     *em)))
      if(a(k).eq.ak)goto 12031                                             zipper
      del=a(k)-ak                                                          zipper
      rsq=rsq+del*(2.0*gk-del*xv(k))                                       zipper
      y(jx(jb:je))=y(jx(jb:je))-del*x(jb:je)/xs(k)                         zipper
      o=o+del*xm(k)/xs(k)                                                  zipper
      dlx=max(xv(k)*del**2,dlx)                                            zipper
12031 continue                                                             zipper
12032 continue                                                             zipper
      if(dlx.lt.thr)goto 12022                                             zipper
      if(nlp .le. maxit)goto 12051                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
12051 continue                                                             zipper
      goto 12021                                                           zipper
12022 continue                                                             zipper
      jz=0                                                                 zipper
      goto 11911                                                           zipper
11912 continue                                                             zipper
      if(nin .le. nx)goto 12071                                            zipper
      jerr=-10000-m                                                        zipper
      goto 11832                                                           zipper
12071 continue                                                             zipper
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                zipper
      kin(m)=nin                                                           zipper
      rsqo(m)=rsq                                                          zipper
      almo(m)=alm                                                          zipper
      lmu=m                                                                zipper
      if(m.lt.mnl)goto 11831                                               zipper
      if(flmin.ge.1.0)goto 11831                                           zipper
      me=0                                                                 zipper
12080 do 12081 j=1,nin                                                     zipper
      if(ao(j,m).ne.0.0) me=me+1                                           zipper
12081 continue                                                             zipper
12082 continue                                                             zipper
      if(me.gt.ne)goto 11832                                               zipper
      if(rsq-rsq0.lt.sml*rsq)goto 11832                                    zipper
      if(rsq.gt.rsqmax)goto 11832                                          zipper
11831 continue                                                             zipper
11832 continue                                                             zipper
      deallocate(a,mm,g,iy)                                                zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spchkvars(no,ni,x,ix,ju)                                  zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*)                                                zipper
      integer ix(*),ju(ni)                                                 zipper
12090 do 12091 j=1,ni                                                      zipper
      ju(j)=0                                                              zipper
      jb=ix(j)                                                             zipper
      nj=ix(j+1)-jb                                                        zipper
      if(nj.eq.0)goto 12091                                                zipper
      je=ix(j+1)-1                                                         zipper
      if(nj .ge. no)goto 12111                                             zipper
12120 do 12121 i=jb,je                                                     zipper
      if(x(i).eq.0.0)goto 12121                                            zipper
      ju(j)=1                                                              zipper
      goto 12122                                                           zipper
12121 continue                                                             zipper
12122 continue                                                             zipper
      goto 12131                                                           zipper
12111 continue                                                             zipper
      t=x(jb)                                                              zipper
12140 do 12141 i=jb+1,je                                                   zipper
      if(x(i).eq.t)goto 12141                                              zipper
      ju(j)=1                                                              zipper
      goto 12142                                                           zipper
12141 continue                                                             zipper
12142 continue                                                             zipper
12131 continue                                                             zipper
12101 continue                                                             zipper
12091 continue                                                             zipper
12092 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine cmodval(a0,ca,ia,nin,x,ix,jx,n,f)                         zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision ca(*),x(*),f(n)                                     zipper
      integer ia(*),ix(*),jx(*)                                            zipper
      f=a0                                                                 zipper
12150 do 12151 j=1,nin                                                     zipper
      k=ia(j)                                                              zipper
      kb=ix(k)                                                             zipper
      ke=ix(k+1)-1                                                         zipper
      f(jx(kb:ke))=f(jx(kb:ke))+ca(j)*x(kb:ke)                             zipper
12151 continue                                                             zipper
12152 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      function row_prod(i,j,ia,ja,ra,w)                                    zipper
      implicit double precision(a-h,o-z)                                   zipper
      integer ia(*),ja(*)                                                  zipper
      double precision ra(*),w(*)                                          zipper
      row_prod=dot(ra(ia(i)),ra(ia(j)),ja(ia(i)),ja(ia(j)),  ia(i+1)-ia(   zipper
     *i),ia(j+1)-ia(j),w)
      return                                                               zipper
      end                                                                  zipper
      function dot(x,y,mx,my,nx,ny,w)                                      zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(*),w(*)                                      zipper
      integer mx(*),my(*)                                                  zipper
      i=1                                                                  zipper
      j=i                                                                  zipper
      s=0.0                                                                zipper
12160 continue                                                             zipper
12161 continue                                                             zipper
12170 continue                                                             zipper
12171 if(mx(i).ge.my(j))goto 12172                                         zipper
      i=i+1                                                                zipper
      if(i.gt.nx) go to 12180                                              zipper
      goto 12171                                                           zipper
12172 continue                                                             zipper
      if(mx(i).eq.my(j)) go to 12190                                       zipper
12200 continue                                                             zipper
12201 if(my(j).ge.mx(i))goto 12202                                         zipper
      j=j+1                                                                zipper
      if(j.gt.ny) go to 12180                                              zipper
      goto 12201                                                           zipper
12202 continue                                                             zipper
      if(mx(i).eq.my(j)) go to 12190                                       zipper
      goto 12161                                                           zipper
12190 continue                                                             zipper
      s=s+w(mx(i))*x(i)*y(j)                                               zipper
      i=i+1                                                                zipper
      if(i.gt.nx)goto 12162                                                zipper
      j=j+1                                                                zipper
      if(j.gt.ny)goto 12162                                                zipper
      goto 12161                                                           zipper
12162 continue                                                             zipper
12180 continue                                                             zipper
      dot=s                                                                zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine lognet (parm,no,ni,nc,x,y,g,jd,vp,cl,ne,nx,nlam,flmin,u   zipper
     *lam,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,je
     *rr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nla   zipper
     *m)
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   zipper
     *(2,ni)
      integer jd(*),ia(nx),nin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 12221                                    zipper
      jerr=10000                                                           zipper
      return                                                               zipper
12221 continue                                                             zipper
      allocate(ww(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(vq(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      if(kopt .ne. 2)goto 12241                                            zipper
      allocate(xv(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
12241 continue                                                             zipper
      if(isd .le. 0)goto 12261                                             zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
12261 continue                                                             zipper
      call chkvars(no,ni,x,ju)                                             zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 12281                                      zipper
      jerr=7777                                                            zipper
      return                                                               zipper
12281 continue                                                             zipper
      vq=max(0d0,vp)                                                       zipper
      vq=vq*ni/sum(vq)                                                     zipper
12290 do 12291 i=1,no                                                      zipper
      ww(i)=sum(y(i,:))                                                    zipper
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 zipper
12291 continue                                                             zipper
12292 continue                                                             zipper
      sw=sum(ww)                                                           zipper
      ww=ww/sw                                                             zipper
      if(nc .ne. 1)goto 12311                                              zipper
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        zipper
      if(isd .le. 0)goto 12331                                             zipper
12340 do 12341 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
12341 continue                                                             zipper
12342 continue                                                             zipper
12331 continue                                                             zipper
      call lognet2n(parm,no,ni,x,y(:,1),g(:,1),ww,ju,vq,cl,ne,nx,nlam,fl   zipper
     *min,ulam,  thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,n
     *lp,jerr)
      goto 12301                                                           zipper
12311 if(kopt .ne. 2)goto 12351                                            zipper
      call multlstandard1(no,ni,x,ww,ju,isd,intr,xm,xs,xv)                 zipper
      if(isd .le. 0)goto 12371                                             zipper
12380 do 12381 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
12381 continue                                                             zipper
12382 continue                                                             zipper
12371 continue                                                             zipper
      call multlognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,   zipper
     *ulam,thr,  intr,maxit,xv,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      goto 12391                                                           zipper
12351 continue                                                             zipper
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        zipper
      if(isd .le. 0)goto 12411                                             zipper
12420 do 12421 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
12421 continue                                                             zipper
12422 continue                                                             zipper
12411 continue                                                             zipper
      call lognetn(parm,no,ni,nc,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam   zipper
     *,thr,  isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
12391 continue                                                             zipper
12301 continue                                                             zipper
      if(jerr.gt.0) return                                                 zipper
      dev0=2.0*sw*dev0                                                     zipper
12430 do 12431 k=1,lmu                                                     zipper
      nk=nin(k)                                                            zipper
12440 do 12441 ic=1,nc                                                     zipper
      if(isd .le. 0)goto 12461                                             zipper
12470 do 12471 l=1,nk                                                      zipper
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      zipper
12471 continue                                                             zipper
12472 continue                                                             zipper
12461 continue                                                             zipper
      if(intr .ne. 0)goto 12491                                            zipper
      a0(ic,k)=0.0                                                         zipper
      goto 12501                                                           zipper
12491 continue                                                             zipper
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            zipper
12501 continue                                                             zipper
12481 continue                                                             zipper
12441 continue                                                             zipper
12442 continue                                                             zipper
12431 continue                                                             zipper
12432 continue                                                             zipper
      deallocate(ww,ju,vq,xm)                                              zipper
      if(isd.gt.0) deallocate(xs)                                          zipper
      if(kopt.eq.2) deallocate(xv)                                         zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine lstandard1 (no,ni,x,w,ju,isd,intr,xm,xs)                  zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),w(no),xm(ni),xs(ni)                        zipper
      integer ju(ni)                                                       zipper
      if(intr .ne. 0)goto 12521                                            zipper
12530 do 12531 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 12531                                             zipper
      xm(j)=0.0                                                            zipper
      if(isd .eq. 0)goto 12551                                             zipper
      vc=dot_product(w,x(:,j)**2)-dot_product(w,x(:,j))**2                 zipper
      xs(j)=sqrt(vc)                                                       zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
12551 continue                                                             zipper
12531 continue                                                             zipper
12532 continue                                                             zipper
      return                                                               zipper
12521 continue                                                             zipper
12560 do 12561 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 12561                                             zipper
      xm(j)=dot_product(w,x(:,j))                                          zipper
      x(:,j)=x(:,j)-xm(j)                                                  zipper
      if(isd .le. 0)goto 12581                                             zipper
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
12581 continue                                                             zipper
12561 continue                                                             zipper
12562 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multlstandard1 (no,ni,x,w,ju,isd,intr,xm,xs,xv)           zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),w(no),xm(ni),xs(ni),xv(ni)                 zipper
      integer ju(ni)                                                       zipper
      if(intr .ne. 0)goto 12601                                            zipper
12610 do 12611 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 12611                                             zipper
      xm(j)=0.0                                                            zipper
      xv(j)=dot_product(w,x(:,j)**2)                                       zipper
      if(isd .eq. 0)goto 12631                                             zipper
      xbq=dot_product(w,x(:,j))**2                                         zipper
      vc=xv(j)-xbq                                                         zipper
      xs(j)=sqrt(vc)                                                       zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
      xv(j)=1.0+xbq/vc                                                     zipper
12631 continue                                                             zipper
12611 continue                                                             zipper
12612 continue                                                             zipper
      return                                                               zipper
12601 continue                                                             zipper
12640 do 12641 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 12641                                             zipper
      xm(j)=dot_product(w,x(:,j))                                          zipper
      x(:,j)=x(:,j)-xm(j)                                                  zipper
      xv(j)=dot_product(w,x(:,j)**2)                                       zipper
      if(isd .le. 0)goto 12661                                             zipper
      xs(j)=sqrt(xv(j))                                                    zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
      xv(j)=1.0                                                            zipper
12661 continue                                                             zipper
12641 continue                                                             zipper
12642 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine lognet2n(parm,no,ni,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin,u   zipper
     *lam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2   zipper
     *,ni)
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             zipper
      integer ju(ni),m(nx),kin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: b,bs,v,r,xv,q,ga          
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      allocate(b(0:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xv(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(bs(0:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ixx(1:ni),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(r(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(v(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(q(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      fmax=log(1.0/pmin-1.0)                                               zipper
      fmin=-fmax                                                           zipper
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      zipper
      bta=parm                                                             zipper
      omb=1.0-bta                                                          zipper
      q0=dot_product(w,y)                                                  zipper
      if(q0 .gt. pmin)goto 12681                                           zipper
      jerr=8001                                                            zipper
      return                                                               zipper
12681 continue                                                             zipper
      if(q0 .lt. 1.0-pmin)goto 12701                                       zipper
      jerr=9001                                                            zipper
      return                                                               zipper
12701 continue                                                             zipper
      if(intr.eq.0.0) q0=0.5                                               zipper
      ixx=0                                                                zipper
      al=0.0                                                               zipper
      bz=0.0                                                               zipper
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    zipper
      if(nonzero(no,g) .ne. 0)goto 12721                                   zipper
      vi=q0*(1.0-q0)                                                       zipper
      b(0)=bz                                                              zipper
      v=vi*w                                                               zipper
      r=w*(y-q0)                                                           zipper
      q=q0                                                                 zipper
      xmz=vi                                                               zipper
      dev1=-(bz*q0+log(1.0-q0))                                            zipper
      goto 12731                                                           zipper
12721 continue                                                             zipper
      b(0)=0.0                                                             zipper
      if(intr .eq. 0)goto 12751                                            zipper
      b(0)=azero(no,y,g,w,jerr)                                            zipper
      if(jerr.ne.0) return                                                 zipper
12751 continue                                                             zipper
      q=1.0/(1.0+exp(-b(0)-g))                                             zipper
      v=w*q*(1.0-q)                                                        zipper
      r=w*(y-q)                                                            zipper
      xmz=sum(v)                                                           zipper
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        zipper
12731 continue                                                             zipper
12711 continue                                                             zipper
      if(kopt .le. 0)goto 12771                                            zipper
      if(isd .le. 0 .or. intr .eq. 0)goto 12791                            zipper
      xv=0.25                                                              zipper
      goto 12801                                                           zipper
12791 continue                                                             zipper
12810 do 12811 j=1,ni                                                      zipper
      if(ju(j).ne.0) xv(j)=0.25*dot_product(w,x(:,j)**2)                   zipper
12811 continue                                                             zipper
12812 continue                                                             zipper
12801 continue                                                             zipper
12781 continue                                                             zipper
12771 continue                                                             zipper
      dev0=dev1                                                            zipper
12820 do 12821 i=1,no                                                      zipper
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        zipper
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              zipper
12821 continue                                                             zipper
12822 continue                                                             zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 12841                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
12841 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      mnl=min(mnlam,nlam)                                                  zipper
      bs=0.0                                                               zipper
      b(1:ni)=0.0                                                          zipper
      shr=shri*dev0                                                        zipper
12850 do 12851 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 12851                                             zipper
      ga(j)=abs(dot_product(r,x(:,j)))                                     zipper
12851 continue                                                             zipper
12852 continue                                                             zipper
12860 do 12861 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 12881                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 12871                                                           zipper
12881 if(ilm .le. 2)goto 12891                                             zipper
      al=al*alf                                                            zipper
      goto 12871                                                           zipper
12891 if(ilm .ne. 1)goto 12901                                             zipper
      al=big                                                               zipper
      goto 12911                                                           zipper
12901 continue                                                             zipper
      al0=0.0                                                              zipper
12920 do 12921 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 12921                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
12921 continue                                                             zipper
12922 continue                                                             zipper
      al0=al0/max(bta,1.0d-3)                                              zipper
      al=alf*al0                                                           zipper
12911 continue                                                             zipper
12871 continue                                                             zipper
      al2=al*omb                                                           zipper
      al1=al*bta                                                           zipper
      tlam=bta*(2.0*al-al0)                                                zipper
12930 do 12931 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 12931                                            zipper
      if(ju(k).eq.0)goto 12931                                             zipper
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     zipper
12931 continue                                                             zipper
12932 continue                                                             zipper
10880 continue                                                             zipper
12940 continue                                                             zipper
12941 continue                                                             zipper
      bs(0)=b(0)                                                           zipper
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                zipper
      if(kopt .ne. 0)goto 12961                                            zipper
12970 do 12971 j=1,ni                                                      zipper
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)                       zipper
12971 continue                                                             zipper
12972 continue                                                             zipper
12961 continue                                                             zipper
12980 continue                                                             zipper
12981 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
12990 do 12991 k=1,ni                                                      zipper
      if(ixx(k).eq.0)goto 12991                                            zipper
      bk=b(k)                                                              zipper
      gk=dot_product(r,x(:,k))                                             zipper
      u=gk+xv(k)*b(k)                                                      zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 13011                                            zipper
      b(k)=0.0                                                             zipper
      goto 13021                                                           zipper
13011 continue                                                             zipper
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          zipper
13021 continue                                                             zipper
13001 continue                                                             zipper
      d=b(k)-bk                                                            zipper
      if(abs(d).le.0.0)goto 12991                                          zipper
      dlx=max(dlx,xv(k)*d**2)                                              zipper
      r=r-d*v*x(:,k)                                                       zipper
      if(mm(k) .ne. 0)goto 13041                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 12992                                              zipper
      mm(k)=nin                                                            zipper
      m(nin)=k                                                             zipper
13041 continue                                                             zipper
12991 continue                                                             zipper
12992 continue                                                             zipper
      if(nin.gt.nx)goto 12982                                              zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=sum(r)/xmz                                           zipper
      if(d .eq. 0.0)goto 13061                                             zipper
      b(0)=b(0)+d                                                          zipper
      dlx=max(dlx,xmz*d**2)                                                zipper
      r=r-d*v                                                              zipper
13061 continue                                                             zipper
      if(dlx.lt.shr)goto 12982                                             zipper
      if(nlp .le. maxit)goto 13081                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
13081 continue                                                             zipper
13090 continue                                                             zipper
13091 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
13100 do 13101 l=1,nin                                                     zipper
      k=m(l)                                                               zipper
      bk=b(k)                                                              zipper
      gk=dot_product(r,x(:,k))                                             zipper
      u=gk+xv(k)*b(k)                                                      zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 13121                                            zipper
      b(k)=0.0                                                             zipper
      goto 13131                                                           zipper
13121 continue                                                             zipper
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          zipper
13131 continue                                                             zipper
13111 continue                                                             zipper
      d=b(k)-bk                                                            zipper
      if(abs(d).le.0.0)goto 13101                                          zipper
      dlx=max(dlx,xv(k)*d**2)                                              zipper
      r=r-d*v*x(:,k)                                                       zipper
13101 continue                                                             zipper
13102 continue                                                             zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=sum(r)/xmz                                           zipper
      if(d .eq. 0.0)goto 13151                                             zipper
      b(0)=b(0)+d                                                          zipper
      dlx=max(dlx,xmz*d**2)                                                zipper
      r=r-d*v                                                              zipper
13151 continue                                                             zipper
      if(dlx.lt.shr)goto 13092                                             zipper
      if(nlp .le. maxit)goto 13171                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
13171 continue                                                             zipper
      goto 13091                                                           zipper
13092 continue                                                             zipper
      goto 12981                                                           zipper
12982 continue                                                             zipper
      if(nin.gt.nx)goto 12942                                              zipper
13180 do 13181 i=1,no                                                      zipper
      fi=b(0)+g(i)                                                         zipper
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))            zipper
      if(fi .ge. fmin)goto 13201                                           zipper
      q(i)=0.0                                                             zipper
      goto 13191                                                           zipper
13201 if(fi .le. fmax)goto 13211                                           zipper
      q(i)=1.0                                                             zipper
      goto 13221                                                           zipper
13211 continue                                                             zipper
      q(i)=1.0/(1.0+exp(-fi))                                              zipper
13221 continue                                                             zipper
13191 continue                                                             zipper
13181 continue                                                             zipper
13182 continue                                                             zipper
      v=w*q*(1.0-q)                                                        zipper
      xmz=sum(v)                                                           zipper
      if(xmz.le.vmin)goto 12942                                            zipper
      r=w*(y-q)                                                            zipper
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 13241                           zipper
      ix=0                                                                 zipper
13250 do 13251 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 13251                           zipper
      ix=1                                                                 zipper
      goto 13252                                                           zipper
13251 continue                                                             zipper
13252 continue                                                             zipper
      if(ix .ne. 0)goto 13271                                              zipper
13280 do 13281 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 13281                                            zipper
      if(ju(k).eq.0)goto 13281                                             zipper
      ga(k)=abs(dot_product(r,x(:,k)))                                     zipper
      if(ga(k) .le. al1*vp(k))goto 13301                                   zipper
      ixx(k)=1                                                             zipper
      ix=1                                                                 zipper
13301 continue                                                             zipper
13281 continue                                                             zipper
13282 continue                                                             zipper
      if(ix.eq.1) go to 10880                                              zipper
      goto 12942                                                           zipper
13271 continue                                                             zipper
13241 continue                                                             zipper
      goto 12941                                                           zipper
12942 continue                                                             zipper
      if(nin .le. nx)goto 13321                                            zipper
      jerr=-10000-ilm                                                      zipper
      goto 12862                                                           zipper
13321 continue                                                             zipper
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                zipper
      kin(ilm)=nin                                                         zipper
      a0(ilm)=b(0)                                                         zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      devi=dev2(no,w,y,q,pmin)                                             zipper
      dev(ilm)=(dev1-devi)/dev0                                            zipper
      if(xmz.le.vmin)goto 12862                                            zipper
      if(ilm.lt.mnl)goto 12861                                             zipper
      if(flmin.ge.1.0)goto 12861                                           zipper
      me=0                                                                 zipper
13330 do 13331 j=1,nin                                                     zipper
      if(a(j,ilm).ne.0.0) me=me+1                                          zipper
13331 continue                                                             zipper
13332 continue                                                             zipper
      if(me.gt.ne)goto 12862                                               zipper
      if(dev(ilm).gt.devmax)goto 12862                                     zipper
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12862                             zipper
12861 continue                                                             zipper
12862 continue                                                             zipper
      g=log(q/(1.0-q))                                                     zipper
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)                                  zipper
      return                                                               zipper
      end                                                                  zipper
      function dev2(n,w,y,p,pmin)                                          zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision w(n),y(n),p(n)                                      zipper
      pmax=1.0-pmin                                                        zipper
      s=0.0                                                                zipper
13340 do 13341 i=1,n                                                       zipper
      pi=min(max(pmin,p(i)),pmax)                                          zipper
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                       zipper
13341 continue                                                             zipper
13342 continue                                                             zipper
      dev2=s                                                               zipper
      return                                                               zipper
      end                                                                  zipper
      function azero(n,y,g,q,jerr)                                         zipper
      implicit double precision(a-h,o-z)                                   zipper
      parameter(eps=1.0d-7)                                                zipper
      double precision y(n),g(n),q(n)                                      zipper
      double precision, dimension (:), allocatable :: e,p,w                     
      azero = 0.0                                                          zipper
      allocate(e(1:n),stat=jerr)                                           zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(p(1:n),stat=jerr)                                           zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(w(1:n),stat=jerr)                                           zipper
      if(jerr.ne.0) return                                                 zipper
      az=0.0                                                               zipper
      e=exp(-g)                                                            zipper
      qy=dot_product(q,y)                                                  zipper
      p=1.0/(1.0+e)                                                        zipper
13350 continue                                                             zipper
13351 continue                                                             zipper
      w=q*p*(1.0-p)                                                        zipper
      d=(qy-dot_product(q,p))/sum(w)                                       zipper
      az=az+d                                                              zipper
      if(abs(d).lt.eps)goto 13352                                          zipper
      ea0=exp(-az)                                                         zipper
      p=1.0/(1.0+ea0*e)                                                    zipper
      goto 13351                                                           zipper
13352 continue                                                             zipper
      azero=az                                                             zipper
      deallocate(e,p,w)                                                    zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine lognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,flmin   zipper
     *,ulam,shri,  isd,intr,maxit,kopt,lmu,a0,a,m,kin,dev0,dev,alm,nlp,j
     *err)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   zipper
     *)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   zipper
     *2,ni)
      integer ju(ni),m(nx),kin(nlam)                                       zipper
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp,sxpl                  
      double precision, dimension (:), allocatable :: di,v,r,ga                 
      double precision, dimension (:,:), allocatable :: b,bs,xv                 
      integer, dimension (:), allocatable :: mm,is,ixx                          
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(xv(1:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return		                                                    
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      exmn=-exmx                                                           zipper
      allocate(r(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(v(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(is(1:max(nc,ni)),stat=jerr)                                 zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sxp(1:no),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sxpl(1:no),stat=jerr)                                       zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(di(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ixx(1:ni),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      pmax=1.0-pmin                                                        zipper
      emin=pmin/pmax                                                       zipper
      emax=1.0/emin                                                        zipper
      pfm=(1.0+pmin)*pmin                                                  zipper
      pfx=(1.0-pmin)*pmax                                                  zipper
      vmin=pfm*pmax                                                        zipper
      bta=parm                                                             zipper
      omb=1.0-bta                                                          zipper
      dev1=0.0                                                             zipper
      dev0=0.0                                                             zipper
13360 do 13361 ic=1,nc                                                     zipper
      q0=dot_product(w,y(:,ic))                                            zipper
      if(q0 .gt. pmin)goto 13381                                           zipper
      jerr =8000+ic                                                        zipper
      return                                                               zipper
13381 continue                                                             zipper
      if(q0 .lt. 1.0-pmin)goto 13401                                       zipper
      jerr =9000+ic                                                        zipper
      return                                                               zipper
13401 continue                                                             zipper
      if(intr .ne. 0)goto 13421                                            zipper
      q0=1.0/nc                                                            zipper
      b(0,ic)=0.0                                                          zipper
      goto 13431                                                           zipper
13421 continue                                                             zipper
      b(0,ic)=log(q0)                                                      zipper
      dev1=dev1-q0*b(0,ic)                                                 zipper
13431 continue                                                             zipper
13411 continue                                                             zipper
      b(1:ni,ic)=0.0                                                       zipper
13361 continue                                                             zipper
13362 continue                                                             zipper
      if(intr.eq.0) dev1=log(float(nc))                                    zipper
      ixx=0                                                                zipper
      al=0.0                                                               zipper
      if(nonzero(no*nc,g) .ne. 0)goto 13451                                zipper
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         zipper
      sxp=0.0                                                              zipper
13460 do 13461 ic=1,nc                                                     zipper
      q(:,ic)=exp(b(0,ic))                                                 zipper
      sxp=sxp+q(:,ic)                                                      zipper
13461 continue                                                             zipper
13462 continue                                                             zipper
      goto 13471                                                           zipper
13451 continue                                                             zipper
13480 do 13481 i=1,no                                                      zipper
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         zipper
13481 continue                                                             zipper
13482 continue                                                             zipper
      sxp=0.0                                                              zipper
      if(intr .ne. 0)goto 13501                                            zipper
      b(0,:)=0.0                                                           zipper
      goto 13511                                                           zipper
13501 continue                                                             zipper
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 zipper
      if(jerr.ne.0) return                                                 zipper
13511 continue                                                             zipper
13491 continue                                                             zipper
      dev1=0.0                                                             zipper
13520 do 13521 ic=1,nc                                                     zipper
      q(:,ic)=b(0,ic)+g(:,ic)                                              zipper
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             zipper
      q(:,ic)=exp(q(:,ic))                                                 zipper
      sxp=sxp+q(:,ic)                                                      zipper
13521 continue                                                             zipper
13522 continue                                                             zipper
      sxpl=w*log(sxp)                                                      zipper
13530 do 13531 ic=1,nc                                                     zipper
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  zipper
13531 continue                                                             zipper
13532 continue                                                             zipper
13471 continue                                                             zipper
13441 continue                                                             zipper
13540 do 13541 ic=1,nc                                                     zipper
13550 do 13551 i=1,no                                                      zipper
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               zipper
13551 continue                                                             zipper
13552 continue                                                             zipper
13541 continue                                                             zipper
13542 continue                                                             zipper
      dev0=dev0+dev1                                                       zipper
      if(kopt .le. 0)goto 13571                                            zipper
      if(isd .le. 0 .or. intr .eq. 0)goto 13591                            zipper
      xv=0.25                                                              zipper
      goto 13601                                                           zipper
13591 continue                                                             zipper
13610 do 13611 j=1,ni                                                      zipper
      if(ju(j).ne.0) xv(j,:)=0.25*dot_product(w,x(:,j)**2)                 zipper
13611 continue                                                             zipper
13612 continue                                                             zipper
13601 continue                                                             zipper
13581 continue                                                             zipper
13571 continue                                                             zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 13631                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
13631 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nin=0                                                                zipper
      nlp=0                                                                zipper
      mnl=min(mnlam,nlam)                                                  zipper
      bs=0.0                                                               zipper
      shr=shri*dev0                                                        zipper
      ga=0.0                                                               zipper
13640 do 13641 ic=1,nc                                                     zipper
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            zipper
13650 do 13651 j=1,ni                                                      zipper
      if(ju(j).ne.0) ga(j)=max(ga(j),abs(dot_product(r,x(:,j))))           zipper
13651 continue                                                             zipper
13652 continue                                                             zipper
13641 continue                                                             zipper
13642 continue                                                             zipper
13660 do 13661 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 13681                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 13671                                                           zipper
13681 if(ilm .le. 2)goto 13691                                             zipper
      al=al*alf                                                            zipper
      goto 13671                                                           zipper
13691 if(ilm .ne. 1)goto 13701                                             zipper
      al=big                                                               zipper
      goto 13711                                                           zipper
13701 continue                                                             zipper
      al0=0.0                                                              zipper
13720 do 13721 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 13721                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
13721 continue                                                             zipper
13722 continue                                                             zipper
      al0=al0/max(bta,1.0d-3)                                              zipper
      al=alf*al0                                                           zipper
13711 continue                                                             zipper
13671 continue                                                             zipper
      al2=al*omb                                                           zipper
      al1=al*bta                                                           zipper
      tlam=bta*(2.0*al-al0)                                                zipper
13730 do 13731 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 13731                                            zipper
      if(ju(k).eq.0)goto 13731                                             zipper
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     zipper
13731 continue                                                             zipper
13732 continue                                                             zipper
10880 continue                                                             zipper
13740 continue                                                             zipper
13741 continue                                                             zipper
      ix=0                                                                 zipper
      jx=ix                                                                zipper
      ig=0                                                                 zipper
13750 do 13751 ic=1,nc                                                     zipper
      bs(0,ic)=b(0,ic)                                                     zipper
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          zipper
      xmz=0.0                                                              zipper
13760 do 13761 i=1,no                                                      zipper
      pic=q(i,ic)/sxp(i)                                                   zipper
      if(pic .ge. pfm)goto 13781                                           zipper
      pic=0.0                                                              zipper
      v(i)=0.0                                                             zipper
      goto 13771                                                           zipper
13781 if(pic .le. pfx)goto 13791                                           zipper
      pic=1.0                                                              zipper
      v(i)=0.0                                                             zipper
      goto 13801                                                           zipper
13791 continue                                                             zipper
      v(i)=w(i)*pic*(1.0-pic)                                              zipper
      xmz=xmz+v(i)                                                         zipper
13801 continue                                                             zipper
13771 continue                                                             zipper
      r(i)=w(i)*(y(i,ic)-pic)                                              zipper
13761 continue                                                             zipper
13762 continue                                                             zipper
      if(xmz.le.vmin)goto 13751                                            zipper
      ig=1                                                                 zipper
      if(kopt .ne. 0)goto 13821                                            zipper
13830 do 13831 j=1,ni                                                      zipper
      if(ixx(j).gt.0) xv(j,ic)=dot_product(v,x(:,j)**2)                    zipper
13831 continue                                                             zipper
13832 continue                                                             zipper
13821 continue                                                             zipper
13840 continue                                                             zipper
13841 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
13850 do 13851 k=1,ni                                                      zipper
      if(ixx(k).eq.0)goto 13851                                            zipper
      bk=b(k,ic)                                                           zipper
      gk=dot_product(r,x(:,k))                                             zipper
      u=gk+xv(k,ic)*b(k,ic)                                                zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 13871                                            zipper
      b(k,ic)=0.0                                                          zipper
      goto 13881                                                           zipper
13871 continue                                                             zipper
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   zipper
     *)
13881 continue                                                             zipper
13861 continue                                                             zipper
      d=b(k,ic)-bk                                                         zipper
      if(abs(d).le.0.0)goto 13851                                          zipper
      dlx=max(dlx,xv(k,ic)*d**2)                                           zipper
      r=r-d*v*x(:,k)                                                       zipper
      if(mm(k) .ne. 0)goto 13901                                           zipper
      nin=nin+1                                                            zipper
      if(nin .le. nx)goto 13921                                            zipper
      jx=1                                                                 zipper
      goto 13852                                                           zipper
13921 continue                                                             zipper
      mm(k)=nin                                                            zipper
      m(nin)=k                                                             zipper
13901 continue                                                             zipper
13851 continue                                                             zipper
13852 continue                                                             zipper
      if(jx.gt.0)goto 13842                                                zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=sum(r)/xmz                                           zipper
      if(d .eq. 0.0)goto 13941                                             zipper
      b(0,ic)=b(0,ic)+d                                                    zipper
      dlx=max(dlx,xmz*d**2)                                                zipper
      r=r-d*v                                                              zipper
13941 continue                                                             zipper
      if(dlx.lt.shr)goto 13842                                             zipper
      if(nlp .le. maxit)goto 13961                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
13961 continue                                                             zipper
13970 continue                                                             zipper
13971 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
13980 do 13981 l=1,nin                                                     zipper
      k=m(l)                                                               zipper
      bk=b(k,ic)                                                           zipper
      gk=dot_product(r,x(:,k))                                             zipper
      u=gk+xv(k,ic)*b(k,ic)                                                zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 14001                                            zipper
      b(k,ic)=0.0                                                          zipper
      goto 14011                                                           zipper
14001 continue                                                             zipper
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   zipper
     *)
14011 continue                                                             zipper
13991 continue                                                             zipper
      d=b(k,ic)-bk                                                         zipper
      if(abs(d).le.0.0)goto 13981                                          zipper
      dlx=max(dlx,xv(k,ic)*d**2)                                           zipper
      r=r-d*v*x(:,k)                                                       zipper
13981 continue                                                             zipper
13982 continue                                                             zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=sum(r)/xmz                                           zipper
      if(d .eq. 0.0)goto 14031                                             zipper
      b(0,ic)=b(0,ic)+d                                                    zipper
      dlx=max(dlx,xmz*d**2)                                                zipper
      r=r-d*v                                                              zipper
14031 continue                                                             zipper
      if(dlx.lt.shr)goto 13972                                             zipper
      if(nlp .le. maxit)goto 14051                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
14051 continue                                                             zipper
      goto 13971                                                           zipper
13972 continue                                                             zipper
      goto 13841                                                           zipper
13842 continue                                                             zipper
      if(jx.gt.0)goto 13752                                                zipper
      if(xmz*(b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                            zipper
      if(ix .ne. 0)goto 14071                                              zipper
14080 do 14081 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 14101                zipper
      ix=1                                                                 zipper
      goto 14082                                                           zipper
14101 continue                                                             zipper
14081 continue                                                             zipper
14082 continue                                                             zipper
14071 continue                                                             zipper
14110 do 14111 i=1,no                                                      zipper
      fi=b(0,ic)+g(i,ic)                                                   zipper
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         zipper
      fi=min(max(exmn,fi),exmx)                                            zipper
      sxp(i)=sxp(i)-q(i,ic)                                                zipper
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    zipper
      sxp(i)=sxp(i)+q(i,ic)                                                zipper
14111 continue                                                             zipper
14112 continue                                                             zipper
13751 continue                                                             zipper
13752 continue                                                             zipper
      s=-sum(b(0,:))/nc                                                    zipper
      b(0,:)=b(0,:)+s                                                      zipper
      di=s                                                                 zipper
14120 do 14121 j=1,nin                                                     zipper
      l=m(j)                                                               zipper
      if(vp(l) .gt. 0.0)goto 14141                                         zipper
      s=sum(b(l,:))/nc                                                     zipper
      goto 14151                                                           zipper
14141 continue                                                             zipper
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     zipper
14151 continue                                                             zipper
14131 continue                                                             zipper
      b(l,:)=b(l,:)-s                                                      zipper
      di=di-s*x(:,l)                                                       zipper
14121 continue                                                             zipper
14122 continue                                                             zipper
      di=exp(di)                                                           zipper
      sxp=sxp*di                                                           zipper
14160 do 14161 ic=1,nc                                                     zipper
      q(:,ic)=q(:,ic)*di                                                   zipper
14161 continue                                                             zipper
14162 continue                                                             zipper
      if(jx.gt.0)goto 13742                                                zipper
      if(ig.eq.0)goto 13742                                                zipper
      if(ix .ne. 0)goto 14181                                              zipper
14190 do 14191 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 14191                                            zipper
      if(ju(k).eq.0)goto 14191                                             zipper
      ga(k)=0.0                                                            zipper
14191 continue                                                             zipper
14192 continue                                                             zipper
14200 do 14201 ic=1,nc                                                     zipper
      r=w*(y(:,ic)-q(:,ic)/sxp)                                            zipper
14210 do 14211 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 14211                                            zipper
      if(ju(k).eq.0)goto 14211                                             zipper
      ga(k)=max(ga(k),abs(dot_product(r,x(:,k))))                          zipper
14211 continue                                                             zipper
14212 continue                                                             zipper
14201 continue                                                             zipper
14202 continue                                                             zipper
14220 do 14221 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 14221                                            zipper
      if(ju(k).eq.0)goto 14221                                             zipper
      if(ga(k) .le. al1*vp(k))goto 14241                                   zipper
      ixx(k)=1                                                             zipper
      ix=1                                                                 zipper
14241 continue                                                             zipper
14221 continue                                                             zipper
14222 continue                                                             zipper
      if(ix.eq.1) go to 10880                                              zipper
      goto 13742                                                           zipper
14181 continue                                                             zipper
      goto 13741                                                           zipper
13742 continue                                                             zipper
      if(jx .le. 0)goto 14261                                              zipper
      jerr=-10000-ilm                                                      zipper
      goto 13662                                                           zipper
14261 continue                                                             zipper
      devi=0.0                                                             zipper
14270 do 14271 ic=1,nc                                                     zipper
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          zipper
      a0(ic,ilm)=b(0,ic)                                                   zipper
14280 do 14281 i=1,no                                                      zipper
      if(y(i,ic).le.0.0)goto 14281                                         zipper
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           zipper
14281 continue                                                             zipper
14282 continue                                                             zipper
14271 continue                                                             zipper
14272 continue                                                             zipper
      kin(ilm)=nin                                                         zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      dev(ilm)=(dev1-devi)/dev0                                            zipper
      if(ig.eq.0)goto 13662                                                zipper
      if(ilm.lt.mnl)goto 13661                                             zipper
      if(flmin.ge.1.0)goto 13661                                           zipper
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 13662             zipper
      if(dev(ilm).gt.devmax)goto 13662                                     zipper
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 13662                             zipper
13661 continue                                                             zipper
13662 continue                                                             zipper
      g=log(q)                                                             zipper
14290 do 14291 i=1,no                                                      zipper
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         zipper
14291 continue                                                             zipper
14292 continue                                                             zipper
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,ga,ixx)                           zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine kazero(kk,n,y,g,q,az,jerr)                                zipper
      implicit double precision(a-h,o-z)                                   zipper
      parameter(eps=1.0d-7)                                                zipper
      double precision y(n,kk),g(n,kk),q(n),az(kk)                         zipper
      double precision, dimension (:), allocatable :: s                         
      double precision, dimension (:,:), allocatable :: e                       
      allocate(e(1:n,1:kk),stat=jerr)                                           
      if(jerr.ne.0) return                                                      
      allocate(s(1:n),stat=jerr)                                           zipper
      if(jerr.ne.0) return                                                 zipper
      az=0.0                                                               zipper
      e=exp(g)                                                             zipper
14300 do 14301 i=1,n                                                       zipper
      s(i)=sum(e(i,:))                                                     zipper
14301 continue                                                             zipper
14302 continue                                                             zipper
14310 continue                                                             zipper
14311 continue                                                             zipper
      dm=0.0                                                               zipper
14320 do 14321 k=1,kk                                                      zipper
      t=0.0                                                                zipper
      u=t                                                                  zipper
14330 do 14331 i=1,n                                                       zipper
      pik=e(i,k)/s(i)                                                      zipper
      t=t+q(i)*(y(i,k)-pik)                                                zipper
      u=u+q(i)*pik*(1.0-pik)                                               zipper
14331 continue                                                             zipper
14332 continue                                                             zipper
      d=t/u                                                                zipper
      az(k)=az(k)+d                                                        zipper
      ed=exp(d)                                                            zipper
      dm=max(dm,abs(d))                                                    zipper
14340 do 14341 i=1,n                                                       zipper
      z=e(i,k)                                                             zipper
      e(i,k)=z*ed                                                          zipper
      s(i)=s(i)-z+e(i,k)                                                   zipper
14341 continue                                                             zipper
14342 continue                                                             zipper
14321 continue                                                             zipper
14322 continue                                                             zipper
      if(dm.lt.eps)goto 14312                                              zipper
      goto 14311                                                           zipper
14312 continue                                                             zipper
      az=az-sum(az)/kk                                                     zipper
      deallocate(e,s)                                                      zipper
      return                                                               zipper
      end                                                                  zipper
      function elc(parm,n,cl,a,m)                                          zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision a(n),cl(2)                                          zipper
      integer m(n)                                                         zipper
      fn=n                                                                 zipper
      am=sum(a)/fn                                                         zipper
      if((parm .ne. 0.0) .and. (n .ne. 2))goto 14361                       zipper
      elc=am                                                               zipper
      go to 14370                                                          zipper
14361 continue                                                             zipper
14380 do 14381 i=1,n                                                       zipper
      m(i)=i                                                               zipper
14381 continue                                                             zipper
14382 continue                                                             zipper
      call psort7(a,m,1,n)                                                 zipper
      if(a(m(1)) .ne. a(m(n)))goto 14401                                   zipper
      elc=a(1)                                                             zipper
      go to 14370                                                          zipper
14401 continue                                                             zipper
      if(mod(n,2) .ne. 1)goto 14421                                        zipper
      ad=a(m(n/2+1))                                                       zipper
      goto 14431                                                           zipper
14421 continue                                                             zipper
      ad=0.5*(a(m(n/2+1))+a(m(n/2)))                                       zipper
14431 continue                                                             zipper
14411 continue                                                             zipper
      if(parm .ne. 1.0)goto 14451                                          zipper
      elc=ad                                                               zipper
      go to 14370                                                          zipper
14451 continue                                                             zipper
      b1=min(am,ad)                                                        zipper
      b2=max(am,ad)                                                        zipper
      k2=1                                                                 zipper
14460 continue                                                             zipper
14461 if(a(m(k2)).gt.b1)goto 14462                                         zipper
      k2=k2+1                                                              zipper
      goto 14461                                                           zipper
14462 continue                                                             zipper
      k1=k2-1                                                              zipper
14470 continue                                                             zipper
14471 if(a(m(k2)).ge.b2)goto 14472                                         zipper
      k2=k2+1                                                              zipper
      goto 14471                                                           zipper
14472 continue                                                             zipper
      r=parm/((1.0-parm)*fn)                                               zipper
      is=0                                                                 zipper
      sm=n-2*(k1-1)                                                        zipper
14480 do 14481 k=k1,k2-1                                                   zipper
      sm=sm-2.0                                                            zipper
      s=r*sm+am                                                            zipper
      if(s .le. a(m(k)) .or. s .gt. a(m(k+1)))goto 14501                   zipper
      is=k                                                                 zipper
      goto 14482                                                           zipper
14501 continue                                                             zipper
14481 continue                                                             zipper
14482 continue                                                             zipper
      if(is .eq. 0)goto 14521                                              zipper
      elc=s                                                                zipper
      go to 14370                                                          zipper
14521 continue                                                             zipper
      r2=2.0*r                                                             zipper
      s1=a(m(k1))                                                          zipper
      am2=2.0*am                                                           zipper
      cri=r2*sum(abs(a-s1))+s1*(s1-am2)                                    zipper
      elc=s1                                                               zipper
14530 do 14531 k=k1+1,k2                                                   zipper
      s=a(m(k))                                                            zipper
      if(s.eq.s1)goto 14531                                                zipper
      c=r2*sum(abs(a-s))+s*(s-am2)                                         zipper
      if(c .ge. cri)goto 14551                                             zipper
      cri=c                                                                zipper
      elc=s                                                                zipper
14551 continue                                                             zipper
      s1=s                                                                 zipper
14531 continue                                                             zipper
14532 continue                                                             zipper
14370 continue                                                             zipper
      elc=max(maxval(a-cl(2)),min(minval(a-cl(1)),elc))                    zipper
      return                                                               zipper
      end                                                                  zipper
      function nintot(ni,nx,nc,a,m,nin,is)                                 zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision a(nx,nc)                                            zipper
      integer m(nx),is(ni)                                                 zipper
      is=0                                                                 zipper
      nintot=0                                                             zipper
14560 do 14561 ic=1,nc                                                     zipper
14570 do 14571 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(is(k).ne.0)goto 14571                                             zipper
      if(a(j,ic).eq.0.0)goto 14571                                         zipper
      is(k)=k                                                              zipper
      nintot=nintot+1                                                      zipper
14571 continue                                                             zipper
14572 continue                                                             zipper
14561 continue                                                             zipper
14562 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine luncomp(ni,nx,nc,ca,ia,nin,a)                             zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision ca(nx,nc),a(ni,nc)                                  zipper
      integer ia(nx)                                                       zipper
      a=0.0                                                                zipper
14580 do 14581 ic=1,nc                                                     zipper
      if(nin.gt.0) a(ia(1:nin),ic)=ca(1:nin,ic)                            zipper
14581 continue                                                             zipper
14582 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine lmodval(nt,x,nc,nx,a0,ca,ia,nin,ans)                      zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision a0(nc),ca(nx,nc),x(nt,*),ans(nc,nt)                 zipper
      integer ia(nx)                                                       zipper
14590 do 14591 i=1,nt                                                      zipper
14600 do 14601 ic=1,nc                                                     zipper
      ans(ic,i)=a0(ic)                                                     zipper
      if(nin.gt.0) ans(ic,i)=ans(ic,i)+dot_product(ca(1:nin,ic),x(i,ia(1   zipper
     *:nin)))
14601 continue                                                             zipper
14602 continue                                                             zipper
14591 continue                                                             zipper
14592 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine splognet (parm,no,ni,nc,x,ix,jx,y,g,jd,vp,cl,ne,nx,nlam   zipper
     *,flmin,  ulam,thr,isd,intr,maxit,kopt,lmu,a0,ca,ia,nin,dev0,dev,al
     *m,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no,max(2,nc)),g(no,nc),vp(ni),ulam(nlam)     zipper
      double precision ca(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl   zipper
     *(2,ni)
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           zipper
      double precision, dimension (:), allocatable :: xm,xs,ww,vq,xv            
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 14621                                    zipper
      jerr=10000                                                           zipper
      return                                                               zipper
14621 continue                                                             zipper
      allocate(ww(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(vq(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      if(kopt .ne. 2)goto 14641                                            zipper
      allocate(xv(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
14641 continue                                                             zipper
      call spchkvars(no,ni,x,ix,ju)                                        zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 14661                                      zipper
      jerr=7777                                                            zipper
      return                                                               zipper
14661 continue                                                             zipper
      vq=max(0d0,vp)                                                       zipper
      vq=vq*ni/sum(vq)                                                     zipper
14670 do 14671 i=1,no                                                      zipper
      ww(i)=sum(y(i,:))                                                    zipper
      if(ww(i).gt.0.0) y(i,:)=y(i,:)/ww(i)                                 zipper
14671 continue                                                             zipper
14672 continue                                                             zipper
      sw=sum(ww)                                                           zipper
      ww=ww/sw                                                             zipper
      if(nc .ne. 1)goto 14691                                              zipper
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                zipper
      if(isd .le. 0)goto 14711                                             zipper
14720 do 14721 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
14721 continue                                                             zipper
14722 continue                                                             zipper
14711 continue                                                             zipper
      call sprlognet2n(parm,no,ni,x,ix,jx,y(:,1),g(:,1),ww,ju,vq,cl,ne,n   zipper
     *x,nlam,  flmin,ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,ia,nin
     *,dev0,dev,  alm,nlp,jerr)
      goto 14681                                                           zipper
14691 if(kopt .ne. 2)goto 14731                                            zipper
      call multsplstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs,xv)         zipper
      if(isd .le. 0)goto 14751                                             zipper
14760 do 14761 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
14761 continue                                                             zipper
14762 continue                                                             zipper
14751 continue                                                             zipper
      call multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nl   zipper
     *am,flmin,  ulam,thr,intr,maxit,xv,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,
     *alm,nlp,jerr)
      goto 14771                                                           zipper
14731 continue                                                             zipper
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                zipper
      if(isd .le. 0)goto 14791                                             zipper
14800 do 14801 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
14801 continue                                                             zipper
14802 continue                                                             zipper
14791 continue                                                             zipper
      call sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,f   zipper
     *lmin,  ulam,thr,isd,intr,maxit,kopt,xm,xs,lmu,a0,ca,  ia,nin,dev0,
     *dev,alm,nlp,jerr)
14771 continue                                                             zipper
14681 continue                                                             zipper
      if(jerr.gt.0) return                                                 zipper
      dev0=2.0*sw*dev0                                                     zipper
14810 do 14811 k=1,lmu                                                     zipper
      nk=nin(k)                                                            zipper
14820 do 14821 ic=1,nc                                                     zipper
      if(isd .le. 0)goto 14841                                             zipper
14850 do 14851 l=1,nk                                                      zipper
      ca(l,ic,k)=ca(l,ic,k)/xs(ia(l))                                      zipper
14851 continue                                                             zipper
14852 continue                                                             zipper
14841 continue                                                             zipper
      if(intr .ne. 0)goto 14871                                            zipper
      a0(ic,k)=0.0                                                         zipper
      goto 14881                                                           zipper
14871 continue                                                             zipper
      a0(ic,k)=a0(ic,k)-dot_product(ca(1:nk,ic,k),xm(ia(1:nk)))            zipper
14881 continue                                                             zipper
14861 continue                                                             zipper
14821 continue                                                             zipper
14822 continue                                                             zipper
14811 continue                                                             zipper
14812 continue                                                             zipper
      deallocate(ww,ju,vq,xm,xs)                                           zipper
      if(kopt.eq.2) deallocate(xv)                                         zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multsplstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs,xv)    zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),w(no),xm(ni),xs(ni),xv(ni)                     zipper
      integer ix(*),jx(*),ju(ni)                                           zipper
      if(intr .ne. 0)goto 14901                                            zipper
14910 do 14911 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 14911                                             zipper
      xm(j)=0.0                                                            zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)                          zipper
      if(isd .eq. 0)goto 14931                                             zipper
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            zipper
      vc=xv(j)-xbq                                                         zipper
      xs(j)=sqrt(vc)                                                       zipper
      xv(j)=1.0+xbq/vc                                                     zipper
      goto 14941                                                           zipper
14931 continue                                                             zipper
      xs(j)=1.0                                                            zipper
14941 continue                                                             zipper
14921 continue                                                             zipper
14911 continue                                                             zipper
14912 continue                                                             zipper
      return                                                               zipper
14901 continue                                                             zipper
14950 do 14951 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 14951                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             zipper
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 zipper
      if(isd .le. 0)goto 14971                                             zipper
      xs(j)=sqrt(xv(j))                                                    zipper
      xv(j)=1.0                                                            zipper
14971 continue                                                             zipper
14951 continue                                                             zipper
14952 continue                                                             zipper
      if(isd.eq.0) xs=1.0                                                  zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine splstandard2(no,ni,x,ix,jx,w,ju,isd,intr,xm,xs)           zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),w(no),xm(ni),xs(ni)                            zipper
      integer ix(*),jx(*),ju(ni)                                           zipper
      if(intr .ne. 0)goto 14991                                            zipper
15000 do 15001 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 15001                                             zipper
      xm(j)=0.0                                                            zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      if(isd .eq. 0)goto 15021                                             zipper
      vc=dot_product(w(jx(jb:je)),x(jb:je)**2)  -dot_product(w(jx(jb:je)   zipper
     *),x(jb:je))**zipper
      xs(j)=sqrt(vc)                                                       zipper
      goto 15031                                                           zipper
15021 continue                                                             zipper
      xs(j)=1.0                                                            zipper
15031 continue                                                             zipper
15011 continue                                                             zipper
15001 continue                                                             zipper
15002 continue                                                             zipper
      return                                                               zipper
14991 continue                                                             zipper
15040 do 15041 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 15041                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             zipper
      if(isd.ne.0) xs(j)=sqrt(dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j   zipper
     *)**2)
15041 continue                                                             zipper
15042 continue                                                             zipper
      if(isd.eq.0) xs=1.0                                                  zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine sprlognet2n (parm,no,ni,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,nl   zipper
     *am,  flmin,ulam,shri,isd,intr,maxit,kopt,xb,xs,  lmu,a0,a,m,kin,de
     *v0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)   zipper
      double precision a(nx,nlam),a0(nlam),dev(nlam),alm(nlam)             zipper
      double precision xb(ni),xs(ni)                                       zipper
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           zipper
      double precision, dimension (:), allocatable :: xm,b,bs,v,r               
      double precision, dimension (:), allocatable :: sc,xv,q,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      allocate(b(0:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xm(0:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xv(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(bs(0:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ixx(1:ni),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(q(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(r(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(v(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sc(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      fmax=log(1.0/pmin-1.0)                                               zipper
      fmin=-fmax                                                           zipper
      vmin=(1.0+pmin)*pmin*(1.0-pmin)                                      zipper
      bta=parm                                                             zipper
      omb=1.0-bta                                                          zipper
      q0=dot_product(w,y)                                                  zipper
      if(q0 .gt. pmin)goto 15061                                           zipper
      jerr=8001                                                            zipper
      return                                                               zipper
15061 continue                                                             zipper
      if(q0 .lt. 1.0-pmin)goto 15081                                       zipper
      jerr=9001                                                            zipper
      return                                                               zipper
15081 continue                                                             zipper
      if(intr.eq.0) q0=0.5                                                 zipper
      bz=0.0                                                               zipper
      if(intr.ne.0) bz=log(q0/(1.0-q0))                                    zipper
      if(nonzero(no,g) .ne. 0)goto 15101                                   zipper
      vi=q0*(1.0-q0)                                                       zipper
      b(0)=bz                                                              zipper
      v=vi*w                                                               zipper
      r=w*(y-q0)                                                           zipper
      q=q0                                                                 zipper
      xm(0)=vi                                                             zipper
      dev1=-(bz*q0+log(1.0-q0))                                            zipper
      goto 15111                                                           zipper
15101 continue                                                             zipper
      b(0)=0.0                                                             zipper
      if(intr .eq. 0)goto 15131                                            zipper
      b(0)=azero(no,y,g,w,jerr)                                            zipper
      if(jerr.ne.0) return                                                 zipper
15131 continue                                                             zipper
      q=1.0/(1.0+exp(-b(0)-g))                                             zipper
      v=w*q*(1.0-q)                                                        zipper
      r=w*(y-q)                                                            zipper
      xm(0)=sum(v)                                                         zipper
      dev1=-(b(0)*q0+dot_product(w,y*g+log(1.0-q)))                        zipper
15111 continue                                                             zipper
15091 continue                                                             zipper
      if(kopt .le. 0)goto 15151                                            zipper
      if(isd .le. 0 .or. intr .eq. 0)goto 15171                            zipper
      xv=0.25                                                              zipper
      goto 15181                                                           zipper
15171 continue                                                             zipper
15190 do 15191 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 15191                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xv(j)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)          zipper
15191 continue                                                             zipper
15192 continue                                                             zipper
15181 continue                                                             zipper
15161 continue                                                             zipper
15151 continue                                                             zipper
      b(1:ni)=0.0                                                          zipper
      dev0=dev1                                                            zipper
15200 do 15201 i=1,no                                                      zipper
      if(y(i).gt.0.0) dev0=dev0+w(i)*y(i)*log(y(i))                        zipper
      if(y(i).lt.1.0) dev0=dev0+w(i)*(1.0-y(i))*log(1.0-y(i))              zipper
15201 continue                                                             zipper
15202 continue                                                             zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 15221                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
15221 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nin=0                                                                zipper
      o=0.0                                                                zipper
      svr=o                                                                zipper
      mnl=min(mnlam,nlam)                                                  zipper
      bs=0.0                                                               zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      shr=shri*dev0                                                        zipper
      al=0.0                                                               zipper
      ixx=0                                                                zipper
15230 do 15231 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 15231                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      jn=ix(j+1)-ix(j)                                                     zipper
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 zipper
      gj=dot_product(sc(1:jn),x(jb:je))                                    zipper
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      zipper
15231 continue                                                             zipper
15232 continue                                                             zipper
15240 do 15241 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 15261                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 15251                                                           zipper
15261 if(ilm .le. 2)goto 15271                                             zipper
      al=al*alf                                                            zipper
      goto 15251                                                           zipper
15271 if(ilm .ne. 1)goto 15281                                             zipper
      al=big                                                               zipper
      goto 15291                                                           zipper
15281 continue                                                             zipper
      al0=0.0                                                              zipper
15300 do 15301 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 15301                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
15301 continue                                                             zipper
15302 continue                                                             zipper
      al0=al0/max(bta,1.0d-3)                                              zipper
      al=alf*al0                                                           zipper
15291 continue                                                             zipper
15251 continue                                                             zipper
      al2=al*omb                                                           zipper
      al1=al*bta                                                           zipper
      tlam=bta*(2.0*al-al0)                                                zipper
15310 do 15311 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 15311                                            zipper
      if(ju(k).eq.0)goto 15311                                             zipper
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     zipper
15311 continue                                                             zipper
15312 continue                                                             zipper
10880 continue                                                             zipper
15320 continue                                                             zipper
15321 continue                                                             zipper
      bs(0)=b(0)                                                           zipper
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))                                zipper
15330 do 15331 j=1,ni                                                      zipper
      if(ixx(j).eq.0)goto 15331                                            zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      jn=ix(j+1)-ix(j)                                                     zipper
      sc(1:jn)=v(jx(jb:je))                                                zipper
      xm(j)=dot_product(sc(1:jn),x(jb:je))                                 zipper
      if(kopt .ne. 0)goto 15351                                            zipper
      xv(j)=dot_product(sc(1:jn),x(jb:je)**2)                              zipper
      xv(j)=(xv(j)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2                zipper
15351 continue                                                             zipper
15331 continue                                                             zipper
15332 continue                                                             zipper
15360 continue                                                             zipper
15361 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
15370 do 15371 k=1,ni                                                      zipper
      if(ixx(k).eq.0)goto 15371                                            zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      jn=ix(k+1)-ix(k)                                                     zipper
      bk=b(k)                                                              zipper
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 zipper
      gk=dot_product(sc(1:jn),x(jb:je))                                    zipper
      gk=(gk-svr*xb(k))/xs(k)                                              zipper
      u=gk+xv(k)*b(k)                                                      zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 15391                                            zipper
      b(k)=0.0                                                             zipper
      goto 15401                                                           zipper
15391 continue                                                             zipper
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          zipper
15401 continue                                                             zipper
15381 continue                                                             zipper
      d=b(k)-bk                                                            zipper
      if(abs(d).le.0.0)goto 15371                                          zipper
      dlx=max(dlx,xv(k)*d**2)                                              zipper
      if(mm(k) .ne. 0)goto 15421                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 15372                                              zipper
      mm(k)=nin                                                            zipper
      m(nin)=k                                                             zipper
      sc(1:jn)=v(jx(jb:je))                                                zipper
      xm(k)=dot_product(sc(1:jn),x(jb:je))                                 zipper
15421 continue                                                             zipper
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              zipper
      o=o+d*(xb(k)/xs(k))                                                  zipper
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  zipper
15371 continue                                                             zipper
15372 continue                                                             zipper
      if(nin.gt.nx)goto 15362                                              zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=svr/xm(0)                                            zipper
      if(d .eq. 0.0)goto 15441                                             zipper
      b(0)=b(0)+d                                                          zipper
      dlx=max(dlx,xm(0)*d**2)                                              zipper
      r=r-d*v                                                              zipper
      svr=svr-d*xm(0)                                                      zipper
15441 continue                                                             zipper
      if(dlx.lt.shr)goto 15362                                             zipper
      if(nlp .le. maxit)goto 15461                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
15461 continue                                                             zipper
15470 continue                                                             zipper
15471 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
15480 do 15481 l=1,nin                                                     zipper
      k=m(l)                                                               zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      jn=ix(k+1)-ix(k)                                                     zipper
      bk=b(k)                                                              zipper
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 zipper
      gk=dot_product(sc(1:jn),x(jb:je))                                    zipper
      gk=(gk-svr*xb(k))/xs(k)                                              zipper
      u=gk+xv(k)*b(k)                                                      zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 15501                                            zipper
      b(k)=0.0                                                             zipper
      goto 15511                                                           zipper
15501 continue                                                             zipper
      b(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(xv(k)+vp(k)*al2)))          zipper
15511 continue                                                             zipper
15491 continue                                                             zipper
      d=b(k)-bk                                                            zipper
      if(abs(d).le.0.0)goto 15481                                          zipper
      dlx=max(dlx,xv(k)*d**2)                                              zipper
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              zipper
      o=o+d*(xb(k)/xs(k))                                                  zipper
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  zipper
15481 continue                                                             zipper
15482 continue                                                             zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=svr/xm(0)                                            zipper
      if(d .eq. 0.0)goto 15531                                             zipper
      b(0)=b(0)+d                                                          zipper
      dlx=max(dlx,xm(0)*d**2)                                              zipper
      r=r-d*v                                                              zipper
      svr=svr-d*xm(0)                                                      zipper
15531 continue                                                             zipper
      if(dlx.lt.shr)goto 15472                                             zipper
      if(nlp .le. maxit)goto 15551                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
15551 continue                                                             zipper
      goto 15471                                                           zipper
15472 continue                                                             zipper
      goto 15361                                                           zipper
15362 continue                                                             zipper
      if(nin.gt.nx)goto 15322                                              zipper
      sc=b(0)                                                              zipper
      b0=0.0                                                               zipper
15560 do 15561 j=1,nin                                                     zipper
      l=m(j)                                                               zipper
      jb=ix(l)                                                             zipper
      je=ix(l+1)-1                                                         zipper
      sc(jx(jb:je))=sc(jx(jb:je))+b(l)*x(jb:je)/xs(l)                      zipper
      b0=b0-b(l)*xb(l)/xs(l)                                               zipper
15561 continue                                                             zipper
15562 continue                                                             zipper
      sc=sc+b0                                                             zipper
15570 do 15571 i=1,no                                                      zipper
      fi=sc(i)+g(i)                                                        zipper
      if(fi .ge. fmin)goto 15591                                           zipper
      q(i)=0.0                                                             zipper
      goto 15581                                                           zipper
15591 if(fi .le. fmax)goto 15601                                           zipper
      q(i)=1.0                                                             zipper
      goto 15611                                                           zipper
15601 continue                                                             zipper
      q(i)=1.0/(1.0+exp(-fi))                                              zipper
15611 continue                                                             zipper
15581 continue                                                             zipper
15571 continue                                                             zipper
15572 continue                                                             zipper
      v=w*q*(1.0-q)                                                        zipper
      xm(0)=sum(v)                                                         zipper
      if(xm(0).lt.vmin)goto 15322                                          zipper
      r=w*(y-q)                                                            zipper
      svr=sum(r)                                                           zipper
      o=0.0                                                                zipper
      if(xm(0)*(b(0)-bs(0))**2 .ge. shr)goto 15631                         zipper
      kx=0                                                                 zipper
15640 do 15641 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 15641                           zipper
      kx=1                                                                 zipper
      goto 15642                                                           zipper
15641 continue                                                             zipper
15642 continue                                                             zipper
      if(kx .ne. 0)goto 15661                                              zipper
15670 do 15671 j=1,ni                                                      zipper
      if(ixx(j).eq.1)goto 15671                                            zipper
      if(ju(j).eq.0)goto 15671                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      jn=ix(j+1)-ix(j)                                                     zipper
      sc(1:jn)=r(jx(jb:je))+v(jx(jb:je))*o                                 zipper
      gj=dot_product(sc(1:jn),x(jb:je))                                    zipper
      ga(j)=abs((gj-svr*xb(j))/xs(j))                                      zipper
      if(ga(j) .le. al1*vp(j))goto 15691                                   zipper
      ixx(j)=1                                                             zipper
      kx=1                                                                 zipper
15691 continue                                                             zipper
15671 continue                                                             zipper
15672 continue                                                             zipper
      if(kx.eq.1) go to 10880                                              zipper
      goto 15322                                                           zipper
15661 continue                                                             zipper
15631 continue                                                             zipper
      goto 15321                                                           zipper
15322 continue                                                             zipper
      if(nin .le. nx)goto 15711                                            zipper
      jerr=-10000-ilm                                                      zipper
      goto 15242                                                           zipper
15711 continue                                                             zipper
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))                                zipper
      kin(ilm)=nin                                                         zipper
      a0(ilm)=b(0)                                                         zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      devi=dev2(no,w,y,q,pmin)                                             zipper
      dev(ilm)=(dev1-devi)/dev0                                            zipper
      if(ilm.lt.mnl)goto 15241                                             zipper
      if(flmin.ge.1.0)goto 15241                                           zipper
      me=0                                                                 zipper
15720 do 15721 j=1,nin                                                     zipper
      if(a(j,ilm).ne.0.0) me=me+1                                          zipper
15721 continue                                                             zipper
15722 continue                                                             zipper
      if(me.gt.ne)goto 15242                                               zipper
      if(dev(ilm).gt.devmax)goto 15242                                     zipper
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 15242                             zipper
      if(xm(0).lt.vmin)goto 15242                                          zipper
15241 continue                                                             zipper
15242 continue                                                             zipper
      g=log(q/(1.0-q))                                                     zipper
      deallocate(xm,b,bs,v,r,sc,xv,q,mm,ga,ixx)                            zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine sprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,nx,n   zipper
     *lam,flmin,  ulam,shri,isd,intr,maxit,kopt,xb,xs,lmu,a0,a,m,kin,dev
     *0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam),xb   zipper
     *(ni),xs(ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   zipper
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           zipper
      double precision, dimension (:,:), allocatable :: q                       
      double precision, dimension (:), allocatable :: sxp,sxpl                  
      double precision, dimension (:), allocatable :: sc,xm,v,r,ga              
      double precision, dimension (:,:), allocatable :: b,bs,xv                 
      integer, dimension (:), allocatable :: mm,is,iy                           
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(xv(1:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return 					                                                
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return					                                                 
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      exmn=-exmx                                                           zipper
      allocate(xm(0:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(r(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(v(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(iy(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(is(1:max(nc,ni)),stat=jerr)                                 zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sxp(1:no),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sxpl(1:no),stat=jerr)                                       zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sc(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      pmax=1.0-pmin                                                        zipper
      emin=pmin/pmax                                                       zipper
      emax=1.0/emin                                                        zipper
      pfm=(1.0+pmin)*pmin                                                  zipper
      pfx=(1.0-pmin)*pmax                                                  zipper
      vmin=pfm*pmax                                                        zipper
      bta=parm                                                             zipper
      omb=1.0-bta                                                          zipper
      dev1=0.0                                                             zipper
      dev0=0.0                                                             zipper
15730 do 15731 ic=1,nc                                                     zipper
      q0=dot_product(w,y(:,ic))                                            zipper
      if(q0 .gt. pmin)goto 15751                                           zipper
      jerr =8000+ic                                                        zipper
      return                                                               zipper
15751 continue                                                             zipper
      if(q0 .lt. 1.0-pmin)goto 15771                                       zipper
      jerr =9000+ic                                                        zipper
      return                                                               zipper
15771 continue                                                             zipper
      if(intr.eq.0) q0=1.0/nc                                              zipper
      b(1:ni,ic)=0.0                                                       zipper
      b(0,ic)=0.0                                                          zipper
      if(intr .eq. 0)goto 15791                                            zipper
      b(0,ic)=log(q0)                                                      zipper
      dev1=dev1-q0*b(0,ic)                                                 zipper
15791 continue                                                             zipper
15731 continue                                                             zipper
15732 continue                                                             zipper
      if(intr.eq.0) dev1=log(float(nc))                                    zipper
      iy=0                                                                 zipper
      al=0.0                                                               zipper
      if(nonzero(no*nc,g) .ne. 0)goto 15811                                zipper
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         zipper
      sxp=0.0                                                              zipper
15820 do 15821 ic=1,nc                                                     zipper
      q(:,ic)=exp(b(0,ic))                                                 zipper
      sxp=sxp+q(:,ic)                                                      zipper
15821 continue                                                             zipper
15822 continue                                                             zipper
      goto 15831                                                           zipper
15811 continue                                                             zipper
15840 do 15841 i=1,no                                                      zipper
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         zipper
15841 continue                                                             zipper
15842 continue                                                             zipper
      sxp=0.0                                                              zipper
      if(intr .ne. 0)goto 15861                                            zipper
      b(0,:)=0.0                                                           zipper
      goto 15871                                                           zipper
15861 continue                                                             zipper
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 zipper
      if(jerr.ne.0) return                                                 zipper
15871 continue                                                             zipper
15851 continue                                                             zipper
      dev1=0.0                                                             zipper
15880 do 15881 ic=1,nc                                                     zipper
      q(:,ic)=b(0,ic)+g(:,ic)                                              zipper
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             zipper
      q(:,ic)=exp(q(:,ic))                                                 zipper
      sxp=sxp+q(:,ic)                                                      zipper
15881 continue                                                             zipper
15882 continue                                                             zipper
      sxpl=w*log(sxp)                                                      zipper
15890 do 15891 ic=1,nc                                                     zipper
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  zipper
15891 continue                                                             zipper
15892 continue                                                             zipper
15831 continue                                                             zipper
15801 continue                                                             zipper
15900 do 15901 ic=1,nc                                                     zipper
15910 do 15911 i=1,no                                                      zipper
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               zipper
15911 continue                                                             zipper
15912 continue                                                             zipper
15901 continue                                                             zipper
15902 continue                                                             zipper
      dev0=dev0+dev1                                                       zipper
      if(kopt .le. 0)goto 15931                                            zipper
      if(isd .le. 0 .or. intr .eq. 0)goto 15951                            zipper
      xv=0.25                                                              zipper
      goto 15961                                                           zipper
15951 continue                                                             zipper
15970 do 15971 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 15971                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xv(j,:)=0.25*(dot_product(w(jx(jb:je)),x(jb:je)**2)-xb(j)**2)        zipper
15971 continue                                                             zipper
15972 continue                                                             zipper
15961 continue                                                             zipper
15941 continue                                                             zipper
15931 continue                                                             zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 15991                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
15991 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nin=0                                                                zipper
      nlp=0                                                                zipper
      mnl=min(mnlam,nlam)                                                  zipper
      bs=0.0                                                               zipper
      svr=0.0                                                              zipper
      o=0.0                                                                zipper
      shr=shri*dev0                                                        zipper
      ga=0.0                                                               zipper
16000 do 16001 ic=1,nc                                                     zipper
      v=q(:,ic)/sxp                                                        zipper
      r=w*(y(:,ic)-v)                                                      zipper
      v=w*v*(1.0-v)                                                        zipper
16010 do 16011 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 16011                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      jn=ix(j+1)-ix(j)                                                     zipper
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 zipper
      gj=dot_product(sc(1:jn),x(jb:je))                                    zipper
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             zipper
16011 continue                                                             zipper
16012 continue                                                             zipper
16001 continue                                                             zipper
16002 continue                                                             zipper
16020 do 16021 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 16041                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 16031                                                           zipper
16041 if(ilm .le. 2)goto 16051                                             zipper
      al=al*alf                                                            zipper
      goto 16031                                                           zipper
16051 if(ilm .ne. 1)goto 16061                                             zipper
      al=big                                                               zipper
      goto 16071                                                           zipper
16061 continue                                                             zipper
      al0=0.0                                                              zipper
16080 do 16081 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 16081                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
16081 continue                                                             zipper
16082 continue                                                             zipper
      al0=al0/max(bta,1.0d-3)                                              zipper
      al=alf*al0                                                           zipper
16071 continue                                                             zipper
16031 continue                                                             zipper
      al2=al*omb                                                           zipper
      al1=al*bta                                                           zipper
      tlam=bta*(2.0*al-al0)                                                zipper
16090 do 16091 k=1,ni                                                      zipper
      if(iy(k).eq.1)goto 16091                                             zipper
      if(ju(k).eq.0)goto 16091                                             zipper
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      zipper
16091 continue                                                             zipper
16092 continue                                                             zipper
10880 continue                                                             zipper
16100 continue                                                             zipper
16101 continue                                                             zipper
      ixx=0                                                                zipper
      jxx=ixx                                                              zipper
      ig=0                                                                 zipper
16110 do 16111 ic=1,nc                                                     zipper
      bs(0,ic)=b(0,ic)                                                     zipper
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          zipper
      xm(0)=0.0                                                            zipper
      svr=0.0                                                              zipper
      o=0.0                                                                zipper
16120 do 16121 i=1,no                                                      zipper
      pic=q(i,ic)/sxp(i)                                                   zipper
      if(pic .ge. pfm)goto 16141                                           zipper
      pic=0.0                                                              zipper
      v(i)=0.0                                                             zipper
      goto 16131                                                           zipper
16141 if(pic .le. pfx)goto 16151                                           zipper
      pic=1.0                                                              zipper
      v(i)=0.0                                                             zipper
      goto 16161                                                           zipper
16151 continue                                                             zipper
      v(i)=w(i)*pic*(1.0-pic)                                              zipper
      xm(0)=xm(0)+v(i)                                                     zipper
16161 continue                                                             zipper
16131 continue                                                             zipper
      r(i)=w(i)*(y(i,ic)-pic)                                              zipper
      svr=svr+r(i)                                                         zipper
16121 continue                                                             zipper
16122 continue                                                             zipper
      if(xm(0).le.vmin)goto 16111                                          zipper
      ig=1                                                                 zipper
16170 do 16171 j=1,ni                                                      zipper
      if(iy(j).eq.0)goto 16171                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(v(jx(jb:je)),x(jb:je))                             zipper
      if(kopt .ne. 0)goto 16191                                            zipper
      xv(j,ic)=dot_product(v(jx(jb:je)),x(jb:je)**2)                       zipper
      xv(j,ic)=(xv(j,ic)-2.0*xb(j)*xm(j)+xm(0)*xb(j)**2)/xs(j)**2          zipper
16191 continue                                                             zipper
16171 continue                                                             zipper
16172 continue                                                             zipper
16200 continue                                                             zipper
16201 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
16210 do 16211 k=1,ni                                                      zipper
      if(iy(k).eq.0)goto 16211                                             zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      jn=ix(k+1)-ix(k)                                                     zipper
      bk=b(k,ic)                                                           zipper
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 zipper
      gk=dot_product(sc(1:jn),x(jb:je))                                    zipper
      gk=(gk-svr*xb(k))/xs(k)                                              zipper
      u=gk+xv(k,ic)*b(k,ic)                                                zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 16231                                            zipper
      b(k,ic)=0.0                                                          zipper
      goto 16241                                                           zipper
16231 continue                                                             zipper
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   zipper
     *)
16241 continue                                                             zipper
16221 continue                                                             zipper
      d=b(k,ic)-bk                                                         zipper
      if(abs(d).le.0.0)goto 16211                                          zipper
      dlx=max(dlx,xv(k,ic)*d**2)                                           zipper
      if(mm(k) .ne. 0)goto 16261                                           zipper
      nin=nin+1                                                            zipper
      if(nin .le. nx)goto 16281                                            zipper
      jxx=1                                                                zipper
      goto 16212                                                           zipper
16281 continue                                                             zipper
      mm(k)=nin                                                            zipper
      m(nin)=k                                                             zipper
      xm(k)=dot_product(v(jx(jb:je)),x(jb:je))                             zipper
16261 continue                                                             zipper
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              zipper
      o=o+d*(xb(k)/xs(k))                                                  zipper
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  zipper
16211 continue                                                             zipper
16212 continue                                                             zipper
      if(jxx.gt.0)goto 16202                                               zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=svr/xm(0)                                            zipper
      if(d .eq. 0.0)goto 16301                                             zipper
      b(0,ic)=b(0,ic)+d                                                    zipper
      dlx=max(dlx,xm(0)*d**2)                                              zipper
      r=r-d*v                                                              zipper
      svr=svr-d*xm(0)                                                      zipper
16301 continue                                                             zipper
      if(dlx.lt.shr)goto 16202                                             zipper
      if(nlp .le. maxit)goto 16321                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
16321 continue                                                             zipper
16330 continue                                                             zipper
16331 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
16340 do 16341 l=1,nin                                                     zipper
      k=m(l)                                                               zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      jn=ix(k+1)-ix(k)                                                     zipper
      bk=b(k,ic)                                                           zipper
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 zipper
      gk=dot_product(sc(1:jn),x(jb:je))                                    zipper
      gk=(gk-svr*xb(k))/xs(k)                                              zipper
      u=gk+xv(k,ic)*b(k,ic)                                                zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 16361                                            zipper
      b(k,ic)=0.0                                                          zipper
      goto 16371                                                           zipper
16361 continue                                                             zipper
      b(k,ic)=max(cl(1,k),min(cl(2,k),sign(au,u)/  (xv(k,ic)+vp(k)*al2))   zipper
     *)
16371 continue                                                             zipper
16351 continue                                                             zipper
      d=b(k,ic)-bk                                                         zipper
      if(abs(d).le.0.0)goto 16341                                          zipper
      dlx=max(dlx,xv(k,ic)*d**2)                                           zipper
      r(jx(jb:je))=r(jx(jb:je))-d*v(jx(jb:je))*x(jb:je)/xs(k)              zipper
      o=o+d*(xb(k)/xs(k))                                                  zipper
      svr=svr-d*(xm(k)-xb(k)*xm(0))/xs(k)                                  zipper
16341 continue                                                             zipper
16342 continue                                                             zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=svr/xm(0)                                            zipper
      if(d .eq. 0.0)goto 16391                                             zipper
      b(0,ic)=b(0,ic)+d                                                    zipper
      dlx=max(dlx,xm(0)*d**2)                                              zipper
      r=r-d*v                                                              zipper
      svr=svr-d*xm(0)                                                      zipper
16391 continue                                                             zipper
      if(dlx.lt.shr)goto 16332                                             zipper
      if(nlp .le. maxit)goto 16411                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
16411 continue                                                             zipper
      goto 16331                                                           zipper
16332 continue                                                             zipper
      goto 16201                                                           zipper
16202 continue                                                             zipper
      if(jxx.gt.0)goto 16112                                               zipper
      if(xm(0)*(b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                         zipper
      if(ixx .ne. 0)goto 16431                                             zipper
16440 do 16441 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(xv(k,ic)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 16461                zipper
      ixx=1                                                                zipper
      goto 16442                                                           zipper
16461 continue                                                             zipper
16441 continue                                                             zipper
16442 continue                                                             zipper
16431 continue                                                             zipper
      sc=b(0,ic)+g(:,ic)                                                   zipper
      b0=0.0                                                               zipper
16470 do 16471 j=1,nin                                                     zipper
      l=m(j)                                                               zipper
      jb=ix(l)                                                             zipper
      je=ix(l+1)-1                                                         zipper
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   zipper
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            zipper
16471 continue                                                             zipper
16472 continue                                                             zipper
      sc=min(max(exmn,sc+b0),exmx)                                         zipper
      sxp=sxp-q(:,ic)                                                      zipper
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          zipper
      sxp=sxp+q(:,ic)                                                      zipper
16111 continue                                                             zipper
16112 continue                                                             zipper
      s=-sum(b(0,:))/nc                                                    zipper
      b(0,:)=b(0,:)+s                                                      zipper
      sc=s                                                                 zipper
      b0=0.0                                                               zipper
16480 do 16481 j=1,nin                                                     zipper
      l=m(j)                                                               zipper
      if(vp(l) .gt. 0.0)goto 16501                                         zipper
      s=sum(b(l,:))/nc                                                     zipper
      goto 16511                                                           zipper
16501 continue                                                             zipper
      s=elc(parm,nc,cl(:,l),b(l,:),is)                                     zipper
16511 continue                                                             zipper
16491 continue                                                             zipper
      b(l,:)=b(l,:)-s                                                      zipper
      jb=ix(l)                                                             zipper
      je=ix(l+1)-1                                                         zipper
      sc(jx(jb:je))=sc(jx(jb:je))-s*x(jb:je)/xs(l)                         zipper
      b0=b0+s*xb(l)/xs(l)                                                  zipper
16481 continue                                                             zipper
16482 continue                                                             zipper
      sc=sc+b0                                                             zipper
      sc=exp(sc)                                                           zipper
      sxp=sxp*sc                                                           zipper
16520 do 16521 ic=1,nc                                                     zipper
      q(:,ic)=q(:,ic)*sc                                                   zipper
16521 continue                                                             zipper
16522 continue                                                             zipper
      if(jxx.gt.0)goto 16102                                               zipper
      if(ig.eq.0)goto 16102                                                zipper
      if(ixx .ne. 0)goto 16541                                             zipper
16550 do 16551 j=1,ni                                                      zipper
      if(iy(j).eq.1)goto 16551                                             zipper
      if(ju(j).eq.0)goto 16551                                             zipper
      ga(j)=0.0                                                            zipper
16551 continue                                                             zipper
16552 continue                                                             zipper
16560 do 16561 ic=1,nc                                                     zipper
      v=q(:,ic)/sxp                                                        zipper
      r=w*(y(:,ic)-v)                                                      zipper
      v=w*v*(1.0-v)                                                        zipper
16570 do 16571 j=1,ni                                                      zipper
      if(iy(j).eq.1)goto 16571                                             zipper
      if(ju(j).eq.0)goto 16571                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      jn=ix(j+1)-ix(j)                                                     zipper
      sc(1:jn)=r(jx(jb:je))+o*v(jx(jb:je))                                 zipper
      gj=dot_product(sc(1:jn),x(jb:je))                                    zipper
      ga(j)=max(ga(j),abs(gj-svr*xb(j))/xs(j))                             zipper
16571 continue                                                             zipper
16572 continue                                                             zipper
16561 continue                                                             zipper
16562 continue                                                             zipper
16580 do 16581 k=1,ni                                                      zipper
      if(iy(k).eq.1)goto 16581                                             zipper
      if(ju(k).eq.0)goto 16581                                             zipper
      if(ga(k) .le. al1*vp(k))goto 16601                                   zipper
      iy(k)=1                                                              zipper
      ixx=1                                                                zipper
16601 continue                                                             zipper
16581 continue                                                             zipper
16582 continue                                                             zipper
      if(ixx.eq.1) go to 10880                                             zipper
      goto 16102                                                           zipper
16541 continue                                                             zipper
      goto 16101                                                           zipper
16102 continue                                                             zipper
      if(jxx .le. 0)goto 16621                                             zipper
      jerr=-10000-ilm                                                      zipper
      goto 16022                                                           zipper
16621 continue                                                             zipper
      devi=0.0                                                             zipper
16630 do 16631 ic=1,nc                                                     zipper
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          zipper
      a0(ic,ilm)=b(0,ic)                                                   zipper
16640 do 16641 i=1,no                                                      zipper
      if(y(i,ic).le.0.0)goto 16641                                         zipper
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           zipper
16641 continue                                                             zipper
16642 continue                                                             zipper
16631 continue                                                             zipper
16632 continue                                                             zipper
      kin(ilm)=nin                                                         zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      dev(ilm)=(dev1-devi)/dev0                                            zipper
      if(ig.eq.0)goto 16022                                                zipper
      if(ilm.lt.mnl)goto 16021                                             zipper
      if(flmin.ge.1.0)goto 16021                                           zipper
      if(nintot(ni,nx,nc,a(1,1,ilm),m,nin,is).gt.ne)goto 16022             zipper
      if(dev(ilm).gt.devmax)goto 16022                                     zipper
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 16022                             zipper
16021 continue                                                             zipper
16022 continue                                                             zipper
      g=log(q)                                                             zipper
16650 do 16651 i=1,no                                                      zipper
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         zipper
16651 continue                                                             zipper
16652 continue                                                             zipper
      deallocate(sxp,b,bs,v,r,xv,q,mm,is,xm,sc,ga,iy)                      zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine lcmodval(nc,nx,a0,ca,ia,nin,x,ix,jx,n,f)                  zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision a0(nc),ca(nx,nc),x(*),f(nc,n)                       zipper
      integer ia(*),ix(*),jx(*)                                            zipper
16660 do 16661 ic=1,nc                                                     zipper
      f(ic,:)=a0(ic)                                                       zipper
16661 continue                                                             zipper
16662 continue                                                             zipper
16670 do 16671 j=1,nin                                                     zipper
      k=ia(j)                                                              zipper
      kb=ix(k)                                                             zipper
      ke=ix(k+1)-1                                                         zipper
16680 do 16681 ic=1,nc                                                     zipper
      f(ic,jx(kb:ke))=f(ic,jx(kb:ke))+ca(j,ic)*x(kb:ke)                    zipper
16681 continue                                                             zipper
16682 continue                                                             zipper
16671 continue                                                             zipper
16672 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine coxnet (parm,no,ni,x,y,d,g,w,jd,vp,cl,ne,nx,nlam,flmin,   zipper
     *ulam,thr,  maxit,isd,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no),d(no),g(no),w(no),vp(ni),ulam(nlam   zipper
     *)
      double precision ca(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)            zipper
      integer jd(*),ia(nx),nin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: xs,ww,vq                  
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 16701                                    zipper
      jerr=10000                                                           zipper
      return                                                               zipper
16701 continue                                                             zipper
      allocate(ww(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(vq(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      if(isd .le. 0)goto 16721                                             zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
16721 continue                                                             zipper
      call chkvars(no,ni,x,ju)                                             zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 16741                                      zipper
      jerr=7777                                                            zipper
      return                                                               zipper
16741 continue                                                             zipper
      vq=max(0d0,vp)                                                       zipper
      vq=vq*ni/sum(vq)                                                     zipper
      ww=max(0d0,w)                                                        zipper
      sw=sum(ww)                                                           zipper
      if(sw .gt. 0.0)goto 16761                                            zipper
      jerr=9999                                                            zipper
      return                                                               zipper
16761 continue                                                             zipper
      ww=ww/sw                                                             zipper
      call cstandard(no,ni,x,ww,ju,isd,xs)                                 zipper
      if(isd .le. 0)goto 16781                                             zipper
16790 do 16791 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
16791 continue                                                             zipper
16792 continue                                                             zipper
16781 continue                                                             zipper
      call coxnet1(parm,no,ni,x,y,d,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,   zipper
     *thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 zipper
      dev0=2.0*sw*dev0                                                     zipper
      if(isd .le. 0)goto 16811                                             zipper
16820 do 16821 k=1,lmu                                                     zipper
      nk=nin(k)                                                            zipper
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                                   zipper
16821 continue                                                             zipper
16822 continue                                                             zipper
16811 continue                                                             zipper
      deallocate(ww,ju,vq)                                                 zipper
      if(isd.gt.0) deallocate(xs)                                          zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine cstandard (no,ni,x,w,ju,isd,xs)                           zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),w(no),xs(ni)                               zipper
      integer ju(ni)                                                       zipper
16830 do 16831 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 16831                                             zipper
      xm=dot_product(w,x(:,j))                                             zipper
      x(:,j)=x(:,j)-xm                                                     zipper
      if(isd .le. 0)goto 16851                                             zipper
      xs(j)=sqrt(dot_product(w,x(:,j)**2))                                 zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
16851 continue                                                             zipper
16831 continue                                                             zipper
16832 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine coxnet1(parm,no,ni,x,y,d,g,q,ju,vp,cl,ne,nx,nlam,flmin,   zipper
     *ulam,cthri,  isd,maxit,lmu,ao,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no),q(no),d(no),g(no),vp(ni),ulam(nlam   zipper
     *)
      double precision ao(nx,nlam),dev(nlam),alm(nlam),cl(2,ni)            zipper
      integer ju(ni),m(nx),kin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: w,dk,v,xs,wr              
      double precision, dimension (:), allocatable :: a,as,f,dq                 
      double precision, dimension (:), allocatable :: e,uu,ga                   
      integer, dimension (:), allocatable :: jp,kp,mm,ixx                       
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      sml=sml*100.0                                                        zipper
      devmax=devmax*0.99/0.999                                             zipper
      allocate(e(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(uu(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(f(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(w(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(v(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(a(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(as(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(ixx(1:ni),stat=jerr)                                        zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(jp(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(kp(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(dk(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(wr(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(dq(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0)go to 12180                                             zipper
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               zipper
      if(jerr.ne.0) go to 12180                                            zipper
      alpha=parm                                                           zipper
      oma=1.0-alpha                                                        zipper
      nlm=0                                                                zipper
      ixx=0                                                                zipper
      al=0.0                                                               zipper
      dq=d*q                                                               zipper
      call died(no,nk,dq,kp,jp,dk)                                         zipper
      a=0.0                                                                zipper
      f(1)=0.0                                                             zipper
      fmax=log(huge(f(1))*0.1)                                             zipper
      if(nonzero(no,g) .eq. 0)goto 16871                                   zipper
      f=g-dot_product(q,g)                                                 zipper
      e=q*exp(sign(min(abs(f),fmax),f))                                    zipper
      goto 16881                                                           zipper
16871 continue                                                             zipper
      f=0.0                                                                zipper
      e=q                                                                  zipper
16881 continue                                                             zipper
16861 continue                                                             zipper
      r0=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                                 zipper
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)                         zipper
      dev0=rr                                                              zipper
16890 do 16891 i=1,no                                                      zipper
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 16911                   zipper
      w(i)=0.0                                                             zipper
      wr(i)=w(i)                                                           zipper
16911 continue                                                             zipper
16891 continue                                                             zipper
16892 continue                                                             zipper
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         zipper
      if(jerr.ne.0) go to 12180                                            zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 16931                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
16931 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      mnl=min(mnlam,nlam)                                                  zipper
      as=0.0                                                               zipper
      cthr=cthri*dev0                                                      zipper
16940 do 16941 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 16941                                             zipper
      ga(j)=abs(dot_product(wr,x(:,j)))                                    zipper
16941 continue                                                             zipper
16942 continue                                                             zipper
16950 do 16951 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 16971                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 16961                                                           zipper
16971 if(ilm .le. 2)goto 16981                                             zipper
      al=al*alf                                                            zipper
      goto 16961                                                           zipper
16981 if(ilm .ne. 1)goto 16991                                             zipper
      al=big                                                               zipper
      goto 17001                                                           zipper
16991 continue                                                             zipper
      al0=0.0                                                              zipper
17010 do 17011 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 17011                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
17011 continue                                                             zipper
17012 continue                                                             zipper
      al0=al0/max(parm,1.0d-3)                                             zipper
      al=alf*al0                                                           zipper
17001 continue                                                             zipper
16961 continue                                                             zipper
      sa=alpha*al                                                          zipper
      omal=oma*al                                                          zipper
      tlam=alpha*(2.0*al-al0)                                              zipper
17020 do 17021 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 17021                                            zipper
      if(ju(k).eq.0)goto 17021                                             zipper
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     zipper
17021 continue                                                             zipper
17022 continue                                                             zipper
10880 continue                                                             zipper
17030 continue                                                             zipper
17031 continue                                                             zipper
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                zipper
      call vars(no,ni,x,w,ixx,v)                                           zipper
17040 continue                                                             zipper
17041 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dli=0.0                                                              zipper
17050 do 17051 j=1,ni                                                      zipper
      if(ixx(j).eq.0)goto 17051                                            zipper
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   zipper
      if(abs(u) .gt. vp(j)*sa)goto 17071                                   zipper
      at=0.0                                                               zipper
      goto 17081                                                           zipper
17071 continue                                                             zipper
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   zipper
     *mal)))
17081 continue                                                             zipper
17061 continue                                                             zipper
      if(at .eq. a(j))goto 17101                                           zipper
      del=at-a(j)                                                          zipper
      a(j)=at                                                              zipper
      dli=max(dli,v(j)*del**2)                                             zipper
      wr=wr-del*w*x(:,j)                                                   zipper
      f=f+del*x(:,j)                                                       zipper
      if(mm(j) .ne. 0)goto 17121                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 17052                                              zipper
      mm(j)=nin                                                            zipper
      m(nin)=j                                                             zipper
17121 continue                                                             zipper
17101 continue                                                             zipper
17051 continue                                                             zipper
17052 continue                                                             zipper
      if(nin.gt.nx)goto 17042                                              zipper
      if(dli.lt.cthr)goto 17042                                            zipper
      if(nlp .le. maxit)goto 17141                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
17141 continue                                                             zipper
17150 continue                                                             zipper
17151 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dli=0.0                                                              zipper
17160 do 17161 l=1,nin                                                     zipper
      j=m(l)                                                               zipper
      u=a(j)*v(j)+dot_product(wr,x(:,j))                                   zipper
      if(abs(u) .gt. vp(j)*sa)goto 17181                                   zipper
      at=0.0                                                               zipper
      goto 17191                                                           zipper
17181 continue                                                             zipper
      at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*sa,u)/  (v(j)+vp(j)*o   zipper
     *mal)))
17191 continue                                                             zipper
17171 continue                                                             zipper
      if(at .eq. a(j))goto 17211                                           zipper
      del=at-a(j)                                                          zipper
      a(j)=at                                                              zipper
      dli=max(dli,v(j)*del**2)                                             zipper
      wr=wr-del*w*x(:,j)                                                   zipper
      f=f+del*x(:,j)                                                       zipper
17211 continue                                                             zipper
17161 continue                                                             zipper
17162 continue                                                             zipper
      if(dli.lt.cthr)goto 17152                                            zipper
      if(nlp .le. maxit)goto 17231                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
17231 continue                                                             zipper
      goto 17151                                                           zipper
17152 continue                                                             zipper
      goto 17041                                                           zipper
17042 continue                                                             zipper
      if(nin.gt.nx)goto 17032                                              zipper
      e=q*exp(sign(min(abs(f),fmax),f))                                    zipper
      call outer(no,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)                         zipper
      if(jerr .eq. 0)goto 17251                                            zipper
      jerr=jerr-ilm                                                        zipper
      go to 12180                                                          zipper
17251 continue                                                             zipper
      ix=0                                                                 zipper
17260 do 17261 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 17261                           zipper
      ix=1                                                                 zipper
      goto 17262                                                           zipper
17261 continue                                                             zipper
17262 continue                                                             zipper
      if(ix .ne. 0)goto 17281                                              zipper
17290 do 17291 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 17291                                            zipper
      if(ju(k).eq.0)goto 17291                                             zipper
      ga(k)=abs(dot_product(wr,x(:,k)))                                    zipper
      if(ga(k) .le. sa*vp(k))goto 17311                                    zipper
      ixx(k)=1                                                             zipper
      ix=1                                                                 zipper
17311 continue                                                             zipper
17291 continue                                                             zipper
17292 continue                                                             zipper
      if(ix.eq.1) go to 10880                                              zipper
      goto 17032                                                           zipper
17281 continue                                                             zipper
      goto 17031                                                           zipper
17032 continue                                                             zipper
      if(nin .le. nx)goto 17331                                            zipper
      jerr=-10000-ilm                                                      zipper
      goto 16952                                                           zipper
17331 continue                                                             zipper
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))                               zipper
      kin(ilm)=nin                                                         zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      dev(ilm)=(risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr                   zipper
      if(ilm.lt.mnl)goto 16951                                             zipper
      if(flmin.ge.1.0)goto 16951                                           zipper
      me=0                                                                 zipper
17340 do 17341 j=1,nin                                                     zipper
      if(ao(j,ilm).ne.0.0) me=me+1                                         zipper
17341 continue                                                             zipper
17342 continue                                                             zipper
      if(me.gt.ne)goto 16952                                               zipper
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 16952              zipper
      if(dev(ilm).gt.devmax)goto 16952                                     zipper
16951 continue                                                             zipper
16952 continue                                                             zipper
      g=f                                                                  zipper
12180 continue                                                             zipper
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)              zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine cxmodval(ca,ia,nin,n,x,f)                                 zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision ca(nin),x(n,*),f(n)                                 zipper
      integer ia(nin)                                                      zipper
      f=0.0                                                                zipper
      if(nin.le.0) return                                                  zipper
17350 do 17351 i=1,n                                                       zipper
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                      zipper
17351 continue                                                             zipper
17352 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine groups(no,y,d,q,nk,kp,jp,t0,jerr)                         zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision y(no),d(no),q(no)                                   zipper
      integer jp(no),kp(*)                                                 zipper
17360 do 17361 j=1,no                                                      zipper
      jp(j)=j                                                              zipper
17361 continue                                                             zipper
17362 continue                                                             zipper
      call psort7(y,jp,1,no)                                               zipper
      nj=0                                                                 zipper
17370 do 17371 j=1,no                                                      zipper
      if(q(jp(j)).le.0.0)goto 17371                                        zipper
      nj=nj+1                                                              zipper
      jp(nj)=jp(j)                                                         zipper
17371 continue                                                             zipper
17372 continue                                                             zipper
      if(nj .ne. 0)goto 17391                                              zipper
      jerr=20000                                                           zipper
      return                                                               zipper
17391 continue                                                             zipper
      j=1                                                                  zipper
17400 continue                                                             zipper
17401 if(d(jp(j)).gt.0.0)goto 17402                                        zipper
      j=j+1                                                                zipper
      if(j.gt.nj)goto 17402                                                zipper
      goto 17401                                                           zipper
17402 continue                                                             zipper
      if(j .lt. nj-1)goto 17421                                            zipper
      jerr=30000                                                           zipper
      return                                                               zipper
17421 continue                                                             zipper
      t0=y(jp(j))                                                          zipper
      j0=j-1                                                               zipper
      if(j0 .le. 0)goto 17441                                              zipper
17450 continue                                                             zipper
17451 if(y(jp(j0)).lt.t0)goto 17452                                        zipper
      j0=j0-1                                                              zipper
      if(j0.eq.0)goto 17452                                                zipper
      goto 17451                                                           zipper
17452 continue                                                             zipper
      if(j0 .le. 0)goto 17471                                              zipper
      nj=nj-j0                                                             zipper
17480 do 17481 j=1,nj                                                      zipper
      jp(j)=jp(j+j0)                                                       zipper
17481 continue                                                             zipper
17482 continue                                                             zipper
17471 continue                                                             zipper
17441 continue                                                             zipper
      jerr=0                                                               zipper
      nk=0                                                                 zipper
      yk=t0                                                                zipper
      j=2                                                                  zipper
17490 continue                                                             zipper
17491 continue                                                             zipper
17500 continue                                                             zipper
17501 if(d(jp(j)).gt.0.0.and.y(jp(j)).gt.yk)goto 17502                     zipper
      j=j+1                                                                zipper
      if(j.gt.nj)goto 17502                                                zipper
      goto 17501                                                           zipper
17502 continue                                                             zipper
      nk=nk+1                                                              zipper
      kp(nk)=j-1                                                           zipper
      if(j.gt.nj)goto 17492                                                zipper
      if(j .ne. nj)goto 17521                                              zipper
      nk=nk+1                                                              zipper
      kp(nk)=nj                                                            zipper
      goto 17492                                                           zipper
17521 continue                                                             zipper
      yk=y(jp(j))                                                          zipper
      j=j+1                                                                zipper
      goto 17491                                                           zipper
17492 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine outer(no,nk,d,dk,kp,jp,e,wr,w,jerr,u)                     zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision d(no),dk(nk),wr(no),w(no)                           zipper
      double precision e(no),u(no),b,c                                     zipper
      integer kp(nk),jp(no)                                                zipper
      call usk(no,nk,kp,jp,e,u)                                            zipper
      b=dk(1)/u(1)                                                         zipper
      c=dk(1)/u(1)**2                                                      zipper
      jerr=0                                                               zipper
17530 do 17531 j=1,kp(1)                                                   zipper
      i=jp(j)                                                              zipper
      w(i)=e(i)*(b-e(i)*c)                                                 zipper
      if(w(i) .gt. 0.0)goto 17551                                          zipper
      jerr=-30000                                                          zipper
      return                                                               zipper
17551 continue                                                             zipper
      wr(i)=d(i)-e(i)*b                                                    zipper
17531 continue                                                             zipper
17532 continue                                                             zipper
17560 do 17561 k=2,nk                                                      zipper
      j1=kp(k-1)+1                                                         zipper
      j2=kp(k)                                                             zipper
      b=b+dk(k)/u(k)                                                       zipper
      c=c+dk(k)/u(k)**2                                                    zipper
17570 do 17571 j=j1,j2                                                     zipper
      i=jp(j)                                                              zipper
      w(i)=e(i)*(b-e(i)*c)                                                 zipper
      if(w(i) .gt. 0.0)goto 17591                                          zipper
      jerr=-30000                                                          zipper
      return                                                               zipper
17591 continue                                                             zipper
      wr(i)=d(i)-e(i)*b                                                    zipper
17571 continue                                                             zipper
17572 continue                                                             zipper
17561 continue                                                             zipper
17562 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine vars(no,ni,x,w,ixx,v)                                     zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),w(no),v(ni)                                zipper
      integer ixx(ni)                                                      zipper
17600 do 17601 j=1,ni                                                      zipper
      if(ixx(j).gt.0) v(j)=dot_product(w,x(:,j)**2)                        zipper
17601 continue                                                             zipper
17602 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine died(no,nk,d,kp,jp,dk)                                    zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision d(no),dk(nk)                                        zipper
      integer kp(nk),jp(no)                                                zipper
      dk(1)=sum(d(jp(1:kp(1))))                                            zipper
17610 do 17611 k=2,nk                                                      zipper
      dk(k)=sum(d(jp((kp(k-1)+1):kp(k))))                                  zipper
17611 continue                                                             zipper
17612 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine usk(no,nk,kp,jp,e,u)                                      zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision e(no),u(nk),h                                       zipper
      integer kp(nk),jp(no)                                                zipper
      h=0.0                                                                zipper
17620 do 17621 k=nk,1,-1                                                   zipper
      j2=kp(k)                                                             zipper
      j1=1                                                                 zipper
      if(k.gt.1) j1=kp(k-1)+1                                              zipper
17630 do 17631 j=j2,j1,-1                                                  zipper
      h=h+e(jp(j))                                                         zipper
17631 continue                                                             zipper
17632 continue                                                             zipper
      u(k)=h                                                               zipper
17621 continue                                                             zipper
17622 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      function risk(no,ni,nk,d,dk,f,e,kp,jp,u)                             zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision d(no),dk(nk),f(no)                                  zipper
      integer kp(nk),jp(no)                                                zipper
      double precision e(no),u(nk),s                                       zipper
      call usk(no,nk,kp,jp,e,u)                                            zipper
      u=log(u)                                                             zipper
      risk=dot_product(d,f)-dot_product(dk,u)                              zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine loglike(no,ni,x,y,d,g,w,nlam,a,flog,jerr)                 zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no),d(no),g(no),w(no),a(ni,nlam),flog(   zipper
     *nlam)
      double precision, dimension (:), allocatable :: dk,f,xm,dq,q              
      double precision, dimension (:), allocatable :: e,uu                      
      integer, dimension (:), allocatable :: jp,kp                              
      allocate(e(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) go to 12180                                            zipper
      allocate(q(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) go to 12180                                            zipper
      allocate(uu(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) go to 12180                                            zipper
      allocate(f(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) go to 12180                                            zipper
      allocate(dk(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) go to 12180                                            zipper
      allocate(jp(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) go to 12180                                            zipper
      allocate(kp(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) go to 12180                                            zipper
      allocate(dq(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) go to 12180                                            zipper
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) go to 12180                                            zipper
      q=max(0d0,w)                                                         zipper
      sw=sum(q)                                                            zipper
      if(sw .gt. 0.0)goto 17651                                            zipper
      jerr=9999                                                            zipper
      go to 12180                                                          zipper
17651 continue                                                             zipper
      call groups(no,y,d,q,nk,kp,jp,t0,jerr)                               zipper
      if(jerr.ne.0) go to 12180                                            zipper
      fmax=log(huge(e(1))*0.1)                                             zipper
      dq=d*q                                                               zipper
      call died(no,nk,dq,kp,jp,dk)                                         zipper
      gm=dot_product(q,g)/sw                                               zipper
17660 do 17661 j=1,ni                                                      zipper
      xm(j)=dot_product(q,x(:,j))/sw                                       zipper
17661 continue                                                             zipper
17662 continue                                                             zipper
17670 do 17671 lam=1,nlam                                                  zipper
17680 do 17681 i=1,no                                                      zipper
      f(i)=g(i)-gm+dot_product(a(:,lam),(x(i,:)-xm))                       zipper
      e(i)=q(i)*exp(sign(min(abs(f(i)),fmax),f(i)))                        zipper
17681 continue                                                             zipper
17682 continue                                                             zipper
      flog(lam)=risk(no,ni,nk,dq,dk,f,e,kp,jp,uu)                          zipper
17671 continue                                                             zipper
17672 continue                                                             zipper
12180 continue                                                             zipper
      deallocate(e,uu,dk,f,jp,kp,dq)                                       zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine fishnet (parm,no,ni,x,y,g,w,jd,vp,cl,ne,nx,nlam,flmin,u   zipper
     *lam,thr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no),g(no),w(no),vp(ni),ulam(nlam)        zipper
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   zipper
      integer jd(*),ia(nx),nin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 17701                                    zipper
      jerr=10000                                                           zipper
      return                                                               zipper
17701 continue                                                             zipper
      if(minval(y) .ge. 0.0)goto 17721                                     zipper
      jerr=8888                                                            zipper
      return                                                               zipper
17721 continue                                                             zipper
      allocate(ww(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(vq(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      if(isd .le. 0)goto 17741                                             zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
17741 continue                                                             zipper
      call chkvars(no,ni,x,ju)                                             zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 17761                                      zipper
      jerr=7777                                                            zipper
      go to 12180                                                          zipper
17761 continue                                                             zipper
      vq=max(0d0,vp)                                                       zipper
      vq=vq*ni/sum(vq)                                                     zipper
      ww=max(0d0,w)                                                        zipper
      sw=sum(ww)                                                           zipper
      if(sw .gt. 0.0)goto 17781                                            zipper
      jerr=9999                                                            zipper
      go to 12180                                                          zipper
17781 continue                                                             zipper
      ww=ww/sw                                                             zipper
      call lstandard1(no,ni,x,ww,ju,isd,intr,xm,xs)                        zipper
      if(isd .le. 0)goto 17801                                             zipper
17810 do 17811 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
17811 continue                                                             zipper
17812 continue                                                             zipper
17801 continue                                                             zipper
      call fishnet1(parm,no,ni,x,y,g,ww,ju,vq,cl,ne,nx,nlam,flmin,ulam,t   zipper
     *hr,  isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr)
      if(jerr.gt.0) go to 12180                                            zipper
      dev0=2.0*sw*dev0                                                     zipper
17820 do 17821 k=1,lmu                                                     zipper
      nk=nin(k)                                                            zipper
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      zipper
      if(intr .ne. 0)goto 17841                                            zipper
      a0(k)=0.0                                                            zipper
      goto 17851                                                           zipper
17841 continue                                                             zipper
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     zipper
17851 continue                                                             zipper
17831 continue                                                             zipper
17821 continue                                                             zipper
17822 continue                                                             zipper
12180 continue                                                             zipper
      deallocate(ww,ju,vq,xm)                                              zipper
      if(isd.gt.0) deallocate(xs)                                          zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine fishnet1(parm,no,ni,x,y,g,q,ju,vp,cl,ne,nx,nlam,flmin,u   zipper
     *lam,shri,  isd,intr,maxit,lmu,a0,ca,m,kin,dev0,dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no),g(no),q(no),vp(ni),ulam(nlam)        zipper
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   zipper
      integer ju(ni),m(nx),kin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: t,w,wr,v,a,f,as,ga        
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      sml=sml*10.0                                                         zipper
      allocate(a(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(as(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(t(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ixx(1:ni),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(wr(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(v(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(w(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(f(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      bta=parm                                                             zipper
      omb=1.0-bta                                                          zipper
      t=q*y                                                                zipper
      yb=sum(t)                                                            zipper
      fmax=log(huge(bta)*0.1)                                              zipper
      if(nonzero(no,g) .ne. 0)goto 17871                                   zipper
      if(intr .eq. 0)goto 17891                                            zipper
      w=q*yb                                                               zipper
      az=log(yb)                                                           zipper
      f=az                                                                 zipper
      dv0=yb*(az-1.0)                                                      zipper
      goto 17901                                                           zipper
17891 continue                                                             zipper
      w=q                                                                  zipper
      az=0.0                                                               zipper
      f=az                                                                 zipper
      dv0=-1.0                                                             zipper
17901 continue                                                             zipper
17881 continue                                                             zipper
      goto 17911                                                           zipper
17871 continue                                                             zipper
      w=q*exp(sign(min(abs(g),fmax),g))                                    zipper
      v0=sum(w)                                                            zipper
      if(intr .eq. 0)goto 17931                                            zipper
      eaz=yb/v0                                                            zipper
      w=eaz*w                                                              zipper
      az=log(eaz)                                                          zipper
      f=az+g                                                               zipper
      dv0=dot_product(t,g)-yb*(1.0-az)                                     zipper
      goto 17941                                                           zipper
17931 continue                                                             zipper
      az=0.0                                                               zipper
      f=g                                                                  zipper
      dv0=dot_product(t,g)-v0                                              zipper
17941 continue                                                             zipper
17921 continue                                                             zipper
17911 continue                                                             zipper
17861 continue                                                             zipper
      a=0.0                                                                zipper
      as=0.0                                                               zipper
      wr=t-w                                                               zipper
      v0=1.0                                                               zipper
      if(intr.ne.0) v0=yb                                                  zipper
      dvr=-yb                                                              zipper
17950 do 17951 i=1,no                                                      zipper
      if(t(i).gt.0.0) dvr=dvr+t(i)*log(y(i))                               zipper
17951 continue                                                             zipper
17952 continue                                                             zipper
      dvr=dvr-dv0                                                          zipper
      dev0=dvr                                                             zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 17971                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
17971 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      mnl=min(mnlam,nlam)                                                  zipper
      shr=shri*dev0                                                        zipper
      ixx=0                                                                zipper
      al=0.0                                                               zipper
17980 do 17981 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 17981                                             zipper
      ga(j)=abs(dot_product(wr,x(:,j)))                                    zipper
17981 continue                                                             zipper
17982 continue                                                             zipper
17990 do 17991 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 18011                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 18001                                                           zipper
18011 if(ilm .le. 2)goto 18021                                             zipper
      al=al*alf                                                            zipper
      goto 18001                                                           zipper
18021 if(ilm .ne. 1)goto 18031                                             zipper
      al=big                                                               zipper
      goto 18041                                                           zipper
18031 continue                                                             zipper
      al0=0.0                                                              zipper
18050 do 18051 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 18051                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
18051 continue                                                             zipper
18052 continue                                                             zipper
      al0=al0/max(bta,1.0d-3)                                              zipper
      al=alf*al0                                                           zipper
18041 continue                                                             zipper
18001 continue                                                             zipper
      al2=al*omb                                                           zipper
      al1=al*bta                                                           zipper
      tlam=bta*(2.0*al-al0)                                                zipper
18060 do 18061 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 18061                                            zipper
      if(ju(k).eq.0)goto 18061                                             zipper
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     zipper
18061 continue                                                             zipper
18062 continue                                                             zipper
10880 continue                                                             zipper
18070 continue                                                             zipper
18071 continue                                                             zipper
      az0=az                                                               zipper
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                zipper
18080 do 18081 j=1,ni                                                      zipper
      if(ixx(j).ne.0) v(j)=dot_product(w,x(:,j)**2)                        zipper
18081 continue                                                             zipper
18082 continue                                                             zipper
18090 continue                                                             zipper
18091 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
18100 do 18101 k=1,ni                                                      zipper
      if(ixx(k).eq.0)goto 18101                                            zipper
      ak=a(k)                                                              zipper
      u=dot_product(wr,x(:,k))+v(k)*ak                                     zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 18121                                            zipper
      a(k)=0.0                                                             zipper
      goto 18131                                                           zipper
18121 continue                                                             zipper
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           zipper
18131 continue                                                             zipper
18111 continue                                                             zipper
      if(a(k).eq.ak)goto 18101                                             zipper
      d=a(k)-ak                                                            zipper
      dlx=max(dlx,v(k)*d**2)                                               zipper
      wr=wr-d*w*x(:,k)                                                     zipper
      f=f+d*x(:,k)                                                         zipper
      if(mm(k) .ne. 0)goto 18151                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 18102                                              zipper
      mm(k)=nin                                                            zipper
      m(nin)=k                                                             zipper
18151 continue                                                             zipper
18101 continue                                                             zipper
18102 continue                                                             zipper
      if(nin.gt.nx)goto 18092                                              zipper
      if(intr .eq. 0)goto 18171                                            zipper
      d=sum(wr)/v0                                                         zipper
      az=az+d                                                              zipper
      dlx=max(dlx,v0*d**2)                                                 zipper
      wr=wr-d*w                                                            zipper
      f=f+d                                                                zipper
18171 continue                                                             zipper
      if(dlx.lt.shr)goto 18092                                             zipper
      if(nlp .le. maxit)goto 18191                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
18191 continue                                                             zipper
18200 continue                                                             zipper
18201 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
18210 do 18211 l=1,nin                                                     zipper
      k=m(l)                                                               zipper
      ak=a(k)                                                              zipper
      u=dot_product(wr,x(:,k))+v(k)*ak                                     zipper
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 18231                                            zipper
      a(k)=0.0                                                             zipper
      goto 18241                                                           zipper
18231 continue                                                             zipper
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           zipper
18241 continue                                                             zipper
18221 continue                                                             zipper
      if(a(k).eq.ak)goto 18211                                             zipper
      d=a(k)-ak                                                            zipper
      dlx=max(dlx,v(k)*d**2)                                               zipper
      wr=wr-d*w*x(:,k)                                                     zipper
      f=f+d*x(:,k)                                                         zipper
18211 continue                                                             zipper
18212 continue                                                             zipper
      if(intr .eq. 0)goto 18261                                            zipper
      d=sum(wr)/v0                                                         zipper
      az=az+d                                                              zipper
      dlx=max(dlx,v0*d**2)                                                 zipper
      wr=wr-d*w                                                            zipper
      f=f+d                                                                zipper
18261 continue                                                             zipper
      if(dlx.lt.shr)goto 18202                                             zipper
      if(nlp .le. maxit)goto 18281                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
18281 continue                                                             zipper
      goto 18201                                                           zipper
18202 continue                                                             zipper
      goto 18091                                                           zipper
18092 continue                                                             zipper
      if(nin.gt.nx)goto 18072                                              zipper
      w=q*exp(sign(min(abs(f),fmax),f))                                    zipper
      v0=sum(w)                                                            zipper
      wr=t-w                                                               zipper
      if(v0*(az-az0)**2 .ge. shr)goto 18301                                zipper
      ix=0                                                                 zipper
18310 do 18311 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 18311                            zipper
      ix=1                                                                 zipper
      goto 18312                                                           zipper
18311 continue                                                             zipper
18312 continue                                                             zipper
      if(ix .ne. 0)goto 18331                                              zipper
18340 do 18341 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 18341                                            zipper
      if(ju(k).eq.0)goto 18341                                             zipper
      ga(k)=abs(dot_product(wr,x(:,k)))                                    zipper
      if(ga(k) .le. al1*vp(k))goto 18361                                   zipper
      ixx(k)=1                                                             zipper
      ix=1                                                                 zipper
18361 continue                                                             zipper
18341 continue                                                             zipper
18342 continue                                                             zipper
      if(ix.eq.1) go to 10880                                              zipper
      goto 18072                                                           zipper
18331 continue                                                             zipper
18301 continue                                                             zipper
      goto 18071                                                           zipper
18072 continue                                                             zipper
      if(nin .le. nx)goto 18381                                            zipper
      jerr=-10000-ilm                                                      zipper
      goto 17992                                                           zipper
18381 continue                                                             zipper
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               zipper
      kin(ilm)=nin                                                         zipper
      a0(ilm)=az                                                           zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      dev(ilm)=(dot_product(t,f)-v0-dv0)/dvr                               zipper
      if(ilm.lt.mnl)goto 17991                                             zipper
      if(flmin.ge.1.0)goto 17991                                           zipper
      me=0                                                                 zipper
18390 do 18391 j=1,nin                                                     zipper
      if(ca(j,ilm).ne.0.0) me=me+1                                         zipper
18391 continue                                                             zipper
18392 continue                                                             zipper
      if(me.gt.ne)goto 17992                                               zipper
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 17992              zipper
      if(dev(ilm).gt.devmax)goto 17992                                     zipper
17991 continue                                                             zipper
17992 continue                                                             zipper
      g=f                                                                  zipper
12180 continue                                                             zipper
      deallocate(t,w,wr,v,a,f,as,mm,ga,ixx)                                zipper
      return                                                               zipper
      end                                                                  zipper
      function nonzero(n,v)                                                zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision v(n)                                                zipper
      nonzero=0                                                            zipper
18400 do 18401 i=1,n                                                       zipper
      if(v(i) .eq. 0.0)goto 18421                                          zipper
      nonzero=1                                                            zipper
      return                                                               zipper
18421 continue                                                             zipper
18401 continue                                                             zipper
18402 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                               zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision a(nx,lmu),b(ni,lmu)                                 zipper
      integer ia(nx),nin(lmu)                                              zipper
18430 do 18431 lam=1,lmu                                                   zipper
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                        zipper
18431 continue                                                             zipper
18432 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine lsolns(ni,nx,nc,lmu,a,ia,nin,b)                           zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision a(nx,nc,lmu),b(ni,nc,lmu)                           zipper
      integer ia(nx),nin(lmu)                                              zipper
18440 do 18441 lam=1,lmu                                                   zipper
      call luncomp(ni,nx,nc,a(1,1,lam),ia,nin(lam),b(1,1,lam))             zipper
18441 continue                                                             zipper
18442 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine deviance(no,ni,x,y,g,q,nlam,a0,a,flog,jerr)               zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no),g(no),q(no),a(ni,nlam),a0(nlam),fl   zipper
     *og(nlam)
      double precision, dimension (:), allocatable :: w                         
      if(minval(y) .ge. 0.0)goto 18461                                     zipper
      jerr=8888                                                            zipper
      return                                                               zipper
18461 continue                                                             zipper
      allocate(w(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      w=max(0d0,q)                                                         zipper
      sw=sum(w)                                                            zipper
      if(sw .gt. 0.0)goto 18481                                            zipper
      jerr=9999                                                            zipper
      go to 12180                                                          zipper
18481 continue                                                             zipper
      yb=dot_product(w,y)/sw                                               zipper
      fmax=log(huge(y(1))*0.1)                                             zipper
18490 do 18491 lam=1,nlam                                                  zipper
      s=0.0                                                                zipper
18500 do 18501 i=1,no                                                      zipper
      if(w(i).le.0.0)goto 18501                                            zipper
      f=g(i)+a0(lam)+dot_product(a(:,lam),x(i,:))                          zipper
      s=s+w(i)*(y(i)*f-exp(sign(min(abs(f),fmax),f)))                      zipper
18501 continue                                                             zipper
18502 continue                                                             zipper
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                zipper
18491 continue                                                             zipper
18492 continue                                                             zipper
12180 continue                                                             zipper
      deallocate(w)                                                        zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spfishnet (parm,no,ni,x,ix,jx,y,g,w,jd,vp,cl,ne,nx,nlam   zipper
     *,flmin,  ulam,thr,isd,intr,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp
     *,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),g(no),w(no),vp(ni),ulam(nlam),cl(2,ni)   zipper
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)            zipper
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           zipper
      double precision, dimension (:), allocatable :: xm,xs,ww,vq               
      integer, dimension (:), allocatable :: ju                                 
      if(maxval(vp) .gt. 0.0)goto 18521                                    zipper
      jerr=10000                                                           zipper
      return                                                               zipper
18521 continue                                                             zipper
      if(minval(y) .ge. 0.0)goto 18541                                     zipper
      jerr=8888                                                            zipper
      return                                                               zipper
18541 continue                                                             zipper
      allocate(ww(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(vq(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      call spchkvars(no,ni,x,ix,ju)                                        zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 18561                                      zipper
      jerr=7777                                                            zipper
      go to 12180                                                          zipper
18561 continue                                                             zipper
      vq=max(0d0,vp)                                                       zipper
      vq=vq*ni/sum(vq)                                                     zipper
      ww=max(0d0,w)                                                        zipper
      sw=sum(ww)                                                           zipper
      if(sw .gt. 0.0)goto 18581                                            zipper
      jerr=9999                                                            zipper
      go to 12180                                                          zipper
18581 continue                                                             zipper
      ww=ww/sw                                                             zipper
      call splstandard2(no,ni,x,ix,jx,ww,ju,isd,intr,xm,xs)                zipper
      if(isd .le. 0)goto 18601                                             zipper
18610 do 18611 j=1,ni                                                      zipper
      cl(:,j)=cl(:,j)*xs(j)                                                zipper
18611 continue                                                             zipper
18612 continue                                                             zipper
18601 continue                                                             zipper
      call spfishnet1(parm,no,ni,x,ix,jx,y,g,ww,ju,vq,cl,ne,nx,nlam,flmi   zipper
     *n,ulam,thr,  isd,intr,maxit,xm,xs,lmu,a0,ca,ia,nin,dev0,dev,alm,nl
     *p,jerr)
      if(jerr.gt.0) go to 12180                                            zipper
      dev0=2.0*sw*dev0                                                     zipper
18620 do 18621 k=1,lmu                                                     zipper
      nk=nin(k)                                                            zipper
      if(isd.gt.0) ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))                      zipper
      if(intr .ne. 0)goto 18641                                            zipper
      a0(k)=0.0                                                            zipper
      goto 18651                                                           zipper
18641 continue                                                             zipper
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))                     zipper
18651 continue                                                             zipper
18631 continue                                                             zipper
18621 continue                                                             zipper
18622 continue                                                             zipper
12180 continue                                                             zipper
      deallocate(ww,ju,vq,xm,xs)                                           zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spfishnet1(parm,no,ni,x,ix,jx,y,g,q,ju,vp,cl,ne,nx,nlam   zipper
     *,flmin,ulam,  shri,isd,intr,maxit,xb,xs,lmu,a0,ca,m,kin,dev0,dev,a
     *lm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),g(no),q(no),vp(ni),ulam(nlam),xb(ni),x   zipper
     *s(ni)
      double precision ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam),cl(2,ni)   zipper
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           zipper
      double precision, dimension (:), allocatable :: qy,t,w,wr,v               
      double precision, dimension (:), allocatable :: a,as,xm,ga                
      integer, dimension (:), allocatable :: mm,ixx                             
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      sml=sml*10.0                                                         zipper
      allocate(a(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(as(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(t(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ixx(1:ni),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(wr(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(v(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(w(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(qy(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      bta=parm                                                             zipper
      omb=1.0-bta                                                          zipper
      fmax=log(huge(bta)*0.1)                                              zipper
      qy=q*y                                                               zipper
      yb=sum(qy)                                                           zipper
      if(nonzero(no,g) .ne. 0)goto 18671                                   zipper
      t=0.0                                                                zipper
      if(intr .eq. 0)goto 18691                                            zipper
      w=q*yb                                                               zipper
      az=log(yb)                                                           zipper
      uu=az                                                                zipper
      xm=yb*xb                                                             zipper
      dv0=yb*(az-1.0)                                                      zipper
      goto 18701                                                           zipper
18691 continue                                                             zipper
      w=q                                                                  zipper
      xm=0.0                                                               zipper
      uu=0.0                                                               zipper
      az=uu                                                                zipper
      dv0=-1.0                                                             zipper
18701 continue                                                             zipper
18681 continue                                                             zipper
      goto 18711                                                           zipper
18671 continue                                                             zipper
      w=q*exp(sign(min(abs(g),fmax),g))                                    zipper
      ww=sum(w)                                                            zipper
      t=g                                                                  zipper
      if(intr .eq. 0)goto 18731                                            zipper
      eaz=yb/ww                                                            zipper
      w=eaz*w                                                              zipper
      az=log(eaz)                                                          zipper
      uu=az                                                                zipper
      dv0=dot_product(qy,g)-yb*(1.0-az)                                    zipper
      goto 18741                                                           zipper
18731 continue                                                             zipper
      uu=0.0                                                               zipper
      az=uu                                                                zipper
      dv0=dot_product(qy,g)-ww                                             zipper
18741 continue                                                             zipper
18721 continue                                                             zipper
18750 do 18751 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 18751                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             zipper
18751 continue                                                             zipper
18752 continue                                                             zipper
18711 continue                                                             zipper
18661 continue                                                             zipper
      tt=yb*uu                                                             zipper
      ww=1.0                                                               zipper
      if(intr.ne.0) ww=yb                                                  zipper
      wr=qy-q*(yb*(1.0-uu))                                                zipper
      a=0.0                                                                zipper
      as=0.0                                                               zipper
      dvr=-yb                                                              zipper
18760 do 18761 i=1,no                                                      zipper
      if(qy(i).gt.0.0) dvr=dvr+qy(i)*log(y(i))                             zipper
18761 continue                                                             zipper
18762 continue                                                             zipper
      dvr=dvr-dv0                                                          zipper
      dev0=dvr                                                             zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 18781                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
18781 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      mnl=min(mnlam,nlam)                                                  zipper
      shr=shri*dev0                                                        zipper
      al=0.0                                                               zipper
      ixx=0                                                                zipper
18790 do 18791 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 18791                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   zipper
     *)-xb(j)*tt)/xs(j)
18791 continue                                                             zipper
18792 continue                                                             zipper
18800 do 18801 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 18821                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 18811                                                           zipper
18821 if(ilm .le. 2)goto 18831                                             zipper
      al=al*alf                                                            zipper
      goto 18811                                                           zipper
18831 if(ilm .ne. 1)goto 18841                                             zipper
      al=big                                                               zipper
      goto 18851                                                           zipper
18841 continue                                                             zipper
      al0=0.0                                                              zipper
18860 do 18861 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 18861                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
18861 continue                                                             zipper
18862 continue                                                             zipper
      al0=al0/max(bta,1.0d-3)                                              zipper
      al=alf*al0                                                           zipper
18851 continue                                                             zipper
18811 continue                                                             zipper
      al2=al*omb                                                           zipper
      al1=al*bta                                                           zipper
      tlam=bta*(2.0*al-al0)                                                zipper
18870 do 18871 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 18871                                            zipper
      if(ju(k).eq.0)goto 18871                                             zipper
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     zipper
18871 continue                                                             zipper
18872 continue                                                             zipper
10880 continue                                                             zipper
18880 continue                                                             zipper
18881 continue                                                             zipper
      az0=az                                                               zipper
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))                                zipper
18890 do 18891 j=1,ni                                                      zipper
      if(ixx(j).eq.0)goto 18891                                            zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             zipper
      v(j)=(dot_product(w(jx(jb:je)),x(jb:je)**2)  -2.0*xb(j)*xm(j)+ww*x   zipper
     *b(j)**2)/xs(j)**zipper
18891 continue                                                             zipper
18892 continue                                                             zipper
18900 continue                                                             zipper
18901 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
18910 do 18911 k=1,ni                                                      zipper
      if(ixx(k).eq.0)goto 18911                                            zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      ak=a(k)                                                              zipper
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   zipper
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 18931                                            zipper
      a(k)=0.0                                                             zipper
      goto 18941                                                           zipper
18931 continue                                                             zipper
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           zipper
18941 continue                                                             zipper
18921 continue                                                             zipper
      if(a(k).eq.ak)goto 18911                                             zipper
      if(mm(k) .ne. 0)goto 18961                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 18912                                              zipper
      mm(k)=nin                                                            zipper
      m(nin)=k                                                             zipper
18961 continue                                                             zipper
      d=a(k)-ak                                                            zipper
      dlx=max(dlx,v(k)*d**2)                                               zipper
      dv=d/xs(k)                                                           zipper
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 zipper
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                zipper
      uu=uu-dv*xb(k)                                                       zipper
      tt=tt-dv*xm(k)                                                       zipper
18911 continue                                                             zipper
18912 continue                                                             zipper
      if(nin.gt.nx)goto 18902                                              zipper
      if(intr .eq. 0)goto 18981                                            zipper
      d=tt/ww-uu                                                           zipper
      az=az+d                                                              zipper
      dlx=max(dlx,ww*d**2)                                                 zipper
      uu=uu+d                                                              zipper
18981 continue                                                             zipper
      if(dlx.lt.shr)goto 18902                                             zipper
      if(nlp .le. maxit)goto 19001                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
19001 continue                                                             zipper
19010 continue                                                             zipper
19011 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
19020 do 19021 l=1,nin                                                     zipper
      k=m(l)                                                               zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      ak=a(k)                                                              zipper
      u=(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(k)-ww*xb(k))-xb(k)   zipper
     **tt)/xs(k)+v(k)*ak
      au=abs(u)-vp(k)*al1                                                  zipper
      if(au .gt. 0.0)goto 19041                                            zipper
      a(k)=0.0                                                             zipper
      goto 19051                                                           zipper
19041 continue                                                             zipper
      a(k)=max(cl(1,k),min(cl(2,k),sign(au,u)/(v(k)+vp(k)*al2)))           zipper
19051 continue                                                             zipper
19031 continue                                                             zipper
      if(a(k).eq.ak)goto 19021                                             zipper
      d=a(k)-ak                                                            zipper
      dlx=max(dlx,v(k)*d**2)                                               zipper
      dv=d/xs(k)                                                           zipper
      wr(jx(jb:je))=wr(jx(jb:je))-dv*w(jx(jb:je))*x(jb:je)                 zipper
      t(jx(jb:je))=t(jx(jb:je))+dv*x(jb:je)                                zipper
      uu=uu-dv*xb(k)                                                       zipper
      tt=tt-dv*xm(k)                                                       zipper
19021 continue                                                             zipper
19022 continue                                                             zipper
      if(intr .eq. 0)goto 19071                                            zipper
      d=tt/ww-uu                                                           zipper
      az=az+d                                                              zipper
      dlx=max(dlx,ww*d**2)                                                 zipper
      uu=uu+d                                                              zipper
19071 continue                                                             zipper
      if(dlx.lt.shr)goto 19012                                             zipper
      if(nlp .le. maxit)goto 19091                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
19091 continue                                                             zipper
      goto 19011                                                           zipper
19012 continue                                                             zipper
      goto 18901                                                           zipper
18902 continue                                                             zipper
      if(nin.gt.nx)goto 18882                                              zipper
      euu=exp(sign(min(abs(uu),fmax),uu))                                  zipper
      w=euu*q*exp(sign(min(abs(t),fmax),t))                                zipper
      ww=sum(w)                                                            zipper
      wr=qy-w*(1.0-uu)                                                     zipper
      tt=sum(wr)                                                           zipper
      if(ww*(az-az0)**2 .ge. shr)goto 19111                                zipper
      kx=0                                                                 zipper
19120 do 19121 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(v(k)*(a(k)-as(k))**2.lt.shr)goto 19121                            zipper
      kx=1                                                                 zipper
      goto 19122                                                           zipper
19121 continue                                                             zipper
19122 continue                                                             zipper
      if(kx .ne. 0)goto 19141                                              zipper
19150 do 19151 j=1,ni                                                      zipper
      if(ixx(j).eq.1)goto 19151                                            zipper
      if(ju(j).eq.0)goto 19151                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             zipper
      ga(j)=abs(dot_product(wr(jx(jb:je)),x(jb:je))  -uu*(xm(j)-ww*xb(j)   zipper
     *)-xb(j)*tt)/xs(j)
      if(ga(j) .le. al1*vp(j))goto 19171                                   zipper
      ixx(j)=1                                                             zipper
      kx=1                                                                 zipper
19171 continue                                                             zipper
19151 continue                                                             zipper
19152 continue                                                             zipper
      if(kx.eq.1) go to 10880                                              zipper
      goto 18882                                                           zipper
19141 continue                                                             zipper
19111 continue                                                             zipper
      goto 18881                                                           zipper
18882 continue                                                             zipper
      if(nin .le. nx)goto 19191                                            zipper
      jerr=-10000-ilm                                                      zipper
      goto 18802                                                           zipper
19191 continue                                                             zipper
      if(nin.gt.0) ca(1:nin,ilm)=a(m(1:nin))                               zipper
      kin(ilm)=nin                                                         zipper
      a0(ilm)=az                                                           zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      dev(ilm)=(dot_product(qy,t)+yb*uu-ww-dv0)/dvr                        zipper
      if(ilm.lt.mnl)goto 18801                                             zipper
      if(flmin.ge.1.0)goto 18801                                           zipper
      me=0                                                                 zipper
19200 do 19201 j=1,nin                                                     zipper
      if(ca(j,ilm).ne.0.0) me=me+1                                         zipper
19201 continue                                                             zipper
19202 continue                                                             zipper
      if(me.gt.ne)goto 18802                                               zipper
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 18802              zipper
      if(dev(ilm).gt.devmax)goto 18802                                     zipper
18801 continue                                                             zipper
18802 continue                                                             zipper
      g=t+uu                                                               zipper
12180 continue                                                             zipper
      deallocate(t,w,wr,v,a,qy,xm,as,mm,ga,ixx)                            zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine spdeviance(no,ni,x,ix,jx,y,g,q,nlam,a0,a,flog,jerr)       zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),g(no),q(no),a(ni,nlam),a0(nlam),flog(n   zipper
     *lam)
      integer ix(*),jx(*)                                                  zipper
      double precision, dimension (:), allocatable :: w,f                       
      if(minval(y) .ge. 0.0)goto 19221                                     zipper
      jerr=8888                                                            zipper
      return                                                               zipper
19221 continue                                                             zipper
      allocate(w(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(f(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      w=max(0d0,q)                                                         zipper
      sw=sum(w)                                                            zipper
      if(sw .gt. 0.0)goto 19241                                            zipper
      jerr=9999                                                            zipper
      go to 12180                                                          zipper
19241 continue                                                             zipper
      yb=dot_product(w,y)/sw                                               zipper
      fmax=log(huge(y(1))*0.1)                                             zipper
19250 do 19251 lam=1,nlam                                                  zipper
      f=a0(lam)                                                            zipper
19260 do 19261 j=1,ni                                                      zipper
      if(a(j,lam).eq.0.0)goto 19261                                        zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      f(jx(jb:je))=f(jx(jb:je))+a(j,lam)*x(jb:je)                          zipper
19261 continue                                                             zipper
19262 continue                                                             zipper
      f=f+g                                                                zipper
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   zipper
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                zipper
19251 continue                                                             zipper
19252 continue                                                             zipper
12180 continue                                                             zipper
      deallocate(w,f)                                                      zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine cspdeviance(no,x,ix,jx,y,g,q,nx,nlam,a0,ca,ia,nin,flog,   zipper
     *jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no),g(no),q(no),ca(nx,nlam),a0(nlam),flog(   zipper
     *nlam)
      integer ix(*),jx(*),nin(nlam),ia(nx)                                 zipper
      double precision, dimension (:), allocatable :: w,f                       
      if(minval(y) .ge. 0.0)goto 19281                                     zipper
      jerr=8888                                                            zipper
      return                                                               zipper
19281 continue                                                             zipper
      allocate(w(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(f(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      w=max(0d0,q)                                                         zipper
      sw=sum(w)                                                            zipper
      if(sw .gt. 0.0)goto 19301                                            zipper
      jerr=9999                                                            zipper
      go to 12180                                                          zipper
19301 continue                                                             zipper
      yb=dot_product(w,y)/sw                                               zipper
      fmax=log(huge(y(1))*0.1)                                             zipper
19310 do 19311 lam=1,nlam                                                  zipper
      f=a0(lam)                                                            zipper
19320 do 19321 k=1,nin(lam)                                                zipper
      j=ia(k)                                                              zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      f(jx(jb:je))=f(jx(jb:je))+ca(k,lam)*x(jb:je)                         zipper
19321 continue                                                             zipper
19322 continue                                                             zipper
      f=f+g                                                                zipper
      s=dot_product(w,y*f-exp(sign(min(abs(f),fmax),f)))                   zipper
      flog(lam)=2.0*(sw*yb*(log(yb)-1.0)-s)                                zipper
19311 continue                                                             zipper
19312 continue                                                             zipper
12180 continue                                                             zipper
      deallocate(w,f)                                                      zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multelnet  (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   zipper
     *in,ulam,thr,isd,jsd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no,nr),w(no),vp(ni),ca(nx,nr,nlam)       zipper
      double precision ulam(nlam),a0(nr,nlam),rsq(nlam),alm(nlam),cl(2,n   zipper
     *i)
      integer jd(*),ia(nx),nin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 19341                                    zipper
      jerr=10000                                                           zipper
      return                                                               zipper
19341 continue                                                             zipper
      allocate(vq(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      vq=max(0d0,vp)                                                       zipper
      vq=vq*ni/sum(vq)                                                     zipper
      call multelnetn(parm,no,ni,nr,x,y,w,jd,vq,cl,ne,nx,nlam,flmin,ulam   zipper
     *,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                       zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multelnetn (parm,no,ni,nr,x,y,w,jd,vp,cl,ne,nx,nlam,flm   zipper
     *in,ulam,thr,  isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr
     *)
      implicit double precision(a-h,o-z)                                   zipper
      double precision vp(ni),x(no,ni),y(no,nr),w(no),ulam(nlam),cl(2,ni   zipper
     *)
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      zipper
      integer jd(*),ia(nx),nin(nlam)                                       zipper
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys            
      integer, dimension (:), allocatable :: ju                                 
      double precision, dimension (:,:,:), allocatable :: clt                   
      allocate(clt(1:2,1:nr,1:ni),stat=jerr);                                   
      if(jerr.ne.0) return                                                      
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ym(1:nr),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ys(1:nr),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xv(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      call chkvars(no,ni,x,ju)                                             zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 19361                                      zipper
      jerr=7777                                                            zipper
      return                                                               zipper
19361 continue                                                             zipper
      call multstandard1(no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym,ys,xv,y   zipper
     *s0,jerr)
      if(jerr.ne.0) return                                                 zipper
19370 do 19371 j=1,ni                                                      zipper
19380 do 19381 k=1,nr                                                      zipper
19390 do 19391 i=1,2                                                       zipper
      clt(i,k,j)=cl(i,j)                                                   zipper
19391 continue                                                             zipper
19392 continue                                                             zipper
19381 continue                                                             zipper
19382 continue                                                             zipper
19371 continue                                                             zipper
19372 continue                                                             zipper
      if(isd .le. 0)goto 19411                                             zipper
19420 do 19421 j=1,ni                                                      zipper
19430 do 19431 k=1,nr                                                      zipper
19440 do 19441 i=1,2                                                       zipper
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          zipper
19441 continue                                                             zipper
19442 continue                                                             zipper
19431 continue                                                             zipper
19432 continue                                                             zipper
19421 continue                                                             zipper
19422 continue                                                             zipper
19411 continue                                                             zipper
      if(jsd .le. 0)goto 19461                                             zipper
19470 do 19471 j=1,ni                                                      zipper
19480 do 19481 k=1,nr                                                      zipper
19490 do 19491 i=1,2                                                       zipper
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          zipper
19491 continue                                                             zipper
19492 continue                                                             zipper
19481 continue                                                             zipper
19482 continue                                                             zipper
19471 continue                                                             zipper
19472 continue                                                             zipper
19461 continue                                                             zipper
      call multelnet2(parm,ni,nr,ju,vp,clt,y,no,ne,nx,x,nlam,flmin,ulam,   zipper
     *thr,maxit,xv,  ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 zipper
19500 do 19501 k=1,lmu                                                     zipper
      nk=nin(k)                                                            zipper
19510 do 19511 j=1,nr                                                      zipper
19520 do 19521 l=1,nk                                                      zipper
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  zipper
19521 continue                                                             zipper
19522 continue                                                             zipper
      if(intr .ne. 0)goto 19541                                            zipper
      a0(j,k)=0.0                                                          zipper
      goto 19551                                                           zipper
19541 continue                                                             zipper
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 zipper
19551 continue                                                             zipper
19531 continue                                                             zipper
19511 continue                                                             zipper
19512 continue                                                             zipper
19501 continue                                                             zipper
19502 continue                                                             zipper
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multstandard1  (no,ni,nr,x,y,w,isd,jsd,intr,ju,xm,xs,ym   zipper
     *,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(n   zipper
     *r),ys(nr)
      integer ju(ni)                                                       zipper
      double precision, dimension (:), allocatable :: v                         
      allocate(v(1:no),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      w=w/sum(w)                                                           zipper
      v=sqrt(w)                                                            zipper
      if(intr .ne. 0)goto 19571                                            zipper
19580 do 19581 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 19581                                             zipper
      xm(j)=0.0                                                            zipper
      x(:,j)=v*x(:,j)                                                      zipper
      z=dot_product(x(:,j),x(:,j))                                         zipper
      if(isd .le. 0)goto 19601                                             zipper
      xbq=dot_product(v,x(:,j))**2                                         zipper
      vc=z-xbq                                                             zipper
      xs(j)=sqrt(vc)                                                       zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
      xv(j)=1.0+xbq/vc                                                     zipper
      goto 19611                                                           zipper
19601 continue                                                             zipper
      xs(j)=1.0                                                            zipper
      xv(j)=z                                                              zipper
19611 continue                                                             zipper
19591 continue                                                             zipper
19581 continue                                                             zipper
19582 continue                                                             zipper
      ys0=0.0                                                              zipper
19620 do 19621 j=1,nr                                                      zipper
      ym(j)=0.0                                                            zipper
      y(:,j)=v*y(:,j)                                                      zipper
      z=dot_product(y(:,j),y(:,j))                                         zipper
      if(jsd .le. 0)goto 19641                                             zipper
      u=z-dot_product(v,y(:,j))**2                                         zipper
      ys0=ys0+z/u                                                          zipper
      ys(j)=sqrt(u)                                                        zipper
      y(:,j)=y(:,j)/ys(j)                                                  zipper
      goto 19651                                                           zipper
19641 continue                                                             zipper
      ys(j)=1.0                                                            zipper
      ys0=ys0+z                                                            zipper
19651 continue                                                             zipper
19631 continue                                                             zipper
19621 continue                                                             zipper
19622 continue                                                             zipper
      go to 10700                                                          zipper
19571 continue                                                             zipper
19660 do 19661 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 19661                                             zipper
      xm(j)=dot_product(w,x(:,j))                                          zipper
      x(:,j)=v*(x(:,j)-xm(j))                                              zipper
      xv(j)=dot_product(x(:,j),x(:,j))                                     zipper
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       zipper
19661 continue                                                             zipper
19662 continue                                                             zipper
      if(isd .ne. 0)goto 19681                                             zipper
      xs=1.0                                                               zipper
      goto 19691                                                           zipper
19681 continue                                                             zipper
19700 do 19701 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 19701                                             zipper
      x(:,j)=x(:,j)/xs(j)                                                  zipper
19701 continue                                                             zipper
19702 continue                                                             zipper
      xv=1.0                                                               zipper
19691 continue                                                             zipper
19671 continue                                                             zipper
      ys0=0.0                                                              zipper
19710 do 19711 j=1,nr                                                      zipper
      ym(j)=dot_product(w,y(:,j))                                          zipper
      y(:,j)=v*(y(:,j)-ym(j))                                              zipper
      z=dot_product(y(:,j),y(:,j))                                         zipper
      if(jsd .le. 0)goto 19731                                             zipper
      ys(j)=sqrt(z)                                                        zipper
      y(:,j)=y(:,j)/ys(j)                                                  zipper
      goto 19741                                                           zipper
19731 continue                                                             zipper
      ys0=ys0+z                                                            zipper
19741 continue                                                             zipper
19721 continue                                                             zipper
19711 continue                                                             zipper
19712 continue                                                             zipper
      if(jsd .ne. 0)goto 19761                                             zipper
      ys=1.0                                                               zipper
      goto 19771                                                           zipper
19761 continue                                                             zipper
      ys0=nr                                                               zipper
19771 continue                                                             zipper
19751 continue                                                             zipper
10700 continue                                                             zipper
      deallocate(v)                                                        zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multelnet2(beta,ni,nr,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,   zipper
     *ulam,thri,  maxit,xv,ys0,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision vp(ni),y(no,nr),x(no,ni),ulam(nlam),ao(nx,nr,nlam   zipper
     *)
      double precision rsqo(nlam),almo(nlam),xv(ni),cl(2,nr,ni)            zipper
      integer ju(ni),ia(nx),kin(nlam)                                      zipper
      double precision, dimension (:), allocatable :: g,gk,del,gj               
      integer, dimension (:), allocatable :: mm,ix,isc                          
      double precision, dimension (:,:), allocatable :: a                       
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               zipper
      allocate(gj(1:nr),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(gk(1:nr),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(del(1:nr),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(g(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ix(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(isc(1:nr),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      bta=beta                                                             zipper
      omb=1.0-bta                                                          zipper
      ix=0                                                                 zipper
      thr=thri*ys0/nr                                                      zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 19791                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
19791 continue                                                             zipper
      rsq=ys0                                                              zipper
      a=0.0                                                                zipper
      mm=0                                                                 zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      iz=0                                                                 zipper
      mnl=min(mnlam,nlam)                                                  zipper
      alm=0.0                                                              zipper
19800 do 19801 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 19801                                             zipper
      g(j)=0.0                                                             zipper
19810 do 19811 k=1,nr                                                      zipper
      g(j)=g(j)+dot_product(y(:,k),x(:,j))**2                              zipper
19811 continue                                                             zipper
19812 continue                                                             zipper
      g(j)=sqrt(g(j))                                                      zipper
19801 continue                                                             zipper
19802 continue                                                             zipper
19820 do 19821 m=1,nlam                                                    zipper
      alm0=alm                                                             zipper
      if(flmin .lt. 1.0)goto 19841                                         zipper
      alm=ulam(m)                                                          zipper
      goto 19831                                                           zipper
19841 if(m .le. 2)goto 19851                                               zipper
      alm=alm*alf                                                          zipper
      goto 19831                                                           zipper
19851 if(m .ne. 1)goto 19861                                               zipper
      alm=big                                                              zipper
      goto 19871                                                           zipper
19861 continue                                                             zipper
      alm0=0.0                                                             zipper
19880 do 19881 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 19881                                             zipper
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           zipper
19881 continue                                                             zipper
19882 continue                                                             zipper
      alm0=alm0/max(bta,1.0d-3)                                            zipper
      alm=alf*alm0                                                         zipper
19871 continue                                                             zipper
19831 continue                                                             zipper
      dem=alm*omb                                                          zipper
      ab=alm*bta                                                           zipper
      rsq0=rsq                                                             zipper
      jz=1                                                                 zipper
      tlam=bta*(2.0*alm-alm0)                                              zipper
19890 do 19891 k=1,ni                                                      zipper
      if(ix(k).eq.1)goto 19891                                             zipper
      if(ju(k).eq.0)goto 19891                                             zipper
      if(g(k).gt.tlam*vp(k)) ix(k)=1                                       zipper
19891 continue                                                             zipper
19892 continue                                                             zipper
19900 continue                                                             zipper
19901 continue                                                             zipper
      if(iz*jz.ne.0) go to 10360                                           zipper
10880 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
19910 do 19911 k=1,ni                                                      zipper
      if(ix(k).eq.0)goto 19911                                             zipper
      gkn=0.0                                                              zipper
19920 do 19921 j=1,nr                                                      zipper
      gj(j)=dot_product(y(:,j),x(:,k))                                     zipper
      gk(j)=gj(j)+a(j,k)*xv(k)                                             zipper
      gkn=gkn+gk(j)**2                                                     zipper
19921 continue                                                             zipper
19922 continue                                                             zipper
      gkn=sqrt(gkn)                                                        zipper
      u=1.0-ab*vp(k)/gkn                                                   zipper
      del=a(:,k)                                                           zipper
      if(u .gt. 0.0)goto 19941                                             zipper
      a(:,k)=0.0                                                           zipper
      goto 19951                                                           zipper
19941 continue                                                             zipper
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      zipper
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   zipper
     *,isc,jerr)
      if(jerr.ne.0) return                                                 zipper
19951 continue                                                             zipper
19931 continue                                                             zipper
      del=a(:,k)-del                                                       zipper
      if(maxval(abs(del)).le.0.0)goto 19911                                zipper
19960 do 19961 j=1,nr                                                      zipper
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              zipper
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          zipper
      dlx=max(dlx,xv(k)*del(j)**2)                                         zipper
19961 continue                                                             zipper
19962 continue                                                             zipper
      if(mm(k) .ne. 0)goto 19981                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 19912                                              zipper
      mm(k)=nin                                                            zipper
      ia(nin)=k                                                            zipper
19981 continue                                                             zipper
19911 continue                                                             zipper
19912 continue                                                             zipper
      if(nin.gt.nx)goto 19902                                              zipper
      if(dlx .ge. thr)goto 20001                                           zipper
      ixx=0                                                                zipper
20010 do 20011 k=1,ni                                                      zipper
      if(ix(k).eq.1)goto 20011                                             zipper
      if(ju(k).eq.0)goto 20011                                             zipper
      g(k)=0.0                                                             zipper
20020 do 20021 j=1,nr                                                      zipper
      g(k)=g(k)+dot_product(y(:,j),x(:,k))**2                              zipper
20021 continue                                                             zipper
20022 continue                                                             zipper
      g(k)=sqrt(g(k))                                                      zipper
      if(g(k) .le. ab*vp(k))goto 20041                                     zipper
      ix(k)=1                                                              zipper
      ixx=1                                                                zipper
20041 continue                                                             zipper
20011 continue                                                             zipper
20012 continue                                                             zipper
      if(ixx.eq.1) go to 10880                                             zipper
      goto 19902                                                           zipper
20001 continue                                                             zipper
      if(nlp .le. maxit)goto 20061                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
20061 continue                                                             zipper
10360 continue                                                             zipper
      iz=1                                                                 zipper
20070 continue                                                             zipper
20071 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
20080 do 20081 l=1,nin                                                     zipper
      k=ia(l)                                                              zipper
      gkn=0.0                                                              zipper
20090 do 20091 j=1,nr                                                      zipper
      gj(j)=dot_product(y(:,j),x(:,k))                                     zipper
      gk(j)=gj(j)+a(j,k)*xv(k)                                             zipper
      gkn=gkn+gk(j)**2                                                     zipper
20091 continue                                                             zipper
20092 continue                                                             zipper
      gkn=sqrt(gkn)                                                        zipper
      u=1.0-ab*vp(k)/gkn                                                   zipper
      del=a(:,k)                                                           zipper
      if(u .gt. 0.0)goto 20111                                             zipper
      a(:,k)=0.0                                                           zipper
      goto 20121                                                           zipper
20111 continue                                                             zipper
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      zipper
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   zipper
     *,isc,jerr)
      if(jerr.ne.0) return                                                 zipper
20121 continue                                                             zipper
20101 continue                                                             zipper
      del=a(:,k)-del                                                       zipper
      if(maxval(abs(del)).le.0.0)goto 20081                                zipper
20130 do 20131 j=1,nr                                                      zipper
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              zipper
      y(:,j)=y(:,j)-del(j)*x(:,k)                                          zipper
      dlx=max(dlx,xv(k)*del(j)**2)                                         zipper
20131 continue                                                             zipper
20132 continue                                                             zipper
20081 continue                                                             zipper
20082 continue                                                             zipper
      if(dlx.lt.thr)goto 20072                                             zipper
      if(nlp .le. maxit)goto 20151                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
20151 continue                                                             zipper
      goto 20071                                                           zipper
20072 continue                                                             zipper
      jz=0                                                                 zipper
      goto 19901                                                           zipper
19902 continue                                                             zipper
      if(nin .le. nx)goto 20171                                            zipper
      jerr=-10000-m                                                        zipper
      goto 19822                                                           zipper
20171 continue                                                             zipper
      if(nin .le. 0)goto 20191                                             zipper
20200 do 20201 j=1,nr                                                      zipper
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         zipper
20201 continue                                                             zipper
20202 continue                                                             zipper
20191 continue                                                             zipper
      kin(m)=nin                                                           zipper
      rsqo(m)=1.0-rsq/ys0                                                  zipper
      almo(m)=alm                                                          zipper
      lmu=m                                                                zipper
      if(m.lt.mnl)goto 19821                                               zipper
      if(flmin.ge.1.0)goto 19821                                           zipper
      me=0                                                                 zipper
20210 do 20211 j=1,nin                                                     zipper
      if(ao(j,1,m).ne.0.0) me=me+1                                         zipper
20211 continue                                                             zipper
20212 continue                                                             zipper
      if(me.gt.ne)goto 19822                                               zipper
      if(rsq0-rsq.lt.sml*rsq)goto 19822                                    zipper
      if(rsqo(m).gt.rsqmax)goto 19822                                      zipper
19821 continue                                                             zipper
19822 continue                                                             zipper
      deallocate(a,mm,g,ix,del,gj,gk)                                      zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine chkbnds(nr,gk,gkn,xv,cl,al1,al2,a,isc,jerr)               zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision gk(nr),cl(2,nr),a(nr)                               zipper
      integer isc(nr)                                                      zipper
      kerr=0                                                               zipper
      al1p=1.0+al1/xv                                                      zipper
      al2p=al2/xv                                                          zipper
      isc=0                                                                zipper
      gsq=gkn**2                                                           zipper
      asq=dot_product(a,a)                                                 zipper
      usq=0.0                                                              zipper
      u=0.0                                                                zipper
      kn=-1                                                                zipper
20220 continue                                                             zipper
20221 continue                                                             zipper
      vmx=0.0                                                              zipper
20230 do 20231 k=1,nr                                                      zipper
      v=max(a(k)-cl(2,k),cl(1,k)-a(k))                                     zipper
      if(v .le. vmx)goto 20251                                             zipper
      vmx=v                                                                zipper
      kn=k                                                                 zipper
20251 continue                                                             zipper
20231 continue                                                             zipper
20232 continue                                                             zipper
      if(vmx.le.0.0)goto 20222                                             zipper
      if(isc(kn).ne.0)goto 20222                                           zipper
      gsq=gsq-gk(kn)**2                                                    zipper
      g=sqrt(gsq)/xv                                                       zipper
      if(a(kn).lt.cl(1,kn)) u=cl(1,kn)                                     zipper
      if(a(kn).gt.cl(2,kn)) u=cl(2,kn)                                     zipper
      usq=usq+u**2                                                         zipper
      if(usq .ne. 0.0)goto 20271                                           zipper
      b=max(0d0,(g-al2p)/al1p)                                             zipper
      goto 20281                                                           zipper
20271 continue                                                             zipper
      b0=sqrt(asq-a(kn)**2)                                                zipper
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     zipper
      if(kerr.ne.0)goto 20222                                              zipper
20281 continue                                                             zipper
20261 continue                                                             zipper
      asq=usq+b**2                                                         zipper
      if(asq .gt. 0.0)goto 20301                                           zipper
      a=0.0                                                                zipper
      goto 20222                                                           zipper
20301 continue                                                             zipper
      a(kn)=u                                                              zipper
      isc(kn)=1                                                            zipper
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     zipper
20310 do 20311 j=1,nr                                                      zipper
      if(isc(j).eq.0) a(j)=f*gk(j)                                         zipper
20311 continue                                                             zipper
20312 continue                                                             zipper
      goto 20221                                                           zipper
20222 continue                                                             zipper
      if(kerr.ne.0) jerr=kerr                                              zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine chkbnds1(nr,gk,gkn,xv,cl1,cl2,al1,al2,a,isc,jerr)         zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision gk(nr),a(nr)                                        zipper
      integer isc(nr)                                                      zipper
      kerr=0                                                               zipper
      al1p=1.0+al1/xv                                                      zipper
      al2p=al2/xv                                                          zipper
      isc=0                                                                zipper
      gsq=gkn**2                                                           zipper
      asq=dot_product(a,a)                                                 zipper
      usq=0.0                                                              zipper
      u=0.0                                                                zipper
      kn=-1                                                                zipper
20320 continue                                                             zipper
20321 continue                                                             zipper
      vmx=0.0                                                              zipper
20330 do 20331 k=1,nr                                                      zipper
      v=max(a(k)-cl2,cl1-a(k))                                             zipper
      if(v .le. vmx)goto 20351                                             zipper
      vmx=v                                                                zipper
      kn=k                                                                 zipper
20351 continue                                                             zipper
20331 continue                                                             zipper
20332 continue                                                             zipper
      if(vmx.le.0.0)goto 20322                                             zipper
      if(isc(kn).ne.0)goto 20322                                           zipper
      gsq=gsq-gk(kn)**2                                                    zipper
      g=sqrt(gsq)/xv                                                       zipper
      if(a(kn).lt.cl1) u=cl1                                               zipper
      if(a(kn).gt.cl2) u=cl2                                               zipper
      usq=usq+u**2                                                         zipper
      if(usq .ne. 0.0)goto 20371                                           zipper
      b=max(0d0,(g-al2p)/al1p)                                             zipper
      goto 20381                                                           zipper
20371 continue                                                             zipper
      b0=sqrt(asq-a(kn)**2)                                                zipper
      b=bnorm(b0,al1p,al2p,g,usq,kerr)                                     zipper
      if(kerr.ne.0)goto 20322                                              zipper
20381 continue                                                             zipper
20361 continue                                                             zipper
      asq=usq+b**2                                                         zipper
      if(asq .gt. 0.0)goto 20401                                           zipper
      a=0.0                                                                zipper
      goto 20322                                                           zipper
20401 continue                                                             zipper
      a(kn)=u                                                              zipper
      isc(kn)=1                                                            zipper
      f=1.0/(xv*(al1p+al2p/sqrt(asq)))                                     zipper
20410 do 20411 j=1,nr                                                      zipper
      if(isc(j).eq.0) a(j)=f*gk(j)                                         zipper
20411 continue                                                             zipper
20412 continue                                                             zipper
      goto 20321                                                           zipper
20322 continue                                                             zipper
      if(kerr.ne.0) jerr=kerr                                              zipper
      return                                                               zipper
      end                                                                  zipper
      function bnorm(b0,al1p,al2p,g,usq,jerr)                              zipper
      implicit double precision(a-h,o-z)                                   zipper
      data thr,mxit /1.0d-10,100/                                          zipper
      b=b0                                                                 zipper
      zsq=b**2+usq                                                         zipper
      if(zsq .gt. 0.0)goto 20431                                           zipper
      bnorm=0.0                                                            zipper
      return                                                               zipper
20431 continue                                                             zipper
      z=sqrt(zsq)                                                          zipper
      f=b*(al1p+al2p/z)-g                                                  zipper
      jerr=0                                                               zipper
20440 do 20441 it=1,mxit                                                   zipper
      b=b-f/(al1p+al2p*usq/(z*zsq))                                        zipper
      zsq=b**2+usq                                                         zipper
      if(zsq .gt. 0.0)goto 20461                                           zipper
      bnorm=0.0                                                            zipper
      return                                                               zipper
20461 continue                                                             zipper
      z=sqrt(zsq)                                                          zipper
      f=b*(al1p+al2p/z)-g                                                  zipper
      if(abs(f).le.thr)goto 20442                                          zipper
      if(b .gt. 0.0)goto 20481                                             zipper
      b=0.0                                                                zipper
      goto 20442                                                           zipper
20481 continue                                                             zipper
20441 continue                                                             zipper
20442 continue                                                             zipper
      bnorm=b                                                              zipper
      if(it.ge.mxit) jerr=90000                                            zipper
      return                                                               zipper
      entry chg_bnorm(arg,irg)                                             zipper
      bnorm = 0.0                                                          zipper
      thr=arg                                                              zipper
      mxit=irg                                                             zipper
      return                                                               zipper
      entry get_bnorm(arg,irg)                                             zipper
      bnorm = 0.0                                                          zipper
      arg=thr                                                              zipper
      irg=mxit                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multsolns(ni,nx,nr,lmu,a,ia,nin,b)                        zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision a(nx,nr,lmu),b(ni,nr,lmu)                           zipper
      integer ia(nx),nin(lmu)                                              zipper
20490 do 20491 lam=1,lmu                                                   zipper
      call multuncomp(ni,nr,nx,a(1,1,lam),ia,nin(lam),b(1,1,lam))          zipper
20491 continue                                                             zipper
20492 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multuncomp(ni,nr,nx,ca,ia,nin,a)                          zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision ca(nx,nr),a(ni,nr)                                  zipper
      integer ia(nx)                                                       zipper
      a=0.0                                                                zipper
      if(nin .le. 0)goto 20511                                             zipper
20520 do 20521 j=1,nr                                                      zipper
      a(ia(1:nin),j)=ca(1:nin,j)                                           zipper
20521 continue                                                             zipper
20522 continue                                                             zipper
20511 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multmodval(nx,nr,a0,ca,ia,nin,n,x,f)                      zipper
      implicit double precision(a-h,o-z)                                   zipper
      double precision a0(nr),ca(nx,nr),x(n,*),f(nr,n)                     zipper
      integer ia(nx)                                                       zipper
20530 do 20531 i=1,n                                                       zipper
      f(:,i)=a0                                                            zipper
20531 continue                                                             zipper
20532 continue                                                             zipper
      if(nin.le.0) return                                                  zipper
20540 do 20541 i=1,n                                                       zipper
20550 do 20551 j=1,nr                                                      zipper
      f(j,i)=f(j,i)+dot_product(ca(1:nin,j),x(i,ia(1:nin)))                zipper
20551 continue                                                             zipper
20552 continue                                                             zipper
20541 continue                                                             zipper
20542 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multspelnet  (parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,   zipper
     *nlam,flmin,ulam,thr,isd,  jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,
     *nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no,nr),w(no),vp(ni),ulam(nlam),cl(2,ni)      zipper
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      zipper
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           zipper
      double precision, dimension (:), allocatable :: vq;                       
      if(maxval(vp) .gt. 0.0)goto 20571                                    zipper
      jerr=10000                                                           zipper
      return                                                               zipper
20571 continue                                                             zipper
      allocate(vq(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      vq=max(0d0,vp)                                                       zipper
      vq=vq*ni/sum(vq)                                                     zipper
      call multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vq,cl,ne,nx,nlam,fl   zipper
     *min,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jer
     *r)
      deallocate(vq)                                                       zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multspelnetn(parm,no,ni,nr,x,ix,jx,y,w,jd,vp,cl,ne,nx,n   zipper
     *lam,flmin,  ulam,thr,isd,jsd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),vp(ni),y(no,nr),w(no),ulam(nlam),cl(2,ni)      zipper
      double precision ca(nx,nr,nlam),a0(nr,nlam),rsq(nlam),alm(nlam)      zipper
      integer ix(*),jx(*),jd(*),ia(nx),nin(nlam)                           zipper
      double precision, dimension (:), allocatable :: xm,xs,xv,ym,ys            
      integer, dimension (:), allocatable :: ju                                 
      double precision, dimension (:,:,:), allocatable :: clt                   
      allocate(clt(1:2,1:nr,1:ni),stat=jerr)                                    
      if(jerr.ne.0) return                                                      
      allocate(xm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xs(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ym(1:nr),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ys(1:nr),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ju(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(xv(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      call spchkvars(no,ni,x,ix,ju)                                        zipper
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                 zipper
      if(maxval(ju) .gt. 0)goto 20591                                      zipper
      jerr=7777                                                            zipper
      return                                                               zipper
20591 continue                                                             zipper
      call multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,  xm,xs,   zipper
     *ym,ys,xv,ys0,jerr)
      if(jerr.ne.0) return                                                 zipper
20600 do 20601 j=1,ni                                                      zipper
20610 do 20611 k=1,nr                                                      zipper
20620 do 20621 i=1,2                                                       zipper
      clt(i,k,j)=cl(i,j)                                                   zipper
20621 continue                                                             zipper
20622 continue                                                             zipper
20611 continue                                                             zipper
20612 continue                                                             zipper
20601 continue                                                             zipper
20602 continue                                                             zipper
      if(isd .le. 0)goto 20641                                             zipper
20650 do 20651 j=1,ni                                                      zipper
20660 do 20661 k=1,nr                                                      zipper
20670 do 20671 i=1,2                                                       zipper
      clt(i,k,j)=clt(i,k,j)*xs(j)                                          zipper
20671 continue                                                             zipper
20672 continue                                                             zipper
20661 continue                                                             zipper
20662 continue                                                             zipper
20651 continue                                                             zipper
20652 continue                                                             zipper
20641 continue                                                             zipper
      if(jsd .le. 0)goto 20691                                             zipper
20700 do 20701 j=1,ni                                                      zipper
20710 do 20711 k=1,nr                                                      zipper
20720 do 20721 i=1,2                                                       zipper
      clt(i,k,j)=clt(i,k,j)/ys(k)                                          zipper
20721 continue                                                             zipper
20722 continue                                                             zipper
20711 continue                                                             zipper
20712 continue                                                             zipper
20701 continue                                                             zipper
20702 continue                                                             zipper
20691 continue                                                             zipper
      call multspelnet2(parm,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,clt,nlam,f   zipper
     *lmin,  ulam,thr,maxit,xm,xs,xv,ys0,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                 zipper
20730 do 20731 k=1,lmu                                                     zipper
      nk=nin(k)                                                            zipper
20740 do 20741 j=1,nr                                                      zipper
20750 do 20751 l=1,nk                                                      zipper
      ca(l,j,k)=ys(j)*ca(l,j,k)/xs(ia(l))                                  zipper
20751 continue                                                             zipper
20752 continue                                                             zipper
      if(intr .ne. 0)goto 20771                                            zipper
      a0(j,k)=0.0                                                          zipper
      goto 20781                                                           zipper
20771 continue                                                             zipper
      a0(j,k)=ym(j)-dot_product(ca(1:nk,j,k),xm(ia(1:nk)))                 zipper
20781 continue                                                             zipper
20761 continue                                                             zipper
20741 continue                                                             zipper
20742 continue                                                             zipper
20731 continue                                                             zipper
20732 continue                                                             zipper
      deallocate(xm,xs,ym,ys,ju,xv,clt)                                    zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multspstandard1(no,ni,nr,x,ix,jx,y,w,ju,isd,jsd,intr,     zipper
     *xm,xs,ym,ys,xv,ys0,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no,nr),w(no),xm(ni),xs(ni),xv(ni),ym(nr),y   zipper
     *s(nr)
      integer ix(*),jx(*),ju(ni)                                           zipper
      w=w/sum(w)                                                           zipper
      if(intr .ne. 0)goto 20801                                            zipper
20810 do 20811 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 20811                                             zipper
      xm(j)=0.0                                                            zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      z=dot_product(w(jx(jb:je)),x(jb:je)**2)                              zipper
      if(isd .le. 0)goto 20831                                             zipper
      xbq=dot_product(w(jx(jb:je)),x(jb:je))**2                            zipper
      vc=z-xbq                                                             zipper
      xs(j)=sqrt(vc)                                                       zipper
      xv(j)=1.0+xbq/vc                                                     zipper
      goto 20841                                                           zipper
20831 continue                                                             zipper
      xs(j)=1.0                                                            zipper
      xv(j)=z                                                              zipper
20841 continue                                                             zipper
20821 continue                                                             zipper
20811 continue                                                             zipper
20812 continue                                                             zipper
      ys0=0.0                                                              zipper
20850 do 20851 j=1,nr                                                      zipper
      ym(j)=0.0                                                            zipper
      z=dot_product(w,y(:,j)**2)                                           zipper
      if(jsd .le. 0)goto 20871                                             zipper
      u=z-dot_product(w,y(:,j))**2                                         zipper
      ys0=ys0+z/u                                                          zipper
      ys(j)=sqrt(u)                                                        zipper
      y(:,j)=y(:,j)/ys(j)                                                  zipper
      goto 20881                                                           zipper
20871 continue                                                             zipper
      ys(j)=1.0                                                            zipper
      ys0=ys0+z                                                            zipper
20881 continue                                                             zipper
20861 continue                                                             zipper
20851 continue                                                             zipper
20852 continue                                                             zipper
      return                                                               zipper
20801 continue                                                             zipper
20890 do 20891 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 20891                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      xm(j)=dot_product(w(jx(jb:je)),x(jb:je))                             zipper
      xv(j)=dot_product(w(jx(jb:je)),x(jb:je)**2)-xm(j)**2                 zipper
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                       zipper
20891 continue                                                             zipper
20892 continue                                                             zipper
      if(isd .ne. 0)goto 20911                                             zipper
      xs=1.0                                                               zipper
      goto 20921                                                           zipper
20911 continue                                                             zipper
      xv=1.0                                                               zipper
20921 continue                                                             zipper
20901 continue                                                             zipper
      ys0=0.0                                                              zipper
20930 do 20931 j=1,nr                                                      zipper
      ym(j)=dot_product(w,y(:,j))                                          zipper
      y(:,j)=y(:,j)-ym(j)                                                  zipper
      z=dot_product(w,y(:,j)**2)                                           zipper
      if(jsd .le. 0)goto 20951                                             zipper
      ys(j)=sqrt(z)                                                        zipper
      y(:,j)=y(:,j)/ys(j)                                                  zipper
      goto 20961                                                           zipper
20951 continue                                                             zipper
      ys0=ys0+z                                                            zipper
20961 continue                                                             zipper
20941 continue                                                             zipper
20931 continue                                                             zipper
20932 continue                                                             zipper
      if(jsd .ne. 0)goto 20981                                             zipper
      ys=1.0                                                               zipper
      goto 20991                                                           zipper
20981 continue                                                             zipper
      ys0=nr                                                               zipper
20991 continue                                                             zipper
20971 continue                                                             zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multspelnet2(beta,ni,nr,y,w,no,ne,nx,x,ix,jx,ju,vp,cl,n   zipper
     *lam,flmin,  ulam,thri,maxit,xm,xs,xv,ys0,lmu,ao,ia,kin,rsqo,almo,n
     *lp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision y(no,nr),w(no),x(*),vp(ni),ulam(nlam),cl(2,nr,ni)   zipper
      double precision ao(nx,nr,nlam),rsqo(nlam),almo(nlam),xm(ni),xs(ni   zipper
     *),xv(ni)
      integer ix(*),jx(*),ju(ni),ia(nx),kin(nlam)                          zipper
      double precision, dimension (:), allocatable :: g,gj,gk,del,o             
      integer, dimension (:), allocatable :: mm,iy,isc                          
      double precision, dimension (:,:), allocatable :: a                       
      allocate(a(1:nr,1:ni),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,rsqmax,pmin,exmx)               zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(g(1:ni),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(gj(1:nr),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(gk(1:nr),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(del(1:nr),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(o(1:nr),stat=jerr)                                          zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(iy(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(isc(1:nr),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      bta=beta                                                             zipper
      omb=1.0-bta                                                          zipper
      alm=0.0                                                              zipper
      iy=0                                                                 zipper
      thr=thri*ys0/nr                                                      zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 21011                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
21011 continue                                                             zipper
      rsq=ys0                                                              zipper
      a=0.0                                                                zipper
      mm=0                                                                 zipper
      o=0.0                                                                zipper
      nlp=0                                                                zipper
      nin=nlp                                                              zipper
      iz=0                                                                 zipper
      mnl=min(mnlam,nlam)                                                  zipper
21020 do 21021 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 21021                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      g(j)=0.0                                                             zipper
21030 do 21031 k=1,nr                                                      zipper
      g(j)=g(j)+(dot_product(y(jx(jb:je),k),w(jx(jb:je))*x(jb:je))/xs(j)   zipper
     *)**zipper
21031 continue                                                             zipper
21032 continue                                                             zipper
      g(j)=sqrt(g(j))                                                      zipper
21021 continue                                                             zipper
21022 continue                                                             zipper
21040 do 21041 m=1,nlam                                                    zipper
      alm0=alm                                                             zipper
      if(flmin .lt. 1.0)goto 21061                                         zipper
      alm=ulam(m)                                                          zipper
      goto 21051                                                           zipper
21061 if(m .le. 2)goto 21071                                               zipper
      alm=alm*alf                                                          zipper
      goto 21051                                                           zipper
21071 if(m .ne. 1)goto 21081                                               zipper
      alm=big                                                              zipper
      goto 21091                                                           zipper
21081 continue                                                             zipper
      alm0=0.0                                                             zipper
21100 do 21101 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 21101                                             zipper
      if(vp(j).gt.0.0) alm0=max(alm0,g(j)/vp(j))                           zipper
21101 continue                                                             zipper
21102 continue                                                             zipper
      alm0=alm0/max(bta,1.0d-3)                                            zipper
      alm=alf*alm0                                                         zipper
21091 continue                                                             zipper
21051 continue                                                             zipper
      dem=alm*omb                                                          zipper
      ab=alm*bta                                                           zipper
      rsq0=rsq                                                             zipper
      jz=1                                                                 zipper
      tlam=bta*(2.0*alm-alm0)                                              zipper
21110 do 21111 k=1,ni                                                      zipper
      if(iy(k).eq.1)goto 21111                                             zipper
      if(ju(k).eq.0)goto 21111                                             zipper
      if(g(k).gt.tlam*vp(k)) iy(k)=1                                       zipper
21111 continue                                                             zipper
21112 continue                                                             zipper
21120 continue                                                             zipper
21121 continue                                                             zipper
      if(iz*jz.ne.0) go to 10360                                           zipper
10880 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
21130 do 21131 k=1,ni                                                      zipper
      if(iy(k).eq.0)goto 21131                                             zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      gkn=0.0                                                              zipper
21140 do 21141 j=1,nr                                                      zipper
      gj(j)=dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(k)   zipper
      gk(j)=gj(j)+a(j,k)*xv(k)                                             zipper
      gkn=gkn+gk(j)**2                                                     zipper
21141 continue                                                             zipper
21142 continue                                                             zipper
      gkn=sqrt(gkn)                                                        zipper
      u=1.0-ab*vp(k)/gkn                                                   zipper
      del=a(:,k)                                                           zipper
      if(u .gt. 0.0)goto 21161                                             zipper
      a(:,k)=0.0                                                           zipper
      goto 21171                                                           zipper
21161 continue                                                             zipper
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      zipper
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   zipper
     *,isc,jerr)
      if(jerr.ne.0) return                                                 zipper
21171 continue                                                             zipper
21151 continue                                                             zipper
      del=a(:,k)-del                                                       zipper
      if(maxval(abs(del)).le.0.0)goto 21131                                zipper
      if(mm(k) .ne. 0)goto 21191                                           zipper
      nin=nin+1                                                            zipper
      if(nin.gt.nx)goto 21132                                              zipper
      mm(k)=nin                                                            zipper
      ia(nin)=k                                                            zipper
21191 continue                                                             zipper
21200 do 21201 j=1,nr                                                      zipper
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              zipper
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  zipper
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         zipper
      dlx=max(xv(k)*del(j)**2,dlx)                                         zipper
21201 continue                                                             zipper
21202 continue                                                             zipper
21131 continue                                                             zipper
21132 continue                                                             zipper
      if(nin.gt.nx)goto 21122                                              zipper
      if(dlx .ge. thr)goto 21221                                           zipper
      ixx=0                                                                zipper
21230 do 21231 j=1,ni                                                      zipper
      if(iy(j).eq.1)goto 21231                                             zipper
      if(ju(j).eq.0)goto 21231                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      g(j)=0.0                                                             zipper
21240 do 21241 k=1,nr                                                      zipper
      g(j)=g(j)+  (dot_product(y(jx(jb:je),k)+o(k),w(jx(jb:je))*x(jb:je)   zipper
     *)/xs(j))**zipper
21241 continue                                                             zipper
21242 continue                                                             zipper
      g(j)=sqrt(g(j))                                                      zipper
      if(g(j) .le. ab*vp(j))goto 21261                                     zipper
      iy(j)=1                                                              zipper
      ixx=1                                                                zipper
21261 continue                                                             zipper
21231 continue                                                             zipper
21232 continue                                                             zipper
      if(ixx.eq.1) go to 10880                                             zipper
      goto 21122                                                           zipper
21221 continue                                                             zipper
      if(nlp .le. maxit)goto 21281                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
21281 continue                                                             zipper
10360 continue                                                             zipper
      iz=1                                                                 zipper
21290 continue                                                             zipper
21291 continue                                                             zipper
      nlp=nlp+1                                                            zipper
      dlx=0.0                                                              zipper
21300 do 21301 l=1,nin                                                     zipper
      k=ia(l)                                                              zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      gkn=0.0                                                              zipper
21310 do 21311 j=1,nr                                                      zipper
      gj(j)=  dot_product(y(jx(jb:je),j)+o(j),w(jx(jb:je))*x(jb:je))/xs(   zipper
     *k)
      gk(j)=gj(j)+a(j,k)*xv(k)                                             zipper
      gkn=gkn+gk(j)**2                                                     zipper
21311 continue                                                             zipper
21312 continue                                                             zipper
      gkn=sqrt(gkn)                                                        zipper
      u=1.0-ab*vp(k)/gkn                                                   zipper
      del=a(:,k)                                                           zipper
      if(u .gt. 0.0)goto 21331                                             zipper
      a(:,k)=0.0                                                           zipper
      goto 21341                                                           zipper
21331 continue                                                             zipper
      a(:,k)=gk*(u/(xv(k)+dem*vp(k)))                                      zipper
      call chkbnds(nr,gk,gkn,xv(k),cl(1,1,k),  dem*vp(k),ab*vp(k),a(:,k)   zipper
     *,isc,jerr)
      if(jerr.ne.0) return                                                 zipper
21341 continue                                                             zipper
21321 continue                                                             zipper
      del=a(:,k)-del                                                       zipper
      if(maxval(abs(del)).le.0.0)goto 21301                                zipper
21350 do 21351 j=1,nr                                                      zipper
      rsq=rsq-del(j)*(2.0*gj(j)-del(j)*xv(k))                              zipper
      y(jx(jb:je),j)=y(jx(jb:je),j)-del(j)*x(jb:je)/xs(k)                  zipper
      o(j)=o(j)+del(j)*xm(k)/xs(k)                                         zipper
      dlx=max(xv(k)*del(j)**2,dlx)                                         zipper
21351 continue                                                             zipper
21352 continue                                                             zipper
21301 continue                                                             zipper
21302 continue                                                             zipper
      if(dlx.lt.thr)goto 21292                                             zipper
      if(nlp .le. maxit)goto 21371                                         zipper
      jerr=-m                                                              zipper
      return                                                               zipper
21371 continue                                                             zipper
      goto 21291                                                           zipper
21292 continue                                                             zipper
      jz=0                                                                 zipper
      goto 21121                                                           zipper
21122 continue                                                             zipper
      if(nin .le. nx)goto 21391                                            zipper
      jerr=-10000-m                                                        zipper
      goto 21042                                                           zipper
21391 continue                                                             zipper
      if(nin .le. 0)goto 21411                                             zipper
21420 do 21421 j=1,nr                                                      zipper
      ao(1:nin,j,m)=a(j,ia(1:nin))                                         zipper
21421 continue                                                             zipper
21422 continue                                                             zipper
21411 continue                                                             zipper
      kin(m)=nin                                                           zipper
      rsqo(m)=1.0-rsq/ys0                                                  zipper
      almo(m)=alm                                                          zipper
      lmu=m                                                                zipper
      if(m.lt.mnl)goto 21041                                               zipper
      if(flmin.ge.1.0)goto 21041                                           zipper
      me=0                                                                 zipper
21430 do 21431 j=1,nin                                                     zipper
      if(ao(j,1,m).ne.0.0) me=me+1                                         zipper
21431 continue                                                             zipper
21432 continue                                                             zipper
      if(me.gt.ne)goto 21042                                               zipper
      if(rsq0-rsq.lt.sml*rsq)goto 21042                                    zipper
      if(rsqo(m).gt.rsqmax)goto 21042                                      zipper
21041 continue                                                             zipper
21042 continue                                                             zipper
      deallocate(a,mm,g,iy,gj,gk,del,o)                                    zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multlognetn(parm,no,ni,nc,x,y,g,w,ju,vp,cl,ne,nx,nlam,f   zipper
     *lmin,ulam,  shri,intr,maxit,xv,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jer
     *r)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(no,ni),y(no,nc),g(no,nc),w(no),vp(ni),ulam(nlam   zipper
     *),cl(2,ni)
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),xv(   zipper
     *ni)
      integer ju(ni),m(nx),kin(nlam)                                       zipper
      double precision, dimension (:,:), allocatable :: q,r,b,bs                
      double precision, dimension (:), allocatable :: sxp,sxpl,ga,gk,del        
      integer, dimension (:), allocatable :: mm,is,ixx,isc                      
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return					                                                 
      allocate(r(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      exmn=-exmx                                                           zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(is(1:max(nc,ni)),stat=jerr)                                 zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sxp(1:no),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sxpl(1:no),stat=jerr)                                       zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ixx(1:ni),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(gk(1:nc),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(del(1:nc),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(isc(1:nc),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      pmax=1.0-pmin                                                        zipper
      emin=pmin/pmax                                                       zipper
      emax=1.0/emin                                                        zipper
      bta=parm                                                             zipper
      omb=1.0-bta                                                          zipper
      dev1=0.0                                                             zipper
      dev0=0.0                                                             zipper
21440 do 21441 ic=1,nc                                                     zipper
      q0=dot_product(w,y(:,ic))                                            zipper
      if(q0 .gt. pmin)goto 21461                                           zipper
      jerr =8000+ic                                                        zipper
      return                                                               zipper
21461 continue                                                             zipper
      if(q0 .lt. pmax)goto 21481                                           zipper
      jerr =9000+ic                                                        zipper
      return                                                               zipper
21481 continue                                                             zipper
      if(intr .ne. 0)goto 21501                                            zipper
      q0=1.0/nc                                                            zipper
      b(0,ic)=0.0                                                          zipper
      goto 21511                                                           zipper
21501 continue                                                             zipper
      b(0,ic)=log(q0)                                                      zipper
      dev1=dev1-q0*b(0,ic)                                                 zipper
21511 continue                                                             zipper
21491 continue                                                             zipper
      b(1:ni,ic)=0.0                                                       zipper
21441 continue                                                             zipper
21442 continue                                                             zipper
      if(intr.eq.0) dev1=log(float(nc))                                    zipper
      ixx=0                                                                zipper
      al=0.0                                                               zipper
      if(nonzero(no*nc,g) .ne. 0)goto 21531                                zipper
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         zipper
      sxp=0.0                                                              zipper
21540 do 21541 ic=1,nc                                                     zipper
      q(:,ic)=exp(b(0,ic))                                                 zipper
      sxp=sxp+q(:,ic)                                                      zipper
21541 continue                                                             zipper
21542 continue                                                             zipper
      goto 21551                                                           zipper
21531 continue                                                             zipper
21560 do 21561 i=1,no                                                      zipper
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         zipper
21561 continue                                                             zipper
21562 continue                                                             zipper
      sxp=0.0                                                              zipper
      if(intr .ne. 0)goto 21581                                            zipper
      b(0,:)=0.0                                                           zipper
      goto 21591                                                           zipper
21581 continue                                                             zipper
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 zipper
      if(jerr.ne.0) return                                                 zipper
21591 continue                                                             zipper
21571 continue                                                             zipper
      dev1=0.0                                                             zipper
21600 do 21601 ic=1,nc                                                     zipper
      q(:,ic)=b(0,ic)+g(:,ic)                                              zipper
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             zipper
      q(:,ic)=exp(q(:,ic))                                                 zipper
      sxp=sxp+q(:,ic)                                                      zipper
21601 continue                                                             zipper
21602 continue                                                             zipper
      sxpl=w*log(sxp)                                                      zipper
21610 do 21611 ic=1,nc                                                     zipper
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  zipper
21611 continue                                                             zipper
21612 continue                                                             zipper
21551 continue                                                             zipper
21521 continue                                                             zipper
21620 do 21621 ic=1,nc                                                     zipper
21630 do 21631 i=1,no                                                      zipper
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               zipper
21631 continue                                                             zipper
21632 continue                                                             zipper
21621 continue                                                             zipper
21622 continue                                                             zipper
      dev0=dev0+dev1                                                       zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 21651                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
21651 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nin=0                                                                zipper
      nlp=0                                                                zipper
      mnl=min(mnlam,nlam)                                                  zipper
      bs=0.0                                                               zipper
      shr=shri*dev0                                                        zipper
      ga=0.0                                                               zipper
21660 do 21661 ic=1,nc                                                     zipper
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      zipper
21670 do 21671 j=1,ni                                                      zipper
      if(ju(j).ne.0) ga(j)=ga(j)+dot_product(r(:,ic),x(:,j))**2            zipper
21671 continue                                                             zipper
21672 continue                                                             zipper
21661 continue                                                             zipper
21662 continue                                                             zipper
      ga=sqrt(ga)                                                          zipper
21680 do 21681 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 21701                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 21691                                                           zipper
21701 if(ilm .le. 2)goto 21711                                             zipper
      al=al*alf                                                            zipper
      goto 21691                                                           zipper
21711 if(ilm .ne. 1)goto 21721                                             zipper
      al=big                                                               zipper
      goto 21731                                                           zipper
21721 continue                                                             zipper
      al0=0.0                                                              zipper
21740 do 21741 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 21741                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
21741 continue                                                             zipper
21742 continue                                                             zipper
      al0=al0/max(bta,1.0d-3)                                              zipper
      al=alf*al0                                                           zipper
21731 continue                                                             zipper
21691 continue                                                             zipper
      al2=al*omb                                                           zipper
      al1=al*bta                                                           zipper
      tlam=bta*(2.0*al-al0)                                                zipper
21750 do 21751 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 21751                                            zipper
      if(ju(k).eq.0)goto 21751                                             zipper
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1                                     zipper
21751 continue                                                             zipper
21752 continue                                                             zipper
10880 continue                                                             zipper
21760 continue                                                             zipper
21761 continue                                                             zipper
      ix=0                                                                 zipper
      jx=ix                                                                zipper
      kx=jx                                                                zipper
      t=0.0                                                                zipper
21770 do 21771 ic=1,nc                                                     zipper
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       zipper
21771 continue                                                             zipper
21772 continue                                                             zipper
      if(t .ge. eps)goto 21791                                             zipper
      kx=1                                                                 zipper
      goto 21762                                                           zipper
21791 continue                                                             zipper
      t=2.0*t                                                              zipper
      alt=al1/t                                                            zipper
      al2t=al2/t                                                           zipper
21800 do 21801 ic=1,nc                                                     zipper
      bs(0,ic)=b(0,ic)                                                     zipper
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          zipper
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    zipper
      d=0.0                                                                zipper
      if(intr.ne.0) d=sum(r(:,ic))                                         zipper
      if(d .eq. 0.0)goto 21821                                             zipper
      b(0,ic)=b(0,ic)+d                                                    zipper
      r(:,ic)=r(:,ic)-d*w                                                  zipper
      dlx=max(dlx,d**2)                                                    zipper
21821 continue                                                             zipper
21801 continue                                                             zipper
21802 continue                                                             zipper
21830 continue                                                             zipper
21831 continue                                                             zipper
      nlp=nlp+nc                                                           zipper
      dlx=0.0                                                              zipper
21840 do 21841 k=1,ni                                                      zipper
      if(ixx(k).eq.0)goto 21841                                            zipper
      gkn=0.0                                                              zipper
21850 do 21851 ic=1,nc                                                     zipper
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     zipper
      gkn=gkn+gk(ic)**2                                                    zipper
21851 continue                                                             zipper
21852 continue                                                             zipper
      gkn=sqrt(gkn)                                                        zipper
      u=1.0-alt*vp(k)/gkn                                                  zipper
      del=b(k,:)                                                           zipper
      if(u .gt. 0.0)goto 21871                                             zipper
      b(k,:)=0.0                                                           zipper
      goto 21881                                                           zipper
21871 continue                                                             zipper
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     zipper
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   zipper
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 zipper
21881 continue                                                             zipper
21861 continue                                                             zipper
      del=b(k,:)-del                                                       zipper
      if(maxval(abs(del)).le.0.0)goto 21841                                zipper
21890 do 21891 ic=1,nc                                                     zipper
      dlx=max(dlx,xv(k)*del(ic)**2)                                        zipper
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     zipper
21891 continue                                                             zipper
21892 continue                                                             zipper
      if(mm(k) .ne. 0)goto 21911                                           zipper
      nin=nin+1                                                            zipper
      if(nin .le. nx)goto 21931                                            zipper
      jx=1                                                                 zipper
      goto 21842                                                           zipper
21931 continue                                                             zipper
      mm(k)=nin                                                            zipper
      m(nin)=k                                                             zipper
21911 continue                                                             zipper
21841 continue                                                             zipper
21842 continue                                                             zipper
      if(jx.gt.0)goto 21832                                                zipper
      if(dlx.lt.shr)goto 21832                                             zipper
      if(nlp .le. maxit)goto 21951                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
21951 continue                                                             zipper
21960 continue                                                             zipper
21961 continue                                                             zipper
      nlp=nlp+nc                                                           zipper
      dlx=0.0                                                              zipper
21970 do 21971 l=1,nin                                                     zipper
      k=m(l)                                                               zipper
      gkn=0.0                                                              zipper
21980 do 21981 ic=1,nc                                                     zipper
      gk(ic)=dot_product(r(:,ic),x(:,k))+b(k,ic)*xv(k)                     zipper
      gkn=gkn+gk(ic)**2                                                    zipper
21981 continue                                                             zipper
21982 continue                                                             zipper
      gkn=sqrt(gkn)                                                        zipper
      u=1.0-alt*vp(k)/gkn                                                  zipper
      del=b(k,:)                                                           zipper
      if(u .gt. 0.0)goto 22001                                             zipper
      b(k,:)=0.0                                                           zipper
      goto 22011                                                           zipper
22001 continue                                                             zipper
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     zipper
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),  cl(2,k),vp(k)*al2t,alt*vp(   zipper
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 zipper
22011 continue                                                             zipper
21991 continue                                                             zipper
      del=b(k,:)-del                                                       zipper
      if(maxval(abs(del)).le.0.0)goto 21971                                zipper
22020 do 22021 ic=1,nc                                                     zipper
      dlx=max(dlx,xv(k)*del(ic)**2)                                        zipper
      r(:,ic)=r(:,ic)-del(ic)*w*x(:,k)                                     zipper
22021 continue                                                             zipper
22022 continue                                                             zipper
21971 continue                                                             zipper
21972 continue                                                             zipper
      if(dlx.lt.shr)goto 21962                                             zipper
      if(nlp .le. maxit)goto 22041                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
22041 continue                                                             zipper
      goto 21961                                                           zipper
21962 continue                                                             zipper
      goto 21831                                                           zipper
21832 continue                                                             zipper
      if(jx.gt.0)goto 21762                                                zipper
22050 do 22051 ic=1,nc                                                     zipper
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ix=1                                zipper
      if(ix .ne. 0)goto 22071                                              zipper
22080 do 22081 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22101                   zipper
      ix=1                                                                 zipper
      goto 22082                                                           zipper
22101 continue                                                             zipper
22081 continue                                                             zipper
22082 continue                                                             zipper
22071 continue                                                             zipper
22110 do 22111 i=1,no                                                      zipper
      fi=b(0,ic)+g(i,ic)                                                   zipper
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin),ic),x(i,m(1:nin)))         zipper
      fi=min(max(exmn,fi),exmx)                                            zipper
      sxp(i)=sxp(i)-q(i,ic)                                                zipper
      q(i,ic)=min(max(emin*sxp(i),exp(fi)),emax*sxp(i))                    zipper
      sxp(i)=sxp(i)+q(i,ic)                                                zipper
22111 continue                                                             zipper
22112 continue                                                             zipper
22051 continue                                                             zipper
22052 continue                                                             zipper
      s=-sum(b(0,:))/nc                                                    zipper
      b(0,:)=b(0,:)+s                                                      zipper
      if(jx.gt.0)goto 21762                                                zipper
      if(ix .ne. 0)goto 22131                                              zipper
22140 do 22141 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 22141                                            zipper
      if(ju(k).eq.0)goto 22141                                             zipper
      ga(k)=0.0                                                            zipper
22141 continue                                                             zipper
22142 continue                                                             zipper
22150 do 22151 ic=1,nc                                                     zipper
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      zipper
22160 do 22161 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 22161                                            zipper
      if(ju(k).eq.0)goto 22161                                             zipper
      ga(k)=ga(k)+dot_product(r(:,ic),x(:,k))**2                           zipper
22161 continue                                                             zipper
22162 continue                                                             zipper
22151 continue                                                             zipper
22152 continue                                                             zipper
      ga=sqrt(ga)                                                          zipper
22170 do 22171 k=1,ni                                                      zipper
      if(ixx(k).eq.1)goto 22171                                            zipper
      if(ju(k).eq.0)goto 22171                                             zipper
      if(ga(k) .le. al1*vp(k))goto 22191                                   zipper
      ixx(k)=1                                                             zipper
      ix=1                                                                 zipper
22191 continue                                                             zipper
22171 continue                                                             zipper
22172 continue                                                             zipper
      if(ix.eq.1) go to 10880                                              zipper
      goto 21762                                                           zipper
22131 continue                                                             zipper
      goto 21761                                                           zipper
21762 continue                                                             zipper
      if(kx .le. 0)goto 22211                                              zipper
      jerr=-20000-ilm                                                      zipper
      goto 21682                                                           zipper
22211 continue                                                             zipper
      if(jx .le. 0)goto 22231                                              zipper
      jerr=-10000-ilm                                                      zipper
      goto 21682                                                           zipper
22231 continue                                                             zipper
      devi=0.0                                                             zipper
22240 do 22241 ic=1,nc                                                     zipper
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          zipper
      a0(ic,ilm)=b(0,ic)                                                   zipper
22250 do 22251 i=1,no                                                      zipper
      if(y(i,ic).le.0.0)goto 22251                                         zipper
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           zipper
22251 continue                                                             zipper
22252 continue                                                             zipper
22241 continue                                                             zipper
22242 continue                                                             zipper
      kin(ilm)=nin                                                         zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      dev(ilm)=(dev1-devi)/dev0                                            zipper
      if(ilm.lt.mnl)goto 21681                                             zipper
      if(flmin.ge.1.0)goto 21681                                           zipper
      me=0                                                                 zipper
22260 do 22261 j=1,nin                                                     zipper
      if(a(j,1,ilm).ne.0.0) me=me+1                                        zipper
22261 continue                                                             zipper
22262 continue                                                             zipper
      if(me.gt.ne)goto 21682                                               zipper
      if(dev(ilm).gt.devmax)goto 21682                                     zipper
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 21682                             zipper
21681 continue                                                             zipper
21682 continue                                                             zipper
      g=log(q)                                                             zipper
22270 do 22271 i=1,no                                                      zipper
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         zipper
22271 continue                                                             zipper
22272 continue                                                             zipper
      deallocate(sxp,b,bs,r,q,mm,is,ga,ixx,gk,del,sxpl)                    zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine multsprlognetn(parm,no,ni,nc,x,ix,jx,y,g,w,ju,vp,cl,ne,   zipper
     *nx,nlam,  flmin,ulam,shri,intr,maxit,xv,xb,xs,lmu,a0,a,m,kin,dev0,
     *dev,alm,nlp,jerr)
      implicit double precision(a-h,o-z)                                   zipper
      double precision x(*),y(no,nc),g(no,nc),w(no),vp(ni)                 zipper
      double precision ulam(nlam),xb(ni),xs(ni),xv(ni)                     zipper
      double precision a(nx,nc,nlam),a0(nc,nlam),dev(nlam),alm(nlam),cl(   zipper
     *2,ni)
      integer ix(*),jx(*),ju(ni),m(nx),kin(nlam)                           zipper
      double precision, dimension (:,:), allocatable :: q,r,b,bs                
      double precision, dimension (:), allocatable :: sxp,sxpl,ga,gk            
      double precision, dimension (:), allocatable :: del,sc,svr                
      integer, dimension (:), allocatable :: mm,is,iy,isc                       
      allocate(b(0:ni,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return                                                      
      allocate(bs(0:ni,1:nc),stat=jerr)                                         
      if(jerr.ne.0) return					                                                 
      allocate(q(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return					                                                 
      allocate(r(1:no,1:nc),stat=jerr)                                          
      if(jerr.ne.0) return					                                                 
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx)               zipper
      exmn=-exmx                                                           zipper
      allocate(mm(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(ga(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(gk(1:nc),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(del(1:nc),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(iy(1:ni),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(is(1:max(nc,ni)),stat=jerr)                                 zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sxp(1:no),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sxpl(1:no),stat=jerr)                                       zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(svr(1:nc),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(sc(1:no),stat=jerr)                                         zipper
      if(jerr.ne.0) return                                                 zipper
      allocate(isc(1:nc),stat=jerr)                                        zipper
      if(jerr.ne.0) return                                                 zipper
      pmax=1.0-pmin                                                        zipper
      emin=pmin/pmax                                                       zipper
      emax=1.0/emin                                                        zipper
      bta=parm                                                             zipper
      omb=1.0-bta                                                          zipper
      dev1=0.0                                                             zipper
      dev0=0.0                                                             zipper
22280 do 22281 ic=1,nc                                                     zipper
      q0=dot_product(w,y(:,ic))                                            zipper
      if(q0 .gt. pmin)goto 22301                                           zipper
      jerr =8000+ic                                                        zipper
      return                                                               zipper
22301 continue                                                             zipper
      if(q0 .lt. pmax)goto 22321                                           zipper
      jerr =9000+ic                                                        zipper
      return                                                               zipper
22321 continue                                                             zipper
      b(1:ni,ic)=0.0                                                       zipper
      if(intr .ne. 0)goto 22341                                            zipper
      q0=1.0/nc                                                            zipper
      b(0,ic)=0.0                                                          zipper
      goto 22351                                                           zipper
22341 continue                                                             zipper
      b(0,ic)=log(q0)                                                      zipper
      dev1=dev1-q0*b(0,ic)                                                 zipper
22351 continue                                                             zipper
22331 continue                                                             zipper
22281 continue                                                             zipper
22282 continue                                                             zipper
      if(intr.eq.0) dev1=log(float(nc))                                    zipper
      iy=0                                                                 zipper
      al=0.0                                                               zipper
      if(nonzero(no*nc,g) .ne. 0)goto 22371                                zipper
      b(0,:)=b(0,:)-sum(b(0,:))/nc                                         zipper
      sxp=0.0                                                              zipper
22380 do 22381 ic=1,nc                                                     zipper
      q(:,ic)=exp(b(0,ic))                                                 zipper
      sxp=sxp+q(:,ic)                                                      zipper
22381 continue                                                             zipper
22382 continue                                                             zipper
      goto 22391                                                           zipper
22371 continue                                                             zipper
22400 do 22401 i=1,no                                                      zipper
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         zipper
22401 continue                                                             zipper
22402 continue                                                             zipper
      sxp=0.0                                                              zipper
      if(intr .ne. 0)goto 22421                                            zipper
      b(0,:)=0.0                                                           zipper
      goto 22431                                                           zipper
22421 continue                                                             zipper
      call kazero(nc,no,y,g,w,b(0,:),jerr)                                 zipper
      if(jerr.ne.0) return                                                 zipper
22431 continue                                                             zipper
22411 continue                                                             zipper
      dev1=0.0                                                             zipper
22440 do 22441 ic=1,nc                                                     zipper
      q(:,ic)=b(0,ic)+g(:,ic)                                              zipper
      dev1=dev1-dot_product(w,y(:,ic)*q(:,ic))                             zipper
      q(:,ic)=exp(q(:,ic))                                                 zipper
      sxp=sxp+q(:,ic)                                                      zipper
22441 continue                                                             zipper
22442 continue                                                             zipper
      sxpl=w*log(sxp)                                                      zipper
22450 do 22451 ic=1,nc                                                     zipper
      dev1=dev1+dot_product(y(:,ic),sxpl)                                  zipper
22451 continue                                                             zipper
22452 continue                                                             zipper
22391 continue                                                             zipper
22361 continue                                                             zipper
22460 do 22461 ic=1,nc                                                     zipper
22470 do 22471 i=1,no                                                      zipper
      if(y(i,ic).gt.0.0) dev0=dev0+w(i)*y(i,ic)*log(y(i,ic))               zipper
22471 continue                                                             zipper
22472 continue                                                             zipper
22461 continue                                                             zipper
22462 continue                                                             zipper
      dev0=dev0+dev1                                                       zipper
      alf=1.0                                                              zipper
      if(flmin .ge. 1.0)goto 22491                                         zipper
      eqs=max(eps,flmin)                                                   zipper
      alf=eqs**(1.0/(nlam-1))                                              zipper
22491 continue                                                             zipper
      m=0                                                                  zipper
      mm=0                                                                 zipper
      nin=0                                                                zipper
      nlp=0                                                                zipper
      mnl=min(mnlam,nlam)                                                  zipper
      bs=0.0                                                               zipper
      shr=shri*dev0                                                        zipper
      ga=0.0                                                               zipper
22500 do 22501 ic=1,nc                                                     zipper
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      zipper
      svr(ic)=sum(r(:,ic))                                                 zipper
22510 do 22511 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 22511                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             zipper
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            zipper
22511 continue                                                             zipper
22512 continue                                                             zipper
22501 continue                                                             zipper
22502 continue                                                             zipper
      ga=sqrt(ga)                                                          zipper
22520 do 22521 ilm=1,nlam                                                  zipper
      al0=al                                                               zipper
      if(flmin .lt. 1.0)goto 22541                                         zipper
      al=ulam(ilm)                                                         zipper
      goto 22531                                                           zipper
22541 if(ilm .le. 2)goto 22551                                             zipper
      al=al*alf                                                            zipper
      goto 22531                                                           zipper
22551 if(ilm .ne. 1)goto 22561                                             zipper
      al=big                                                               zipper
      goto 22571                                                           zipper
22561 continue                                                             zipper
      al0=0.0                                                              zipper
22580 do 22581 j=1,ni                                                      zipper
      if(ju(j).eq.0)goto 22581                                             zipper
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))                            zipper
22581 continue                                                             zipper
22582 continue                                                             zipper
      al0=al0/max(bta,1.0d-3)                                              zipper
      al=alf*al0                                                           zipper
22571 continue                                                             zipper
22531 continue                                                             zipper
      al2=al*omb                                                           zipper
      al1=al*bta                                                           zipper
      tlam=bta*(2.0*al-al0)                                                zipper
22590 do 22591 k=1,ni                                                      zipper
      if(iy(k).eq.1)goto 22591                                             zipper
      if(ju(k).eq.0)goto 22591                                             zipper
      if(ga(k).gt.tlam*vp(k)) iy(k)=1                                      zipper
22591 continue                                                             zipper
22592 continue                                                             zipper
10880 continue                                                             zipper
22600 continue                                                             zipper
22601 continue                                                             zipper
      ixx=0                                                                zipper
      jxx=ixx                                                              zipper
      kxx=jxx                                                              zipper
      t=0.0                                                                zipper
22610 do 22611 ic=1,nc                                                     zipper
      t=max(t,maxval(q(:,ic)*(1.0-q(:,ic)/sxp)/sxp))                       zipper
22611 continue                                                             zipper
22612 continue                                                             zipper
      if(t .ge. eps)goto 22631                                             zipper
      kxx=1                                                                zipper
      goto 22602                                                           zipper
22631 continue                                                             zipper
      t=2.0*t                                                              zipper
      alt=al1/t                                                            zipper
      al2t=al2/t                                                           zipper
22640 do 22641 ic=1,nc                                                     zipper
      bs(0,ic)=b(0,ic)                                                     zipper
      if(nin.gt.0) bs(m(1:nin),ic)=b(m(1:nin),ic)                          zipper
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)/t                                    zipper
      svr(ic)=sum(r(:,ic))                                                 zipper
      if(intr .eq. 0)goto 22661                                            zipper
      b(0,ic)=b(0,ic)+svr(ic)                                              zipper
      r(:,ic)=r(:,ic)-svr(ic)*w                                            zipper
      dlx=max(dlx,svr(ic)**2)                                              zipper
22661 continue                                                             zipper
22641 continue                                                             zipper
22642 continue                                                             zipper
22670 continue                                                             zipper
22671 continue                                                             zipper
      nlp=nlp+nc                                                           zipper
      dlx=0.0                                                              zipper
22680 do 22681 k=1,ni                                                      zipper
      if(iy(k).eq.0)goto 22681                                             zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      del=b(k,:)                                                           zipper
      gkn=0.0                                                              zipper
22690 do 22691 ic=1,nc                                                     zipper
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))-svr(ic)*xb(k))/xs(k)        zipper
      gk(ic)=u+del(ic)*xv(k)                                               zipper
      gkn=gkn+gk(ic)**2                                                    zipper
22691 continue                                                             zipper
22692 continue                                                             zipper
      gkn=sqrt(gkn)                                                        zipper
      u=1.0-alt*vp(k)/gkn                                                  zipper
      if(u .gt. 0.0)goto 22711                                             zipper
      b(k,:)=0.0                                                           zipper
      goto 22721                                                           zipper
22711 continue                                                             zipper
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     zipper
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   zipper
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 zipper
22721 continue                                                             zipper
22701 continue                                                             zipper
      del=b(k,:)-del                                                       zipper
      if(maxval(abs(del)).le.0.0)goto 22681                                zipper
22730 do 22731 ic=1,nc                                                     zipper
      dlx=max(dlx,xv(k)*del(ic)**2)                                        zipper
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   zipper
     *b(k))/xs(k)
22731 continue                                                             zipper
22732 continue                                                             zipper
      if(mm(k) .ne. 0)goto 22751                                           zipper
      nin=nin+1                                                            zipper
      if(nin .le. nx)goto 22771                                            zipper
      jxx=1                                                                zipper
      goto 22682                                                           zipper
22771 continue                                                             zipper
      mm(k)=nin                                                            zipper
      m(nin)=k                                                             zipper
22751 continue                                                             zipper
22681 continue                                                             zipper
22682 continue                                                             zipper
      if(jxx.gt.0)goto 22672                                               zipper
      if(dlx.lt.shr)goto 22672                                             zipper
      if(nlp .le. maxit)goto 22791                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
22791 continue                                                             zipper
22800 continue                                                             zipper
22801 continue                                                             zipper
      nlp=nlp+nc                                                           zipper
      dlx=0.0                                                              zipper
22810 do 22811 l=1,nin                                                     zipper
      k=m(l)                                                               zipper
      jb=ix(k)                                                             zipper
      je=ix(k+1)-1                                                         zipper
      del=b(k,:)                                                           zipper
      gkn=0.0                                                              zipper
22820 do 22821 ic=1,nc                                                     zipper
      u=(dot_product(r(jx(jb:je),ic),x(jb:je))  -svr(ic)*xb(k))/xs(k)      zipper
      gk(ic)=u+del(ic)*xv(k)                                               zipper
      gkn=gkn+gk(ic)**2                                                    zipper
22821 continue                                                             zipper
22822 continue                                                             zipper
      gkn=sqrt(gkn)                                                        zipper
      u=1.0-alt*vp(k)/gkn                                                  zipper
      if(u .gt. 0.0)goto 22841                                             zipper
      b(k,:)=0.0                                                           zipper
      goto 22851                                                           zipper
22841 continue                                                             zipper
      b(k,:)=gk*(u/(xv(k)+vp(k)*al2t))                                     zipper
      call chkbnds1(nc,gk,gkn,xv(k),cl(1,k),cl(2,k),  vp(k)*al2t,alt*vp(   zipper
     *k),b(k,:),isc,jerr)
      if(jerr.ne.0) return                                                 zipper
22851 continue                                                             zipper
22831 continue                                                             zipper
      del=b(k,:)-del                                                       zipper
      if(maxval(abs(del)).le.0.0)goto 22811                                zipper
22860 do 22861 ic=1,nc                                                     zipper
      dlx=max(dlx,xv(k)*del(ic)**2)                                        zipper
      r(jx(jb:je),ic)=r(jx(jb:je),ic)  -del(ic)*w(jx(jb:je))*(x(jb:je)-x   zipper
     *b(k))/xs(k)
22861 continue                                                             zipper
22862 continue                                                             zipper
22811 continue                                                             zipper
22812 continue                                                             zipper
      if(dlx.lt.shr)goto 22802                                             zipper
      if(nlp .le. maxit)goto 22881                                         zipper
      jerr=-ilm                                                            zipper
      return                                                               zipper
22881 continue                                                             zipper
      goto 22801                                                           zipper
22802 continue                                                             zipper
      goto 22671                                                           zipper
22672 continue                                                             zipper
      if(jxx.gt.0)goto 22602                                               zipper
22890 do 22891 ic=1,nc                                                     zipper
      if((b(0,ic)-bs(0,ic))**2.gt.shr) ixx=1                               zipper
      if(ixx .ne. 0)goto 22911                                             zipper
22920 do 22921 j=1,nin                                                     zipper
      k=m(j)                                                               zipper
      if(xv(k)*(b(k,ic)-bs(k,ic))**2 .le. shr)goto 22941                   zipper
      ixx=1                                                                zipper
      goto 22922                                                           zipper
22941 continue                                                             zipper
22921 continue                                                             zipper
22922 continue                                                             zipper
22911 continue                                                             zipper
      sc=b(0,ic)+g(:,ic)                                                   zipper
      b0=0.0                                                               zipper
22950 do 22951 j=1,nin                                                     zipper
      l=m(j)                                                               zipper
      jb=ix(l)                                                             zipper
      je=ix(l+1)-1                                                         zipper
      sc(jx(jb:je))=sc(jx(jb:je))+b(l,ic)*x(jb:je)/xs(l)                   zipper
      b0=b0-b(l,ic)*xb(l)/xs(l)                                            zipper
22951 continue                                                             zipper
22952 continue                                                             zipper
      sc=min(max(exmn,sc+b0),exmx)                                         zipper
      sxp=sxp-q(:,ic)                                                      zipper
      q(:,ic)=min(max(emin*sxp,exp(sc)),emax*sxp)                          zipper
      sxp=sxp+q(:,ic)                                                      zipper
22891 continue                                                             zipper
22892 continue                                                             zipper
      s=sum(b(0,:))/nc                                                     zipper
      b(0,:)=b(0,:)-s                                                      zipper
      if(jxx.gt.0)goto 22602                                               zipper
      if(ixx .ne. 0)goto 22971                                             zipper
22980 do 22981 j=1,ni                                                      zipper
      if(iy(j).eq.1)goto 22981                                             zipper
      if(ju(j).eq.0)goto 22981                                             zipper
      ga(j)=0.0                                                            zipper
22981 continue                                                             zipper
22982 continue                                                             zipper
22990 do 22991 ic=1,nc                                                     zipper
      r(:,ic)=w*(y(:,ic)-q(:,ic)/sxp)                                      zipper
23000 do 23001 j=1,ni                                                      zipper
      if(iy(j).eq.1)goto 23001                                             zipper
      if(ju(j).eq.0)goto 23001                                             zipper
      jb=ix(j)                                                             zipper
      je=ix(j+1)-1                                                         zipper
      gj=dot_product(r(jx(jb:je),ic),x(jb:je))                             zipper
      ga(j)=ga(j)+((gj-svr(ic)*xb(j))/xs(j))**2                            zipper
23001 continue                                                             zipper
23002 continue                                                             zipper
22991 continue                                                             zipper
22992 continue                                                             zipper
      ga=sqrt(ga)                                                          zipper
23010 do 23011 k=1,ni                                                      zipper
      if(iy(k).eq.1)goto 23011                                             zipper
      if(ju(k).eq.0)goto 23011                                             zipper
      if(ga(k) .le. al1*vp(k))goto 23031                                   zipper
      iy(k)=1                                                              zipper
      ixx=1                                                                zipper
23031 continue                                                             zipper
23011 continue                                                             zipper
23012 continue                                                             zipper
      if(ixx.eq.1) go to 10880                                             zipper
      goto 22602                                                           zipper
22971 continue                                                             zipper
      goto 22601                                                           zipper
22602 continue                                                             zipper
      if(kxx .le. 0)goto 23051                                             zipper
      jerr=-20000-ilm                                                      zipper
      goto 22522                                                           zipper
23051 continue                                                             zipper
      if(jxx .le. 0)goto 23071                                             zipper
      jerr=-10000-ilm                                                      zipper
      goto 22522                                                           zipper
23071 continue                                                             zipper
      devi=0.0                                                             zipper
23080 do 23081 ic=1,nc                                                     zipper
      if(nin.gt.0) a(1:nin,ic,ilm)=b(m(1:nin),ic)                          zipper
      a0(ic,ilm)=b(0,ic)                                                   zipper
23090 do 23091 i=1,no                                                      zipper
      if(y(i,ic).le.0.0)goto 23091                                         zipper
      devi=devi-w(i)*y(i,ic)*log(q(i,ic)/sxp(i))                           zipper
23091 continue                                                             zipper
23092 continue                                                             zipper
23081 continue                                                             zipper
23082 continue                                                             zipper
      kin(ilm)=nin                                                         zipper
      alm(ilm)=al                                                          zipper
      lmu=ilm                                                              zipper
      dev(ilm)=(dev1-devi)/dev0                                            zipper
      if(ilm.lt.mnl)goto 22521                                             zipper
      if(flmin.ge.1.0)goto 22521                                           zipper
      me=0                                                                 zipper
23100 do 23101 j=1,nin                                                     zipper
      if(a(j,1,ilm).ne.0.0) me=me+1                                        zipper
23101 continue                                                             zipper
23102 continue                                                             zipper
      if(me.gt.ne)goto 22522                                               zipper
      if(dev(ilm).gt.devmax)goto 22522                                     zipper
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 22522                             zipper
22521 continue                                                             zipper
22522 continue                                                             zipper
      g=log(q)                                                             zipper
23110 do 23111 i=1,no                                                      zipper
      g(i,:)=g(i,:)-sum(g(i,:))/nc                                         zipper
23111 continue                                                             zipper
23112 continue                                                             zipper
      deallocate(sxp,b,bs,r,q,mm,is,sc,ga,iy,gk,del,sxpl)                  zipper
      return                                                               zipper
      end                                                                  zipper
      subroutine psort7 (v,a,ii,jj)                                             
      implicit double precision(a-h,o-z)                                        
c                                                                               
c     puts into a the permutation vector which sorts v into                     
c     increasing order. the array v is not modified.                            
c     only elements from ii to jj are considered.                               
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements           
c                                                                               
c     this is a modification of cacm algorithm #347 by r. c. singleton,         
c     which is a modified hoare quicksort.                                      
c                                                                               
      dimension a(jj),v(jj),iu(20),il(20)                                       
      integer t,tt                                                              
      integer a                                                                 
      double precision v                                                        
      m=1                                                                       
      i=ii                                                                      
      j=jj                                                                      
 10   if (i.ge.j) go to 80                                                      
 20   k=i                                                                       
      ij=(j+i)/2                                                                
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 30                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
 30   l=j                                                                       
      if (v(a(j)).ge.vt) go to 50                                               
      a(ij)=a(j)                                                                
      a(j)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 50                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      go to 50                                                                  
 40   a(l)=a(k)                                                                 
      a(k)=tt                                                                   
 50   l=l-1                                                                     
      if (v(a(l)).gt.vt) go to 50                                               
      tt=a(l)                                                                   
      vtt=v(tt)                                                                 
 60   k=k+1                                                                     
      if (v(a(k)).lt.vt) go to 60                                               
      if (k.le.l) go to 40                                                      
      if (l-i.le.j-k) go to 70                                                  
      il(m)=i                                                                   
      iu(m)=l                                                                   
      i=k                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 70   il(m)=k                                                                   
      iu(m)=j                                                                   
      j=l                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 80   m=m-1                                                                     
      if (m.eq.0) return                                                        
      i=il(m)                                                                   
      j=iu(m)                                                                   
 90   if (j-i.gt.10) go to 20                                                   
      if (i.eq.ii) go to 10                                                     
      i=i-1                                                                     
 100  i=i+1                                                                     
      if (i.eq.j) go to 80                                                      
      t=a(i+1)                                                                  
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 100                                              
      k=i                                                                       
 110  a(k+1)=a(k)                                                               
      k=k-1                                                                     
      if (vt.lt.v(a(k))) go to 110                                              
      a(k+1)=t                                                                  
      go to 100                                                                 
      end                                                                       
