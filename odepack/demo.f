c-----------------------------------------------------------------------
c demonstration program for the lsode package.
c this is the version of 13 august, 1981.
c
c this version is in double precision.
c
c for computer systems requiring a program card, the following (with
c the c in column 1 removed) may be used..
c     program lsdem(lsout,tape6=lsout)
c
c the package is used to solve two simple problems,
c one with a full jacobian, the other with a banded jacobian,
c with all 8 of the appropriate values of mf in each case.
c if the errors are too large, or other difficulty occurs,
c a warning message is printed.  all output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, jac1, f2, jac2
      integer i, iopar, iopt, iout, istate, itask, itol, iwork,
     1   leniw, lenrw, liw, lout, lrw, mband, meth, mf, miter, miter1,
     2   ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst
      double precision atol, dtout, er, erm, ero, hu, rtol, rwork, t,
     1   tout, tout1, y
      dimension y(25), rwork(697), iwork(45)
      data lout/6/, tout1/1.39283880203d0/, dtout/2.214773875d0/
c
      nerr = 0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-6
      lrw = 697
      liw = 45
      iopt = 0
c
c first problem
c
      neq = 2
      nout = 4
      write (lout,110) neq,itol,rtol,atol
 110  format(1h1/1x,40h demonstration program for lsode package////
     1  1x,39h problem 1..   van der pol oscillator../
     2  1x,38h  xdotdot - 3*(1 - x**2)*xdot + x = 0, ,
     3  24h   x(0) = 2, xdot(0) = 0/
     4  1x,6h neq =,i2/
     5  1x,7h itol =,i3,9h   rtol =,e10.1,9h   atol =,e10.1//)
c
      do 195 meth = 1,2
      do 190 miter1 = 1,4
      miter = miter1 - 1
      mf = 10*meth + miter
      write (lout,120) mf
 120  format(////1x,5h mf =,i3///
     1  6x,1ht,15x,1hx,15x,4hxdot,7x,2hnq,6x,1hh//)
      t = 0.0d0
      y(1) = 2.0d0
      y(2) = 0.0d0
      itask = 1
      istate = 1
      tout = tout1
      ero = 0.0d0
      do 170 iout = 1,nout
        call lsode(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac1,mf)
        hu = rwork(11)
        nqu = iwork(14)
        write (lout,140) t,y(1),y(2),nqu,hu
 140    format(1x,e15.5,e16.5,e14.3,i5,e14.3)
        if (istate .lt. 0) go to 175
        iopar = iout - 2*(iout/2)
        if (iopar .ne. 0) go to 170
        er = dabs(y(1))/atol
        ero = dmax1(ero,er)
        if (er .lt. 1000.0d0) go to 170
        write (lout,150)
 150    format(//1x,41h warning.. error exceeds 1000 * tolerance//)
        nerr = nerr + 1
 170    tout = tout + dtout
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (miter .eq. 2) nfea = nfe - neq*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(//1x,32h final statistics for this run../
     1  1x,13h rwork size =,i4,15h   iwork size =,i4/
     2  1x,18h number of steps =,i5/
     3  1x,18h number of f-s   =,i5/
     4  1x,18h (excluding j-s) =,i5/
     5  1x,18h number of j-s   =,i5/
     6  1x,16h error overrun =,e10.2)
 190  continue
 195  continue
c
c second problem
c
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      nout = 5
      write (lout,210) neq,ml,mu,itol,rtol,atol
 210  format(1h1/1x,40h demonstration program for lsode package////
     1  1x,33h problem 2.. ydot = a * y , where,
     2  39h  a is a banded lower triangular matrix/
     3  1x,32h  derived from 2-d advection pde/
     4  1x,6h neq =,i3,7h   ml =,i2,7h   mu =,i2/
     5  1x,7h itol =,i3,9h   rtol =,e10.1,9h   atol =,e10.1//)
      do 295 meth = 1,2
      do 290 miter1 = 1,6
      miter = miter1 - 1
      if (miter .eq. 1 .or. miter .eq. 2) go to 290
      mf = 10*meth + miter
      write (lout,220) mf
 220  format(////1x,5h mf =,i3///
     1  6x,1ht,13x,8hmax.err.,5x,2hnq,6x,1hh//)
      t = 0.0d0
      do 230 i = 2,neq
 230    y(i) = 0.0d0
      y(1) = 1.0d0
      itask = 1
      istate = 1
      tout = 0.01d0
      ero = 0.0d0
      do 270 iout = 1,nout
        call lsode(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac2,mf)
        call edit2(y,t,erm)
        hu = rwork(11)
        nqu = iwork(14)
        write (lout,240) t,erm,nqu,hu
 240    format(1x,e15.5,e14.3,i5,e14.3)
        if (istate .lt. 0) go to 275
        er = erm/atol
        ero = dmax1(ero,er)
        if (er .le. 1000.0d0) go to 270
        write (lout,150)
        nerr = nerr + 1
 270    tout = tout*10.0d0
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (miter .eq. 5) nfea = nfe - mband*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 290  continue
 295  continue
      write (lout,300) nerr
 300  format(////1x,31h number of errors encountered =,i3)
      stop
      end
      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(2), ydot(2)
      ydot(1) = y(2)
      ydot(2) = 3.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end
      subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(2), pd(nrowpd,2)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -6.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 3.0d0*(1.0d0 - y(1)*y(1))
      return
      end
      subroutine f2 (neq, t, y, ydot)
      integer neq, i, j, k, ng
      double precision t, y, ydot, alph1, alph2, d
      dimension y(1), ydot(1)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      do 10 j = 1,ng
      do 10 i = 1,ng
        k = i + (j - 1)*ng
        d = -2.0d0*y(k)
        if (i .ne. 1) d = d + y(k-1)*alph1
        if (j .ne. 1) d = d + y(k-ng)*alph2
 10     ydot(k) = d
      return
      end
      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd, j, mband, mu1, mu2, ng
      double precision t, y, pd, alph1, alph2
      dimension y(1), pd(nrowpd,1)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do 10 j = 1,neq
        pd(mu1,j) = -2.0d0
        pd(mu2,j) = alph1
 10     pd(mband,j) = alph2
      do 20 j = ng,neq,ng
 20     pd(mu2,j) = 0.0d0
      return
      end
      subroutine edit2 (y, t, erm)
      integer i, j, k, ng
      double precision y, t, erm, alph1, alph2, a1, a2, er, ex, yt
      dimension y(25)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      erm = 0.0d0
      if (t .eq. 0.0d0) return
      ex = 0.0d0
      if (t .le. 30.0d0) ex = dexp(-2.0d0*t)
      a2 = 1.0d0
      do 60 j = 1,ng
        a1 = 1.0d0
        do 50 i = 1,ng
          k = i + (j - 1)*ng
          yt = t**(i+j-2)*ex*a1*a2
          er = dabs(y(k)-yt)
          erm = dmax1(erm,er)
          a1 = a1*alph1/dfloat(i)
 50       continue
        a2 = a2*alph2/dfloat(j)
 60     continue
      return
      end
c-----------------------------------------------------------------------
c demonstration program for the lsodes package.
c this is the version of may 3, 1983.
c
c this version is in double precision.
c
c for computer systems requiring a program card, the following (with
c the c in column 1 removed) may be used..
c
c     program lssdem(lssout,tape6=lssout)
c
c the package is used for each of the relevant values of mf to solve
c the problem ydot = a * y, where a is the 9 by 9 sparse matrix
c
c               -4  1     1
c                1 -4  1     1
c                   1 -4        1
c                        -4  1     1
c       a =               1 -4  1     1
c                            1 -4        1
c                                 -4  1
c                                  1 -4  1
c                                     1 -4
c
c the initial conditions are  y(0) = (1, 2, 3, ..., 9).
c output is printed at t = 1., 2., and 3.
c each case is solved first with nominal (large) values of lrw and liw,
c and then with values given by lenrw and leniw (optional outputs)
c on the first run, as a check on these computed work array lengths.
c if the errors are too large, or other difficulty occurs,
c a warning message is printed.
c all output is on unit lout, which is data-loaded to 6 below.
c-----------------------------------------------------------------------
      external fdem, jdem
      integer i, ia, igrid, iopt, iout, irun, istate, itask, itol,
     1  iwork, j, ja, k, l, leniw, lenrw, liw, lout, lrw,
     2  m, meth, mf, miter, miter1, moss, moss1, neq, nerr, nfe, nfea,
     3  ngp, nje, nlu, nnz, nout, nqu, nst, nzl, nzu
      double precision atol, erm, ero, hu, rtol, rwork, t, tout, y
      dimension y(9), ia(10), ja(50), iwork(90), rwork(1000)
      equivalence (ia(1),iwork(31)), (ja(1),iwork(41))
      data lout/6/
c
c write heading and set fixed parameters.
      write(lout,10)
 10   format(1x/45h demonstration problem for the lsodes package////)
      nerr = 0
      igrid = 3
      neq = igrid**2
      t = 0.0d0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-5
      itask = 1
      iopt = 0
      do 20 i = 1,neq
 20     y(i) = dfloat(i)
      ia(1) = 1
      k = 1
      do 60 m = 1,igrid
        do 50 l = 1,igrid
          j = l + (m - 1)*igrid
          if (m .eq. 1) go to 30
          ja(k) = j - igrid
          k = k + 1
 30       if (l .eq. 1) go to 35
          ja(k) = j - 1
          k = k + 1
 35       ja(k) = j
          k = k + 1
          if (l .eq. igrid) go to 40
          ja(k) = j + 1
          k = k + 1
 40       ia(j+1) = k
 50       continue
 60     continue
      write (lout,80)neq,t,rtol,atol,(y(i),i=1,neq)
 80   format(1x,6h neq =,i4,5x,4ht0 =,f4.1,5x,6hrtol =,e12.3,5x,
     1   6hatol =,e12.3//1x,21h initial y vector =  ,9f5.1////)
c
c loop over all relevant values of mf.
      do 193 moss1 = 1,3
      moss = moss1 - 1
      do 192 meth = 1,2
      do 191 miter1 = 1,4
      miter = miter1 - 1
      if ( (miter.eq.0 .or. miter.eq.3) .and. moss.ne.0) go to 191
      mf = 100*moss + 10*meth + miter
      if (mf .ne. 10) write (lout,100)
 100  format(1h1)
c first run.. nominal work array lengths, 3 output points.
      irun = 1
      lrw = 1000
      liw = 90
      nout = 3
 110  continue
      write (lout,120)mf,lrw,liw
 120  format(///1x,4hmf =,i4,5x,29hinput work lengths lrw, liw =,2i6)
      do 125 i = 1,neq
 125    y(i) = dfloat(i)
      t = 0.0d0
      tout = 1.0d0
      istate = 1
      ero = 0.0d0
c loop over output points.  do output and accuracy check at each.
      do 170 iout = 1,nout
        call lsodes (fdem, neq, y, t, tout, itol, rtol, atol,
     1     itask, istate, iopt, rwork, lrw, iwork, liw, jdem, mf)
        nst = iwork(11)
        hu = rwork(11)
        nqu = iwork(14)
        call edit (y, iout, erm)
        write(lout,140)t,nst,hu,nqu,erm,(y(i),i=1,neq)
 140    format(//1x,7h at t =,f5.1,5x,5hnst =,i4,5x,4hhu =,e12.3,5x,
     1    5hnqu =,i3,5x,12h max. err. =,e12.4/
     2    1x,15h  y array =    ,4e15.6/1x,5e15.6)
        if (istate .lt. 0) go to 175
        erm = erm/atol
        ero = dmax1(ero,erm)
        if (erm .lt. 100.0d0) go to 160
        write (lout,150)
 150    format(//1x,40h warning.. error exceeds 100 * tolerance//)
        nerr = nerr + 1
 160    tout = tout + 1.0d0
 170    continue
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      if (irun .eq. 2) go to 191
c print final statistics (first run only)
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nnz = iwork(19)
      ngp = iwork(20)
      nlu = iwork(21)
      nzl = iwork(25)
      nzu = iwork(26)
      nfea = nfe
      if (miter .eq. 2) nfea = nfe - ngp*nje
      if (miter .eq. 3) nfea = nfe - nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(//1x,32h final statistics for this run../
     1  1x,13h rwork size =,i4,15h   iwork size =,i4/
     2  1x,18h number of steps =,i5/
     3  1x,18h number of f-s   =,i5/
     4  1x,18h (excluding j-s) =,i5/
     5  1x,18h number of j-s   =,i5/
     6  1x,16h error overrun =,e10.2)
      if (miter .eq. 1 .or. miter .eq. 2)
     1   write (lout,185)nnz,ngp,nlu,nzl,nzu
 185  format(1x,27h number of nonzeros in j = ,i5/
     1  1x,27h number of j index groups =,i5/
     2  1x,27h number of lu decomp-s    =,i5/
     3  1x,34h nonzeros in strict lower factor =,i5/
     4  1x,34h nonzeros in strict upper factor =,i5)
      if (istate .lt. 0) go to 191
      if (miter .eq. 1 .or. miter .eq. 2)
     1   call ssout (neq, rwork(21), iwork, lout)
c return for second run.. minimal work array lengths, 1 output point.
      irun = irun + 1
      lrw = lenrw
      liw = leniw
      nout = 1
      go to 110
 191  continue
 192  continue
 193  continue
c
      write (lout,200) nerr
 200  format(////1x,31h number of errors encountered =,i3)
      stop
      end
      subroutine fdem (neq, t, y, ydot)
      integer neq,  i, igrid, j, l, m
      double precision t, y, ydot
      dimension y(neq), ydot(neq)
      data igrid/3/
      do 5 i = 1,neq
 5      ydot(i) = 0.0d0
      do 20 m = 1,igrid
        do 10 l = 1,igrid
          j = l + (m - 1)*igrid
          if (m .ne. 1) ydot(j-igrid) = ydot(j-igrid) + y(j)
          if (l .ne. 1) ydot(j-1) = ydot(j-1) + y(j)
          ydot(j) = ydot(j) - 4.0d0*y(j)
          if (l .ne. igrid) ydot(j+1) = ydot(j+1) + y(j)
 10       continue
 20     continue
      return
      end
      subroutine jdem (neq, t, y, j, ia, ja, pdj)
      integer neq, j, ia, ja,  igrid, l, m
      double precision t, y, pdj
      dimension y(1), ia(1), ja(1), pdj(1)
      data igrid/3/
      m = (j - 1)/igrid + 1
      l = j - (m - 1)*igrid
      pdj(j) = -4.0d0
      if (m .ne. 1) pdj(j-igrid) = 1.0d0
      if (l .ne. 1) pdj(j-1) = 1.0d0
      if (l .ne. igrid) pdj(j+1) = 1.0d0
      return
      end
      subroutine edit (y, iout, erm)
      integer iout,  i, neq
      double precision y, erm,   er, yex
      dimension y(1),yex(9,3)
      data neq /9/
      data yex /6.687279d-01, 9.901910d-01, 7.603061d-01,
     1   8.077979d-01, 1.170226d+00, 8.810605d-01, 5.013331d-01,
     2   7.201389d-01, 5.379644d-01, 1.340488d-01, 1.917157d-01,
     3   1.374034d-01, 1.007882d-01, 1.437868d-01, 1.028010d-01,
     4   3.844343d-02, 5.477593d-02, 3.911435d-02, 1.929166d-02,
     5   2.735444d-02, 1.939611d-02, 1.055981d-02, 1.496753d-02,
     6   1.060897d-02, 2.913689d-03, 4.128975d-03, 2.925977d-03/
      erm = 0.0d0
      do 10 i = 1,neq
        er = dabs(y(i) - yex(i,iout))
 10     erm = dmax1(erm,er)
      return
      end
      subroutine ssout (neq, iwk, iwork, lout)
      integer neq, iwk, iwork, lout
      integer i, i1, i2, ipian, ipjan, nnz
      dimension iwk(1), iwork(1)
      ipian = iwork(23)
      ipjan = iwork(24)
      nnz = iwork(19)
      i1 = ipian
      i2 = i1 + neq
      write (lout,10)(iwk(i),i=i1,i2)
 10   format(//1x,33h structure descriptor array ian =/(20i4))
      i1 = ipjan
      i2 = i1 + nnz - 1
      write (lout,20)(iwk(i),i=i1,i2)
 20   format(/1x,33h structure descriptor array jan =/(20i4))
      return
      end
c-----------------------------------------------------------------------
c demonstration program for the lsoda package.
c this is the january 26, 1982 version.
c this version is in double precision.
c
c for computer systems requiring a program card, the following (with
c the c in column 1 removed) may be used..
c     program lsadem(laout,tape6=laout)
c
c the package is used to solve two simple problems,
c one with a full jacobian, the other with a banded jacobian,
c with the 2 appropriate values of jt in each case.
c if the errors are too large, or other difficulty occurs,
c a warning message is printed.  all output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, jac1, f2, jac2
      integer i, iopar, iopt, iout, istate, itask, itol, iwork,
     1   jt, leniw, lenrw, liw, lout, lrw, mband, mused,
     2   ml, mu, neq, nerr, nfe, nfea, nje, nout, nqu, nst
      double precision atol, dtout, er, erm, ero, hu, rtol, rwork, t,
     1   tout, tout1, tsw, y
      dimension y(25), rwork(522), iwork(45)
      data lout/6/, tout1/1.39283880203d0/, dtout/2.214773875d0/
c
      nerr = 0
      itol = 1
      rtol = 0.0d0
      atol = 1.0d-6
      lrw = 522
      liw = 45
      iopt = 0
c
c first problem
c
      neq = 2
      nout = 4
      write (lout,110) neq,itol,rtol,atol
 110  format(1h1/1x,40h demonstration program for lsoda package////
     1  1x,39h problem 1..   van der pol oscillator../
     2  1x,38h  xdotdot - 3*(1 - x**2)*xdot + x = 0, ,
     3  24h   x(0) = 2, xdot(0) = 0/
     4  1x,6h neq =,i2/
     5  1x,7h itol =,i3,9h   rtol =,e10.1,9h   atol =,e10.1//)
c
      do 190 jt = 1,2
      write (lout,120) jt
 120  format(////1x,5h jt =,i3///
     1  6x,1ht,15x,1hx,15x,4hxdot,7x,4hmeth,3x,2hnq,6x,1hh,12x,3htsw//)
      t = 0.0d0
      y(1) = 2.0d0
      y(2) = 0.0d0
      itask = 1
      istate = 1
      tout = tout1
      ero = 0.0d0
      do 170 iout = 1,nout
        call lsoda(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac1,jt)
        hu = rwork(11)
        tsw = rwork(15)
        nqu = iwork(14)
        mused = iwork(19)
        write (lout,140) t,y(1),y(2),mused,nqu,hu,tsw
 140    format(1x,e15.5,e16.5,e14.3,2i6,2e14.3)
        if (istate .lt. 0) go to 175
        iopar = iout - 2*(iout/2)
        if (iopar .ne. 0) go to 170
        er = dabs(y(1))/atol
        ero = dmax1(ero,er)
        if (er .lt. 1000.0d0) go to 170
        write (lout,150)
 150    format(//1x,41h warning.. error exceeds 1000 * tolerance//)
        nerr = nerr + 1
 170    tout = tout + dtout
 175  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 180  format(//1x,32h final statistics for this run../
     1  1x,13h rwork size =,i4,15h   iwork size =,i4/
     2  1x,18h number of steps =,i5/
     3  1x,18h number of f-s   =,i5/
     4  1x,18h (excluding j-s) =,i5/
     5  1x,18h number of j-s   =,i5/
     6  1x,16h error overrun =,e10.2)
 190  continue
c
c second problem
c
      neq = 25
      ml = 5
      mu = 0
      iwork(1) = ml
      iwork(2) = mu
      mband = ml + mu + 1
      nout = 5
      write (lout,210) neq,ml,mu,itol,rtol,atol
 210  format(1h1/1x,40h demonstration program for lsoda package////
     1  1x,33h problem 2.. ydot = a * y , where,
     2  39h  a is a banded lower triangular matrix/
     2  1x,32h  derived from 2-d advection pde/
     3  1x,6h neq =,i3,7h   ml =,i2,7h   mu =,i2/
     4  1x,7h itol =,i3,9h   rtol =,e10.1,9h   atol =,e10.1//)
      do 290 jt = 4,5
      write (lout,220) jt
 220  format(////1x,5h jt =,i3///
     1  6x,1ht,13x,8hmax.err.,5x,4hmeth,3x,2hnq,6x,1hh,12x,3htsw//)
      t = 0.0d0
      do 230 i = 2,neq
 230    y(i) = 0.0d0
      y(1) = 1.0d0
      itask = 1
      istate = 1
      tout = 0.01d0
      ero = 0.0d0
      do 270 iout = 1,nout
        call lsoda(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac2,jt)
        call edit2(y,t,erm)
        hu = rwork(11)
        tsw = rwork(15)
        nqu = iwork(14)
        mused = iwork(19)
        write (lout,240) t,erm,mused,nqu,hu,tsw
 240    format(1x,e15.5,e14.3,2i6,2e14.3)
        if (istate .lt. 0) go to 275
        er = erm/atol
        ero = dmax1(ero,er)
        if (er .le. 1000.0d0) go to 270
        write (lout,150)
        nerr = nerr + 1
 270    tout = tout*10.0d0
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 5) nfea = nfe - mband*nje
      write (lout,180) lenrw,leniw,nst,nfe,nfea,nje,ero
 290  continue
      write (lout,300) nerr
 300  format(////1x,31h number of errors encountered =,i3)
      stop
      end
      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(2), ydot(2)
      ydot(1) = y(2)
      ydot(2) = 3.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end
      subroutine jac1 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(2), pd(nrowpd,2)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -6.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 3.0d0*(1.0d0 - y(1)*y(1))
      return
      end
      subroutine f2 (neq, t, y, ydot)
      integer neq, i, j, k, ng
      double precision t, y, ydot, alph1, alph2, d
      dimension y(1), ydot(1)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      do 10 j = 1,ng
      do 10 i = 1,ng
        k = i + (j - 1)*ng
        d = -2.0d0*y(k)
        if (i .ne. 1) d = d + y(k-1)*alph1
        if (j .ne. 1) d = d + y(k-ng)*alph2
 10     ydot(k) = d
      return
      end
      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd, j, mband, mu1, mu2, ng
      double precision t, y, pd, alph1, alph2
      dimension y(1), pd(nrowpd,1)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      mband = ml + mu + 1
      mu1 = mu + 1
      mu2 = mu + 2
      do 10 j = 1,neq
        pd(mu1,j) = -2.0d0
        pd(mu2,j) = alph1
 10     pd(mband,j) = alph2
      do 20 j = ng,neq,ng
 20     pd(mu2,j) = 0.0d0
      return
      end
      subroutine edit2 (y, t, erm)
      integer i, j, k, ng
      double precision y, t, erm, alph1, alph2, a1, a2, er, ex, yt
      dimension y(25)
      data alph1/1.0d0/, alph2/1.0d0/, ng/5/
      erm = 0.0d0
      if (t .eq. 0.0d0) return
      ex = 0.0d0
      if (t .le. 30.0d0) ex = dexp(-2.0d0*t)
      a2 = 1.0d0
      do 60 j = 1,ng
        a1 = 1.0d0
        do 50 i = 1,ng
          k = i + (j - 1)*ng
          yt = t**(i+j-2)*ex*a1*a2
          er = dabs(y(k)-yt)
          erm = dmax1(erm,er)
          a1 = a1*alph1/dfloat(i)
 50       continue
        a2 = a2*alph2/dfloat(j)
 60     continue
      return
      end
c-----------------------------------------------------------------------
c demonstration program for the lsodar package.
c this is the january 27, 1982 version.
c this version is in double precision.
c
c for computer systems requiring a program card, the following (with
c the c in column 1 removed) may be used..
c     program lardem(larout,tape6=larout)
c
c the lsodar package is used to solve two simple problems,
c one nonstiff and one intermittently stiff.
c if the errors are too large, or other difficulty occurs,
c a warning message is printed.  all output is on unit lout = 6.
c-----------------------------------------------------------------------
      external f1, gr1, f2, jac2, gr2
      integer iopt, iout, istate, itask, itol, iwork, jroot, jt,
     1   kroot, leniw, lenrw, liw, lrw, lout, neq, nerr, ng,
     2   nfe, nfea, nge, nje, nst
      double precision atol, er, ero, errt, rtol, rwork,
     1   t, tout, tzero, y, yt
      dimension y(2), atol(2), rwork(57), iwork(22), jroot(2)
      data lout/6/
c
      nerr = 0
c-----------------------------------------------------------------------
c first problem.
c the initial value problem is..
c   dy/dt = ((2*log(y) + 8)/t - 5)*y,  y(1) = 1,  1 .le. t .le. 6
c the solution is  y(t) = exp(-t**2 + 5*t - 4)
c the two root functions are..
c   g1 = ((2*log(y)+8)/t - 5)*y (= dy/dt)  (with root at t = 2.5),
c   g2 = log(y) - 2.2491  (with roots at t = 2.47 and 2.53)
c-----------------------------------------------------------------------
c set all input parameters and print heading.
      neq = 1
      y(1) = 1.0d0
      t = 1.0d0
      tout = 2.0d0
      itol = 1
      rtol = 0.0d0
      atol(1) = 1.0d-6
      itask = 1
      istate = 1
      iopt = 0
      lrw = 44
      liw = 21
      jt = 2
      ng = 2
      write (lout,110) itol,rtol,atol(1),jt
 110  format(1h1/1x,41h demonstration program for lsodar package/////
     1  1x,14h first problem///
     2  1x,54h problem is  dy/dt = ((2*log(y)+8)/t - 5)*y,  y(1) = 1//
     3  1x,41h solution is  y(t) = exp(-t**2 + 5*t - 4)//
     4  1x,21h root functions are../
     5  10x,30h g1 = dy/dt  (root at t = 2.5)/
     6  10x,55h g2 = log(y) - 6.2491  (roots at t = 2.47 and t = 2.53)//
     7  1x,7h itol =,i3,9h   rtol =,e10.1,9h   atol =,e10.1/
     8  1x,5h jt =,i3/////)
c
c call lsodar in loop over tout values 2,3,4,5,6.
      ero = 0.0d0
      do 180 iout = 1,5
 120    continue
        call lsodar(f1,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jdum,jt,gr1,ng,jroot)
c
c print y and error in y, and print warning if error too large.
        yt = dexp(-t*t + 5.0d0*t - 4.0d0)
        er = y(1) - yt
        write (lout,130) t,y(1),er
 130    format(1x,7h at t =,e15.7,5x,3hy =,e15.7,5x,7herror =,e12.4)
        if (istate .lt. 0) go to 185
        er = dabs(er)/atol(1)
        ero = dmax1(ero,er)
        if (er .lt. 1000.0d0) go to 140
        write (lout,135)
 135    format(//1x,41h warning.. error exceeds 1000 * tolerance//)
        nerr = nerr + 1
 140    continue
        if (istate .ne. 3) go to 175
c
c if a root was found, write results and check root location.
c then reset istate to 2 and return to lsodar call.
        write (lout,150) t,jroot(1),jroot(2)
 150    format(/1x,18h root found at t =,e15.7,5x,7hjroot =,2i5)
        if (jroot(1) .eq. 1) errt = t - 2.5d0
        if (jroot(2) .eq. 1 .and. t .le. 2.5d0) errt = t - 2.47d0
        if (jroot(2) .eq. 1 .and. t .gt. 2.5d0) errt = t - 2.53d0
        write (lout,160) errt
 160    format(1x,31h error in t location of root is,e12.4//)
        if (dabs(errt) .lt. 1.0d-3) go to 170
        write (lout,165)
 165    format(//1x,36h warning.. root error exceeds 1.0e-3//)
        nerr = nerr + 1
 170    istate = 2
        go to 120
c
c if no root found, increment tout and loop back.
 175    tout = tout + 1.0d0
 180    continue
c
c problem complete.  print final statistics.
 185  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nge = iwork(10)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,190) lenrw,leniw,nst,nfe,nfea,nje,nge,ero
 190  format(//1x,32h final statistics for this run../
     1  1x,13h rwork size =,i4,15h   iwork size =,i4/
     2  1x,18h number of steps =,i5/
     3  1x,18h number of f-s   =,i5/
     4  1x,18h (excluding j-s) =,i5/
     5  1x,18h number of j-s   =,i5/
     6  1x,18h number of g-s   =,i5/
     7  1x,16h error overrun =,e10.2)
c
c-----------------------------------------------------------------------
c second problem (van der pol oscillator).
c the initial value problem is..
c   dy1/dt = y2,  dy2/dt = 100*(1 - y1**2)*y2 - y1,
c   y1(0) = 2,  y2(0) = 0,  0 .le. t .le. 200
c the root function is  g = y1.
c an analytic solution is not known, but the zeros of y1 are known
c to 15 figures for purposes of checking the accuracy.
c-----------------------------------------------------------------------
c set tolerance parameters and print heading.
      itol = 2
      rtol = 1.0d-6
      atol(1) = 1.0d-6
      atol(2) = 1.0d-4
      write (lout,200) itol,rtol,atol(1),atol(2)
 200  format(1h1/1x,40h second problem (van der pol oscillator)///
     1  1x,56h problem is dy1/dt = y2,  dy2/dt = 100*(1-y1**2)*y2 - y1/
     2  1x,33h            y1(0) = 2,  y2(0) = 0//
     3  1x,25h root function is  g = y1//
     4  1x,7h itol =,i3,9h   rtol =,e10.1,9h   atol =,2e10.1///)
c
c loop over jt = 1, 2.  set remaining parameters and print jt.
      do 290 jt = 1,2
      neq = 2
      y(1) = 2.0d0
      y(2) = 0.0d0
      t = 0.0d0
      tout = 20.0d0
      itask = 1
      istate = 1
      iopt = 0
      lrw = 57
      liw = 22
      ng = 1
      write (lout,210) jt
 210  format(////////1x,19h solution with jt =,i2////)
c
c call lsodar in loop over tout values 20,40,...,200.
      do 270 iout = 1,10
 220    continue
        call lsodar(f2,neq,y,t,tout,itol,rtol,atol,itask,istate,
     1     iopt,rwork,lrw,iwork,liw,jac2,jt,gr2,ng,jroot)
c
c print y1 and y2.
        write (lout,230) t,y(1),y(2)
 230    format(1x,7h at t =,e15.7,5x,4hy1 =,e15.7,5x,4hy2 =,e15.7)
        if (istate .lt. 0) go to 275
        if (istate .ne. 3) go to 265
c
c if a root was found, write results and check root location.
c then reset istate to 2 and return to lsodar call.
        write (lout,240) t
 240    format(/1x,18h root found at t =,e15.7)
        kroot = ifix(sngl(t/81.2d0 + 0.5d0))
        tzero = 81.17237787055d0 + dfloat(kroot-1)*81.41853556212d0
        errt = t - tzero
        write (lout,250) errt
 250    format(1x,31h error in t location of root is,e12.4//)
        if (errt .lt. 1.0d0) go to 260
        write (lout,255)
 255    format(//1x,33h warning.. root error exceeds 1.0//)
        nerr = nerr + 1
 260    istate = 2
        go to 220
c
c if no root found, increment tout and loop back.
 265    tout = tout + 20.0d0
 270    continue
c
c problem complete.  print final statistics.
 275  continue
      if (istate .lt. 0) nerr = nerr + 1
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nge = iwork(10)
      lenrw = iwork(17)
      leniw = iwork(18)
      nfea = nfe
      if (jt .eq. 2) nfea = nfe - neq*nje
      write (lout,280) lenrw,leniw,nst,nfe,nfea,nje,nge
 280  format(//1x,32h final statistics for this run../
     1  1x,13h rwork size =,i4,15h   iwork size =,i4/
     2  1x,18h number of steps =,i5/
     3  1x,18h number of f-s   =,i5/
     4  1x,18h (excluding j-s) =,i5/
     5  1x,18h number of j-s   =,i5/
     6  1x,18h number of g-s   =,i5)
 290  continue
c
c
      write (lout,300) nerr
 300  format(////1x,31h number of errors encountered =,i3)
      stop
      end
      subroutine f1 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(1), ydot(1)
      ydot(1) = ((2.0d0*dlog(y(1)) + 8.0d0)/t - 5.0d0)*y(1)
      return
      end
      subroutine gr1 (neq, t, y, ng, groot)
      integer neq, ng
      double precision t, y, groot
      dimension y(1), groot(2)
      groot(1) = ((2.0d0*dlog(y(1)) + 8.0d0)/t - 5.0d0)*y(1)
      groot(2) = dlog(y(1)) - 2.2491d0
      return
      end
      subroutine f2 (neq, t, y, ydot)
      integer neq
      double precision t, y, ydot
      dimension y(2), ydot(2)
      ydot(1) = y(2)
      ydot(2) = 100.0d0*(1.0d0 - y(1)*y(1))*y(2) - y(1)
      return
      end
      subroutine jac2 (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(2), pd(nrowpd,2)
      pd(1,1) = 0.0d0
      pd(1,2) = 1.0d0
      pd(2,1) = -200.0d0*y(1)*y(2) - 1.0d0
      pd(2,2) = 100.0d0*(1.0d0 - y(1)*y(1))
      return
      end
      subroutine gr2 (neq, t, y, ng, groot)
      integer neq, ng
      double precision t, y, groot
      dimension y(2), groot(1)
      groot(1) = y(1)
      return
      end
c-----------------------------------------------------------------------
c demonstration program for the lsodi package.
c this is the version of 18 august, 1981.
c
c this version is in double precision.
c
c for computer systems requiring a program card, the following (with
c the c in column 1 removed) may be used..
c     program lsidem(lsiout,tape6=lsiout)
c
c this program solves a semi-discretized form of the burgers equation,
c
c     u  = -(u*u/2)  + eta * u
c      t           x          xx
c
c for a = -1 .le. x .le. 1 = b, t .ge. 0.
c here eta = 0.05.
c boundary conditions.. u(-1,t) = u(1,t) = 0.
c initial profile.. square wave
c     u(0,x) = 0    for 1/2 .lt. abs(x) .le. 1
c     u(0,x) = 1/2  for abs(x) = 1/2
c     u(0,x) = 1    for 0 .le. abs(x) .lt. 1/2
c
c an ode system is generated by a simplified galerkin treatment
c of the spatial variable x.
c
c reference..
c r. c. y. chin, g. w. hedstrom, and k. e. karlsson,
c a simplified galerkin method for hyperbolic equations,
c math. comp., vol. 33, no. 146 (april 1979), pp. 647-658.
c
c the problem is run with the lsodi package with a 10-point mesh
c and a 100-point mesh.  in each case, it is run with two tolerances
c and for various appropriate values of the method flag mf.
c output is on unit lout, set to 6 in a data statement below.
c-----------------------------------------------------------------------
      external res, addabd, addafl, jacbd, jacfl
      integer i, io, istate, itol, iwork, j, k,
     1   l, lout, liw, lrw, meth, miter, mf, ml, mu,
     2   n, nout, npts, nerr,
     3   nptsm1, n14, n34, n14m1, n14p1, n34m1, n34p1
      integer nm1
      double precision a, b, eta, delta,
     1   zero, fourth, half, one, hun,
     2   t, tout, tlast, tinit, errfac,
     3   atol, rtol, rwork, y, ydoti, elkup
      double precision eodsq, r4d
      dimension y(99), ydoti(99), tout(4), atol(2), rtol(2)
      dimension rwork(2002), iwork(125)
c pass problem parameters in the common block test1.
      common /test1/ r4d, eodsq, nm1
c
c set problem parameters and run parameters
      data eta/0.05d0/, a/-1.0d0/, b/1.0d0/
      data zero/0.0d0/, fourth/0.25d0/, half/.5d0/, one/1.0d0/,
     1   hun/100.0d0/
      data tinit/0.0d0/, tlast/0.4d0/
      data tout/.10d0,.20d0,.30d0,.40d0/
      data ml/1/, mu/1/, lout/6/
      data nout/4/, lrw/2002/, liw/125/
      data itol/1/, rtol/1.0d-3, 1.0d-6/, atol/1.0d-3, 1.0d-6/
c
      iwork(1) = ml
      iwork(2) = mu
      nerr = 0
c
c loop over two values of npts.
      do 300  npts = 10, 100, 90
c
c compute the mesh width delta and other parameters.
      delta = (b - a)/dfloat(npts)
      r4d = fourth/delta
      eodsq = eta/delta**2
      nptsm1 = npts - 1
      n14 = npts/4
      n34 = 3 * n14
      n14m1 = n14 - 1
      n14p1 = n14m1 + 2
      n34m1 = n34 - 1
      n34p1 = n34m1 + 2
      n = nptsm1
      nm1 = n - 1
c
c set the initial profile (for output purposes only).
c
      do 10 i = 1,n14m1
   10   y(i) = zero
      y(n14) = half
      do 20 i = n14p1,n34m1
   20   y(i) = one
      y(n34) = half
      do 30 i = n34p1,nptsm1
   30   y(i) = zero
c
      write (lout,1000)
      write (lout,1100) eta,a,b,tinit,tlast,ml,mu,n
      write (lout,1200) zero, (y(i), i=1,n), zero
c
c the j loop is over error tolerances.
c
      do 200 j = 1,2
c
c   this method flag loop is for demonstration only.
c
      do 100 meth = 1,2
       do 100 miter = 1,5
        if ( miter .eq. 3 )  go to 100
        if ( miter .le. 2 .and. npts .gt. 10 )  go to 100
        if ( miter .eq. 5 .and. npts .lt. 100 )  go to 100
        mf = 10*meth + miter
c
c set the initial profile.
c
        do 40 i = 1,n14m1
   40     y(i) = zero
        y(n14) = half
        do 50 i = n14p1,n34m1
   50     y(i) = one
        y(n34) = half
        do 60 i = n34p1,nptsm1
   60     y(i) = zero
c
        t = tinit
        istate = 0
c
        write (lout,1500) itol, rtol(j), atol(j), mf, npts
c
c  output loop for each case
c
        do 80 io = 1,nout
c
c         call lsodi
          if ( miter .le. 2 )  call lsodi ( res, addafl, jacfl,
     1      n, y, ydoti, t, tout(io), itol, rtol(j), atol(j), 1,
     2      istate, 0, rwork, lrw, iwork, liw, mf )
          if ( miter .ge. 4 )  call lsodi ( res, addabd, jacbd,
     1      n, y, ydoti, t, tout(io), itol, rtol(j), atol(j), 1,
     2      istate, 0, rwork, lrw, iwork, liw, mf )
          write (lout,2000) t, rwork(11), iwork(14),zero,
     1                      (y(i), i=1,n), zero
c
c if istate is not 2 on return, go to error section.
          if (istate .ne. 2)  go to 90
c
   80     continue
c
c
        write (lout,3000) mf, iwork(11), iwork(12), iwork(13),
     1                iwork(17), iwork(18)
c
c estimate final error and print result.
        errfac = elkup( n, y, rwork(21), itol, rtol(j), atol(j) )
        if ( errfac .gt. hun )  go to 85
          write (lout,5000)  errfac
          go to 100
   85     write (lout,5001)  errfac
          nerr = nerr + 1
          go to 100
   90   write (lout,4000) mf, t, istate
        nerr = nerr + 1
  100   continue
  200 continue
  300 continue
c
      write (lout,6000) nerr
      stop
c
 1000 format(1h1///20x,32h demonstration problem for lsodi )
 1100 format(//10x,35h-- simplified galerkin solution of ,
     1       19hburgers equation --///
     1       13x,30hdiffusion coefficient is eta =,e10.2/
     1       13x,24huniform mesh on interval,e12.3,4h to ,e12.3/
     2       13x,24hzero boundary conditions/
     2       13x,30hinitial data were as follows..//20x,5ht0 = ,e12.5/
     2       20x,8htlast = ,e12.5/20x,5hml = ,i2/
     4       20x,5hmu = ,i2/20x,5hn  = ,i3//)
c
 1200 format(//1x,18h initial profile..,89x,e12.4/10(10e12.4/))
c
 1500 format(1h1/1x,15hrun with itol =,i2,8h  rtol =,e12.2,
     1       8h  atol =,e12.2,7h   mf =,i3,9h   npts =,i4///)
c
 2000 format(/1x,20houtput for time t = ,e12.5,16h     current h =,
     1       e12.5,18h   current order =,i2,37x,e12.5/10(10e12.4/))
c
 3000 format(1x,80(1h*)/1x,26hfinal statistics for mf = ,i2,3h.. ,
     1       i5,7h steps,,i6,5h res,,i6,11h jacobians,,
     2       16h    rwork size =,i6,17h,    iwork size =,i6)
c
 4000 format(/1x,80(1h*)//20x,28hfinal time reached for mf = ,i2,
     1       9h was t = ,e12.5/25x,18hat which istate = ,i2//1x,80(1h*))
 5000 format(1x,34hfinal output is correct to within ,e8.1,
     1       30h  times local error tolerance. )
 5001 format(1x,25hfinal output is wrong by ,e8.1,
     1       30h  times local error tolerance. )
 6000 format(////1x,47h run completed.. number of errors encountered =,
     1       i3)
c
c end of main program for the lsodi demonstration problem.
      end
      subroutine gfun (n, t, y, g)
c this subroutine computes the function g(y,t) for the
c lsodi demonstration problem.
c it uses r4d = 1/(4*delta), eodsq = eta/delta**2, and nm1 = n - 1
c from the common block test1.
c
      integer i, n, nm1
      double precision t, y, g, r4d, eodsq, two
      dimension g(n), y(n)
      common /test1/ r4d, eodsq, nm1
      data two/2.0d0/
c
      g(1) = -r4d*y(2)**2 + eodsq*(y(2) - two*y(1))
c
      do 20 i = 2,nm1
        g(i) = r4d*(y(i-1)**2 - y(i+1)**2)
     1        + eodsq*(y(i+1) - two*y(i) + y(i-1))
   20   continue
c
      g(n) = r4d*y(nm1)**2 + eodsq*(y(nm1) - two*y(n))
c
      return
c end of subroutine gfun for the lsodi demonstration problem.
      end
      subroutine addabd (n, t, y, ml, mu, pa, m0)
c this subroutine computes the matrix a in band form, adds it to pa,
c and returns the sum in pa.
c this coding is for the lsodi demonstration problem.
c the matrix a is tridiagonal, of order n, with nonzero elements
c (reading across) of  1/6, 4/6, 1/6.
c
      integer i, n, m0, ml, mu, mup1, mup2
      double precision t, y, pa, fact1, fact4, one, four, six
      dimension y(n), pa(m0,n)
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
c set the pointers.
      mup1 = mu + 1
      mup2 = mu + 2
c compute the elements of a.
      fact1 = one/six
      fact4 = four/six
c add the matrix a to the matrix pa (banded).
      do 10 i = 1,n
        pa(mu,i) = pa(mu,i) + fact1
        pa(mup1,i) = pa(mup1,i) + fact4
        pa(mup2,i) = pa(mup2,i) + fact1
   10   continue
      return
c end of subroutine addabd for the lsodi demonstration problem.
      end
      subroutine addafl (n, t, y, ml, mu, pa, m0)
c this subroutine computes the matrix a in full form, adds it to
c pa, and returns the sum in pa.
c it uses nm1 = n - 1 from common.
c this coding is for the lsodi demonstration problem.
c the matrix a is tridiagonal, of order n, with nonzero elements
c (reading across) of  1/6, 4/6, 1/6.
c
      integer i, n, m0, ml, mu, nm1
      double precision t, y, pa, r4d, eodsq, one, four, six,
     1   fact1, fact4
      dimension y(n), pa(m0,n)
      common /test1/ r4d, eodsq, nm1
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
c compute the elements of a.
      fact1 = one/six
      fact4 = four/six
c
c add the matrix a to the matrix pa (full).
c
      do 110  i = 2, nm1
         pa(i,i+1) = pa(i,i+1) + fact1
         pa(i,i) = pa(i,i) + fact4
         pa(i,i-1) = pa(i,i-1) + fact1
  110    continue
      pa(1,2) = pa(1,2) + fact1
      pa(1,1) = pa(1,1) + fact4
      pa(n,n) = pa(n,n) + fact4
      pa(n,nm1) = pa(n,nm1) + fact1
      return
c end of subroutine addafl for the lsodi demonstration problem.
      end
      subroutine jacbd (n, t, y, s, ml, mu, pa, m0)
c this subroutine computes the jacobian dg/dy = d(g-a*s)/dy
c and stores
c
c   i   j
c dg /dy   in  pa(i-j+mu+1,j)
c
c for the lsodi demonstration problem (band matrix case).
c it uses r4d = 1/(4*delta), eodsq = eta/delta**2, and nm1 = n - 1
c from the common block test1.
c
      integer i, j, n, m0, ml, mu, mup1, mup2, nm1
      double precision t, y, s, pa, diag, r4d, eodsq, two, r2d
      dimension y(n), s(n), pa(m0,n)
      common /test1/ r4d, eodsq, nm1
      data two/2.0d0/
c
      mup1 = mu + 1
      mup2 = mu + 2
      diag = -two*eodsq
      r2d = two*r4d
c                     1   1
c compute and store dg /dy
      pa(mup1,1) = diag
c
c                     1   2
c compute and store dg /dy
      pa(mu,2) = -r2d*y(2) + eodsq
c
      do 20 i = 2,nm1
c
c                     i   i-1
c compute and store dg /dy
        pa(mup2,i-1) = r2d*y(i-1) + eodsq
c
c                     i   i
c compute and store dg /dy
      pa(mup1,i) = diag
c
c                     i   i+1
c compute and store dg /dy
        pa(mu,i+1) = -r2d*y(i+1) + eodsq
   20   continue
c
c                     n   n-1
c compute and store dg /dy
      pa(mup2,nm1) = r2d*y(nm1) + eodsq
c
c                     n   n
c compute and store dg /dy
      pa(mup1,n) = diag
c
      return
c end of subroutine jacbd for the lsodi demonstration problem.
      end
      subroutine jacfl (n, t, y, s, ml, mu, pa, m0)
c this subroutine computes the jacobian dg/dy = d(g-a*s)/dy
c and stores
c
c   i   j
c dg /dy   in  pa(i,j)
c
c for the lsodi demonstration problem (full matrix case).
c it uses r4d = 1/(4*delta), eodsq = eta/delta**2, and nm1 = n - 1
c from the common block test1.
c
      integer i, j, n, m0, ml, mu, nm1
      double precision t, y, s, pa, diag, r4d, eodsq, two, r2d
      dimension y(n), s(n), pa(m0,n)
      common /test1/ r4d, eodsq, nm1
      data two/2.0d0/
c
      diag = -two*eodsq
      r2d = two*r4d
c
c                     1   1
c compute and store dg /dy
      pa(1,1) = diag
c
c                     1   2
c compute and store dg /dy
      pa(1,2) = -r2d*y(2) + eodsq
c
      do 120  i = 2,nm1
c
c                     i   i-1
c compute and store dg /dy
        pa(i,i-1) = r2d*y(i-1) + eodsq
c
c                     i   i
c compute and store dg /dy
      pa(i,i) = diag
c
c                     i   i+1
c compute and store dg /dy
        pa(i,i+1) = -r2d*y(i+1) + eodsq
  120   continue
c
c                     n   n-1
c compute and store dg /dy
      pa(n,nm1) = r2d*y(nm1) + eodsq
c
c                     n   n
c compute and store dg /dy
      pa(n,n) = diag
c
      return
c end of subroutine jacfl for the lsodi demonstration problem.
      end
      subroutine res (n, t, y, v, r, ires)
c this subroutine computes the residual vector
c   r = g(t,y) - a(t,y)*v
c for the lsodi demonstration problem.
c it uses nm1 = n - 1 from common.
c if ires = -1, only g(t,y) is returned in r, since a(t,y) does
c not depend on y.
c no changes need to be made to this routine if n is changed.
c
      integer i, ires, n, nm1
      double precision t, y, v, r, r4d, eodsq, one, four, six,
     1   fact1, fact4
      dimension y(n), v(n), r(n)
      common /test1/ r4d, eodsq, nm1
      data one /1.0d0/, four /4.0d0/, six /6.0d0/
c
      call gfun (n, t, y, r)
      if (ires .eq. -1) return
c
      fact1 = one/six
      fact4 = four/six
      r(1) = r(1) - (fact4*v(1) + fact1*v(2))
      do 10 i = 2, nm1
  10   r(i) = r(i) - (fact1*v(i-1) + fact4*v(i) + fact1*v(i+1))
      r(n) = r(n) - (fact1*v(nm1) + fact4*v(n))
      return
c end of subroutine res for the lsodi demonstration problem.
      end
      double precision function elkup (n, y, ewt, itol, rtol, atol)
c elkup looks up approximately correct values of y at t = 0.4,
c ytrue = y9 or y99 depending on whether n = 9 or 99.  these were
c obtained by running lsodi with very tight tolerances.
c the returned value is
c elkup  =  norm of  ( y - ytrue ) / ( rtol*abs(ytrue) + atol ).
c
      integer n, itol, i
      double precision y, ewt, rtol, atol, y9, y99, y99a, y99b, y99c,
     1   y99d, y99e, y99f, y99g, vnorm
      dimension y(n), ewt(n), y9(9), y99(99)
      dimension y99a(16), y99b(16), y99c(16), y99d(16), y99e(16),
     1          y99f(16), y99g(3)
      equivalence (y99a(1),y99(1)), (y99b(1),y99(17)),
     1      (y99c(1),y99(33)), (y99d(1),y99(49)), (y99e(1),y99(65)),
     1      (y99f(1),y99(81)), (y99g(1),y99(97))
      data y9 /
     1 1.07001457d-01, 2.77432492d-01, 5.02444616d-01, 7.21037157d-01,
     1 9.01670441d-01, 8.88832048d-01, 4.96572850d-01, 9.46924362d-02,
     1-6.90855199d-03 /
      data y99a /
     1 2.05114384d-03, 4.19527452d-03, 6.52533872d-03, 9.13412751d-03,
     1 1.21140191d-02, 1.55565301d-02, 1.95516488d-02, 2.41869487d-02,
     1 2.95465081d-02, 3.57096839d-02, 4.27498067d-02, 5.07328729d-02,
     1 5.97163151d-02, 6.97479236d-02, 8.08649804d-02, 9.30936515d-02 /
      data y99b /
     1 1.06448659d-01, 1.20933239d-01, 1.36539367d-01, 1.53248227d-01,
     1 1.71030869d-01, 1.89849031d-01, 2.09656044d-01, 2.30397804d-01,
     1 2.52013749d-01, 2.74437805d-01, 2.97599285d-01, 3.21423708d-01,
     1 3.45833531d-01, 3.70748792d-01, 3.96087655d-01, 4.21766871d-01 /
      data y99c /
     1 4.47702161d-01, 4.73808532d-01, 5.00000546d-01, 5.26192549d-01,
     1 5.52298887d-01, 5.78234121d-01, 6.03913258d-01, 6.29252015d-01,
     1 6.54167141d-01, 6.78576790d-01, 7.02400987d-01, 7.25562165d-01,
     1 7.47985803d-01, 7.69601151d-01, 7.90342031d-01, 8.10147715d-01 /
      data y99d /
     1 8.28963844d-01, 8.46743353d-01, 8.63447369d-01, 8.79046021d-01,
     1 8.93519106d-01, 9.06856541d-01, 9.19058529d-01, 9.30135374d-01,
     1 9.40106872d-01, 9.49001208d-01, 9.56853318d-01, 9.63702661d-01,
     1 9.69590361d-01, 9.74555682d-01, 9.78631814d-01, 9.81840924d-01 /
      data y99e /
     1 9.84188430d-01, 9.85656465d-01, 9.86196496d-01, 9.85721098d-01,
     1 9.84094964d-01, 9.81125395d-01, 9.76552747d-01, 9.70041743d-01,
     1 9.61175143d-01, 9.49452051d-01, 9.34294085d-01, 9.15063568d-01,
     1 8.91098383d-01, 8.61767660d-01, 8.26550038d-01, 7.85131249d-01 /
      data y99f /
     1 7.37510044d-01, 6.84092540d-01, 6.25748369d-01, 5.63802368d-01,
     1 4.99946558d-01, 4.36077986d-01, 3.74091566d-01, 3.15672765d-01,
     1 2.62134958d-01, 2.14330497d-01, 1.72640946d-01, 1.37031155d-01,
     1 1.07140815d-01, 8.23867920d-02, 6.20562432d-02, 4.53794321d-02 /
      data y99g / 3.15789227d-02, 1.98968820d-02, 9.60472135d-03 /
c
      if ( n - 9 )  9,9, 99
c
c compute local error tolerance using correct y (n=9).
c
    9 call ewset( n, itol, rtol, atol, y9, ewt )
c
c invert ewt and replace y by the error, y - ytrue.
c
      do 19  i = 1, 9
        ewt(i) = 1.0d0/ewt(i)
   19   y(i) = y(i) - y9(i)
      go to 200
c
c compute local error tolerance using correct y (n=99).
c
   99 call ewset( n, itol, rtol, atol, y99, ewt )
c
c invert ewt and replace y by the error, y - ytrue.
c
      do 199  i = 1, 99
        ewt(i) = 1.0d0/ewt(i)
  199   y(i) = y(i) - y99(i)
c
c find weighted norm of the error.
c
  200 elkup = vnorm (n, y, ewt)
      return
c end of function elkup for the lsodi demonstration program.
      end
c-----------------------------------------------------------------------
c demonstration program for the lsoibt package.
c this the version of june 22, 1984.
c this version is in double precision.
c
c this program solves a semi-discretized form of the following system
c of three pdes (each similar to a burgers equation)..
c
c   u(i)   =  -(u(1)+u(2)+u(3)) u(i)   +  eta(i) u(i)    (i=1,2,3),
c       t                           x                xx
c
c on the interval  -1 .le. x .le. 1, and with time t .ge. 0.
c the diffusion coefficients are eta(*) = .1, .02, .01.
c the boundary conditions are u(i) = 0 at x = -1 and x = 1 for all i.
c the initial profile for each u(i) is a square wave..
c     u(i) = 0         on 1/2 .lt. abs(x) .le. 1
c     u(i) = amp(i)/2  on abs(x) = 1/2
c     u(i) = amp(i)    on 0 .le. abs(x) .lt. 1/2
c where the amplitudes are amp(*) = .2, .3, .5.
c
c a simplified galerkin treatment of the spatial variable x is used,
c with piecewise linear basis functions on a uniform mesh of 100
c intervals.  the result is a system of odes in the discrete values
c u(i,k) approximating u(i)  (i=1,2,3) at the interior points
c (k = 1,...,99).  the ode-s are..
c
c    .            .        .
c   (u(i,k-1) + 4 u(i,k) + u(i,k+1))/6  =
c
c     -(1/6dx) (c(k-1)dul(i) + 2c(k)(dul(i)+dur(i)) + c(k+1)dur(i))
c
c     + (eta(i)/dx**2) (dur(i) - dul(i))     (i=1,2,3,  k=1,...,99),
c
c where
c     c(j) = u(1,j)+u(2,j)+u(3,j),   dx = .02 = the interval size,
c     dul(i) = u(i,k) - u(i,k-1),   dur(i) = u(i,k+1) - u(i,k).
c terms involving boundary values (subscripts 0 or 100) are dropped
c from the equations for k = 1 and k = 99 above.
c
c the problem is run for each of the 4 values of mf, and for two values
c of the tolerances.  output is taken at t = .1, .2, .3, .4.
c output is on unit lout, set to 6 in a data statement below.
c-----------------------------------------------------------------------
      external res, addabt, jacbt
      integer ncomp, nip, nm1
      integer i, io, istate, itol, iwork, jtol, lout, liw, lrw,
     1   meth, miter, mf, neq, nerr, nint, nout
      double precision eodsq, r6d
      double precision abermx, atol, dx, errfac, eta, hun, one,
     1   rtol, rwork, six, t, tinit, tlast, tout, tols, two, y, ydoti
      dimension eta(3), y(297), ydoti(297), tout(4), tols(2)
      dimension rwork(7447), iwork(317)
c pass problem parameters in the common block par.
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
c
c set problem parameters and run parameters
      data eta/0.1d0,0.02d0,0.01d0/, tinit/0.0d0/, tlast/0.4d0/
      data one/1.0d0/, two/2.0d0/, six/6.0d0/, hun/100.0d0/
      data tout/.10d0,.20d0,.30d0,.40d0/
      data lout/6/, nout/4/, lrw/7447/, liw/317/
      data itol/1/, tols/1.0d-3, 1.0d-6/
c
c set mesh parameters nint, dxc etc.
      nint = 100
      ncomp = 3
      dx = two/dfloat(nint)
      r6d = one/(six*dx)
      do 10 i = 1,ncomp
 10     eodsq(i) = eta(i)/dx**2
      nip = nint - 1
      neq = ncomp*nip
      nm1 = nip - 1
      iwork(1) = ncomp
      iwork(2) = nip
c
      nerr = 0
c
c set the initial conditions (for output purposes only).
      call setic (nint, ncomp, y)
c
      write (lout,1000)
      write (lout,1100) (eta(i),i=1,ncomp), tinit, tlast, nint,
     1   ncomp, nip, neq
      write (lout,1200)
      call edit (y, ncomp, nip, lout)
c
c the jtol loop is over error tolerances.
      do 200 jtol = 1,2
      rtol = tols(jtol)
      atol = rtol
c
c the meth/miter loops cover 4 values of method flag mf.
      do 100 meth = 1,2
       do 100 miter = 1,2
        mf = 10*meth + miter
c
c set the initial conditions.
        call setic (nint, ncomp, y)
        t = tinit
        istate = 0
c
        write (lout,1500)  rtol, atol, mf
c
c loop over output times for each case
        do 80 io = 1,nout
c
          call lsoibt (res, addabt,jacbt, neq, y, ydoti, t, tout(io),
     1     itol,rtol,atol, 1, istate, 0, rwork,lrw,iwork,liw, mf)
c
          write (lout,2000) t, rwork(11), iwork(14), iwork(11)
          if (io .eq. nout) call edit (y, ncomp, nip, lout)
c
c if istate is not 2 on return, go to error section.
          if (istate .ne. 2)  go to 90
c
 80       continue
c
c print final statistics.
        write (lout,3000) mf, iwork(11), iwork(12), iwork(13),
     1                iwork(17), iwork(18)
c
c estimate final error and print result.
        call maxerr (y, ncomp, nip, abermx)
        errfac = abermx/tols(jtol)
        if (errfac .gt. hun)  go to 85
          write (lout,5000) errfac
          go to 100
 85       write (lout,5100) errfac
          nerr = nerr + 1
          go to 100
 90     write (lout,4000) mf, t, istate
        nerr = nerr + 1
 100    continue
 200  continue
c
      write (lout,6000) nerr
      stop
c
 1000 format(1h1/20x,33h demonstration problem for lsoibt//
     1   10x,47hgalerkin method solution of system of 3 pde-s..//
     2   10x,51h  u(i)   =  -(u(1)+u(2)+u(3)) u(i)   +  eta(i) u(i),
     3   5x,9h(i=1,2,3)/16x,1ht,27x,1hx,16x,2hxx//
     4   10x,48hx interval is -1 to 1,  zero boundary conditions/
     5   10x,52hx discretized using piecewise linear basis functions)
 1100 format(/10x,33hfixed parameters are as follows../
     1       13x,32hdiffusion coefficients are eta =,3e10.2/
     2       13x,5ht0 = ,e12.5/13x,8htlast = ,e12.5/
     3       13x,35huniform mesh, number of intervals =,i4/
     4       13x,15hblock size mb =,i2/13x,15hno. blocks nb =,i4/
     5       13x,21hode system size neq =,i5//)
c
 1200 format(/1x,19h initial profiles../)
c
 1500 format(1h1/1x,15hrun with rtol =,e9.1,8h  atol =,e9.1,
     1       7h   mf =,i3///)
c
 2000 format(1x,12hat time t = ,e12.5,16h     current h =,e12.5,
     1       18h   current order =,i2,16h   current nst =,i5/)
c
 3000 format(//1x,26hfinal statistics for mf = ,i2,3h.. ,
     1       i5,7h steps,,i6,5h res,,i6,11h jacobians,,
     2       16h    rwork size =,i6,17h,    iwork size =,i6)
c
 4000 format(//1x,20x,28hfinal time reached for mf = ,i2,
     1       9h was t = ,e12.5/25x,18hat which istate = ,i2//1x,80(1h*))
 5000 format(1x,34hfinal output is correct to within ,e8.1,
     1       30h  times local error tolerance. )
 5100 format(1x,25hfinal output is wrong by ,e8.1,
     1       30h  times local error tolerance. )
 6000 format(////1x,17h run completed.. ,i3,19h errors encountered)
c
c end of main program for the lsoibt demonstration problem.
      end
      subroutine setic (nint, mb, y)
c this routine loads the y array with initial data based on a
c square wave profile for each of the mb pde variables.
c
      integer nint, mb, i, k, nip, n14, n14m1, n14p1, n34, n34m1, n34p1
      double precision y,  amp, half, zero
      dimension y(mb,1), amp(3)
      data zero/0.0d0/, half/0.5d0/, amp/0.2d0,0.3d0,0.5d0/
c
      nip = nint - 1
      n14 = nint/4
      n34 = 3*n14
      n14m1 = n14 - 1
      n14p1 = n14 + 1
      n34m1 = n34 - 1
      n34p1 = n34 + 1
c
      do 15 k = 1,n14m1
        do 10 i = 1,mb
 10       y(i,k) = zero
 15     continue
c
      do 20 i = 1,mb
 20     y(i,n14) = half*amp(i)
c
      do 35 k = n14p1,n34m1
        do 30 i = 1,mb
 30       y(i,k) = amp(i)
 35     continue
c
      do 40 i = 1,mb
 40     y(i,n34) = half*amp(i)
c
      do 55 k = n34p1,nip
        do 50 i = 1,mb
 50       y(i,k) = zero
 55     continue
c
      return
c end of subroutine setic
      end
      subroutine res (n, t, y, v, r, ires)
c this subroutine computes the residual vector
c   r = g(t,y) - a(t,y)*v
c for the lsoibt demonstration problem, using routines gfun and subav.
c if ires = -1, only g(t,y) is returned in r, since a(t,y) does
c not depend on y.
c no changes need to be made to this routine if nip is changed.
c
      integer ires, n,  ncomp, nip, nm1
      double precision t, y, v, r, r6d, eodsq
      dimension y(n), v(n), r(n)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
c
      call gfun (t, y, r, ncomp)
      if (ires .eq. -1) return
c
      call subav (r, v, ncomp)
c
      return
c end of subroutine res
      end
      subroutine gfun (t, y, g, mb)
c this subroutine computes the function g(y,t) for the
c lsoibt demonstration problem.
c it uses r6d = 1/(6*dx), eodsq(*) = eta(*)/dx**2, nip,
c and nm1 = nip - 1 from the common block par.
c
      integer mb,  ncomp, nip, nm1,  i, k
      double precision t, y, g,  r6d, eodsq,  cc, cl, cr, dli, dri, two
      dimension g(mb,1), y(mb,1)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data two/2.0d0/
c
c left-most interior point (k = 1)
      cc = y(1,1) + y(2,1) + y(3,1)
      cr = y(1,2) + y(2,2) + y(3,2)
      do 10 i = 1,mb
        dri = y(i,2) - y(i,1)
        g(i,1) = -r6d*(two*cc*y(i,2) + cr*dri)
     1         + eodsq(i)*(dri - y(i,1))
 10     continue
c
c interior points k = 2 to nip-1
      do 20 k = 2,nm1
        cl = y(1,k-1) + y(2,k-1) + y(3,k-1)
        cc = y(1,k) + y(2,k) + y(3,k)
        cr = y(1,k+1) + y(2,k+1) + y(3,k+1)
        do 15 i = 1,mb
          dli = y(i,k) - y(i,k-1)
          dri = y(i,k+1) - y(i,k)
          g(i,k) = -r6d*(cl*dli + two*cc*(dli + dri) + cr*dri)
     1           + eodsq(i)*(dri - dli)
 15       continue
 20     continue
c
c right-most interior point (k = nip)
      cl = y(1,nm1) + y(2,nm1) + y(3,nm1)
      cc = y(1,nip) + y(2,nip) + y(3,nip)
      do 30 i = 1,mb
        dli = y(i,nip) - y(i,nm1)
        g(i,nip) = -r6d*(cl*dli - two*cc*y(i,nm1))
     1           - eodsq(i)*(y(i,nip) + dli)
 30     continue
c
      return
c end of subroutine gfun
      end
      subroutine subav (r, v, mb)
c this routine subtracts the matrix a time the vector v from r,
c in order to form the residual vector, stored in r.
c
      integer mb,  ncomp, nip, nm1,  i, k
      double precision r, v,  r6d, eodsq,  aa1, aa4, four, one, six
      dimension r(mb,1), v(mb,1)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data one /1.0d0/, four /4.0d0/, six /6.0d0/
c
      aa1 = one/six
      aa4 = four/six
c
      do 10 i = 1,mb
 10     r(i,1) = r(i,1) - (aa4*v(i,1) + aa1*v(i,2))
c
      do 20 k = 2,nm1
        do 15 i = 1,mb
 15       r(i,k) = r(i,k) - (aa1*v(i,k-1) + aa4*v(i,k) + aa1*v(i,k+1))
 20     continue
c
      do 30 i = 1,mb
 30     r(i,nip) = r(i,nip) - (aa1*v(i,nm1) + aa4*v(i,nip))
c
      return
c end of subroutine subav
      end
      subroutine addabt (n, t, y, mb, nb, pa, pb, pc)
c this subroutine computes the elements of the matrix a,
c and adds them to pa, pb, and pc in the appropriate manner.
c this coding is for the lsoibt demonstration problem.
c the matrix a is tridiagonal, of order n,with
c nonzero elements (reading across) of 1/6,4/6,1/6.
c
      integer n, mb, nb,  i, k
      double precision pa, pb, pc, t, y,  aa1, aa4, four, one, six
      dimension y(mb,nb),pa(mb,mb,nb),pb(mb,mb,nb),pc(mb,mb,nb)
      data one/1.0d0/, four/4.0d0/, six/6.0d0/
c
      aa1 = one/six
      aa4 = four/six
      do 50 k = 1,nb
        do 10 i = 1,mb
 10       pa(i,i,k) = pa(i,i,k) + aa4
        if (k .eq. nb) go to 25
        do 20 i = 1,mb
 20       pb(i,i,k) = pb(i,i,k) + aa1
 25     continue
        if (k .eq. 1) go to 50
        do 30 i = 1,mb
 30       pc(i,i,k) = pc(i,i,k) + aa1
 50     continue
c
      return
c end of subroutine addabt
      end
      subroutine jacbt (n, t, y, s, mb, nb, pa, pb, pc)
c this subroutine computes the jacobian dg/dy=d(g-a*s)/dy
c which has block tridiagonal structure. the main, upper,
c and lower diagonals are stored in pa, pb, and pc respectively.
c
      integer n, mb, nb,  ncomp, nip, nm1,  i, j, k
      double precision t, y, s, pa, pb, pc,  r6d, eodsq,  cc, cl, cr,
     1   dlj, drj, paij, pbij, pcij, terma, termb, termc, two
      dimension y(mb,nb),s(n),pa(mb,mb,nb),pb(mb,mb,nb),pc(mb,mb,nb)
      common /par/ r6d, eodsq(3), ncomp, nip, nm1
      data two/2.0d0/
c
c left-most interior point (k = 1)
      cc = y(1,1) + y(2,1) + y(3,1)
      cr = y(1,2) + y(2,2) + y(3,2)
      terma = r6d*cr
      termb = -r6d*(two*cc + cr)
      do 20 j = 1,mb
        drj = y(j,2) - y(j,1)
        paij = -r6d*two*y(j,2)
        pbij = -r6d*drj
        do 10 i = 1,mb
          pa(i,j,1) = paij
 10       pb(i,j,1) = pbij
        pa(j,j,1) = pa(j,j,1) + terma - two*eodsq(j)
        pb(j,j,1) = pb(j,j,1) + termb + eodsq(j)
 20     continue
c
c interior points k = 2 to nip-1
      do 50 k = 2,nm1
        cl = y(1,k-1) + y(2,k-1) + y(3,k-1)
        cc = y(1,k) + y(2,k) + y(3,k)
        cr = y(1,k+1) + y(2,k+1) + y(3,k+1)
        terma = r6d*(cr - cl)
        termb = -r6d*(two*cc + cr)
        termc = r6d*(two*cc + cl)
        do 40 j = 1,mb
          dlj = y(j,k) - y(j,k-1)
          drj = y(j,k+1) - y(j,k)
          paij = -r6d*two*(dlj + drj)
          pbij = -r6d*drj
          pcij = -r6d*dlj
          do 30 i = 1,mb
            pa(i,j,k) = paij
            pb(i,j,k) = pbij
 30         pc(i,j,k) = pcij
          pa(j,j,k) = pa(j,j,k) + terma - two*eodsq(j)
          pb(j,j,k) = pb(j,j,k) + termb + eodsq(j)
          pc(j,j,k) = pc(j,j,k) + termc + eodsq(j)
 40       continue
 50     continue
c
c right-most interior point (k = nip)
      cl = y(1,nm1) + y(2,nm1) + y(3,nm1)
      cc = y(1,nip) + y(2,nip) + y(3,nip)
      terma = -r6d*cl
      termc = r6d*(two*cc + cl)
      do 70 j = 1,mb
        dlj = y(j,nip) - y(j,nm1)
        paij = r6d*two*y(j,nm1)
        pcij = -r6d*dlj
        do 60 i = 1,mb
          pa(i,j,nip) = paij
 60       pc(i,j,nip) = pcij
        pa(j,j,nip) = pa(j,j,nip) + terma - two*eodsq(j)
        pc(j,j,nip) = pc(j,j,nip) + termc + eodsq(j)
 70     continue
c
      return
c end of subroutine jacbt
      end
      subroutine edit (y, mb, nip, lout)
c this routine prints output.  for each of the mb pde components,
c the values at the nip points are printed.
c all output is on unit lout.
c
      integer mb, nip, lout,  i, k
      double precision y
      dimension y(mb,nip)
c
      do 10 i = 1,mb
        write (lout,20) i, (y(i,k),k=1,nip)
 10     continue
c
 20   format(1x,28h values of pde component i =,i3/10(10e12.4/))
c
      return
c end of subroutine edit
      end
      subroutine maxerr (y, mb, nb, abermx)
c this routine computes the maximum absolute error in the
c y array, as a computed solution at t = 0.4, using data-loaded values
c for accurate answers (from running with much smaller tolerances).
c
      integer mb, nb,  k
      double precision y, abermx,  ae1, ae2, ae3, u1, u2, u3, zero,
     1   u1a, u1b, u1c, u1d, u1e, u1f, u1g,
     2   u2a, u2b, u2c, u2d, u2e, u2f, u2g,
     3   u3a, u3b, u3c, u3d, u3e, u3f, u3g
      dimension y(mb,nb), u1(99), u2(99), u3(99)
      dimension u1a(16),u1b(16),u1c(16),u1d(16),u1e(16),u1f(16),u1g(3),
     2          u2a(16),u2b(16),u2c(16),u2d(16),u2e(16),u2f(16),u2g(3),
     3          u3a(16),u3b(16),u3c(16),u3d(16),u3e(16),u3f(16),u3g(3)
      equivalence (u1a(1),u1(1)), (u1b(1),u1(17)),
     1      (u1c(1),u1(33)), (u1d(1),u1(49)), (u1e(1),u1(65)),
     1      (u1f(1),u1(81)), (u1g(1),u1(97))
      equivalence (u2a(1),u2(1)), (u2b(1),u2(17)),
     1      (u2c(1),u2(33)), (u2d(1),u2(49)), (u2e(1),u2(65)),
     1      (u2f(1),u2(81)), (u2g(1),u2(97))
      equivalence (u3a(1),u3(1)), (u3b(1),u3(17)),
     1      (u3c(1),u3(33)), (u3d(1),u3(49)), (u3e(1),u3(65)),
     1      (u3f(1),u3(81)), (u3g(1),u3(97))
c
      data u1a /
     1  1.70956682d-03, 3.43398445d-03, 5.18783349d-03, 6.98515842d-03,
     1  8.83921016d-03, 1.07622016d-02, 1.27650806d-02, 1.48573251d-02,
     1  1.70467655d-02, 1.93394396d-02, 2.17394852d-02, 2.42490773d-02,
     1  2.68684152d-02, 2.95957660d-02, 3.24275691d-02, 3.53586054d-02/
      data u1b /
     1  3.83822285d-02, 4.14906520d-02, 4.46752791d-02, 4.79270545d-02,
     1  5.12368132d-02, 5.45956048d-02, 5.79949684d-02, 6.14271460d-02,
     1  6.48852271d-02, 6.83632267d-02, 7.18561029d-02, 7.53597274d-02,
     1  7.88708192d-02, 8.23868545d-02, 8.59059616d-02, 8.94268082d-02/
      data u1c /
     1  9.29484864d-02, 9.64703968d-02, 9.99921344d-02, 1.03513375d-01,
     1  1.07033760d-01, 1.10552783d-01, 1.14069668d-01, 1.17583246d-01,
     1  1.21091827d-01, 1.24593066d-01, 1.28083828d-01, 1.31560049d-01,
     1  1.35016617d-01, 1.38447256d-01, 1.41844451d-01, 1.45199401d-01/
      data u1d /
     1  1.48502033d-01, 1.51741065d-01, 1.54904135d-01, 1.57977973d-01,
     1  1.60948623d-01, 1.63801670d-01, 1.66522463d-01, 1.69096305d-01,
     1  1.71508595d-01, 1.73744902d-01, 1.75790974d-01, 1.77632682d-01,
     1  1.79255895d-01, 1.80646319d-01, 1.81789276d-01, 1.82669470d-01/
      data u1e /
     1  1.83270725d-01, 1.83575716d-01, 1.83565712d-01, 1.83220322d-01,
     1  1.82517279d-01, 1.81432251d-01, 1.79938706d-01, 1.78007835d-01,
     1  1.75608540d-01, 1.72707519d-01, 1.69269456d-01, 1.65257378d-01,
     1  1.60633244d-01, 1.55358941d-01, 1.49398029d-01, 1.42718981d-01/
      data u1f /
     1  1.35301474d-01, 1.27148627d-01, 1.18308730d-01, 1.08905085d-01,
     1  9.91559295d-02, 8.93515884d-02, 7.97824293d-02, 7.06663514d-02,
     1  6.21244732d-02, 5.41994827d-02, 4.68848207d-02, 4.01465202d-02,
     1  3.39357642d-02, 2.81954415d-02, 2.28635569d-02, 1.78750916d-02/
      data u1g / 1.31630892d-02, 8.65933391d-03, 4.29480447d-03/
      data u2a /
     1  7.17416019d-06, 1.70782645d-05, 3.31245126d-05, 6.01588363d-05,
     1  1.05339286d-04, 1.79174771d-04, 2.96719122d-04, 4.78862606d-04,
     1  7.53598916d-04, 1.15707860d-03, 1.73420412d-03, 2.53849668d-03,
     1  3.63099110d-03, 5.07800919d-03, 6.94782549d-03, 9.30645443d-03/
      data u2b /
     1  1.22130079d-02, 1.57152366d-02, 1.98459102d-02, 2.46205841d-02,
     1  3.00370492d-02, 3.60764461d-02, 4.27057301d-02, 4.98809820d-02,
     1  5.75510102d-02, 6.56607602d-02, 7.41541974d-02, 8.29764928d-02,
     1  9.20754824d-02, 1.01402468d-01, 1.10912474d-01, 1.20564094d-01/
      data u2c /
     1  1.30319039d-01, 1.40141489d-01, 1.49997326d-01, 1.59853293d-01,
     1  1.69676126d-01, 1.79431680d-01, 1.89084097d-01, 1.98595037d-01,
     1  2.07923034d-01, 2.17023055d-01, 2.25846345d-01, 2.34340694d-01,
     1  2.42451240d-01, 2.50121934d-01, 2.57297724d-01, 2.63927433d-01/
      data u2d /
     1  2.69967170d-01, 2.75383917d-01, 2.80158840d-01, 2.84289739d-01,
     1  2.87792167d-01, 2.90698875d-01, 2.93057586d-01, 2.94927384d-01,
     1  2.96374262d-01, 2.97466488d-01, 2.98270390d-01, 2.98847025d-01,
     1  2.99249945d-01, 2.99524080d-01, 2.99705593d-01, 2.99822450d-01/
      data u2e /
     1  2.99895431d-01, 2.99939301d-01, 2.99963931d-01, 2.99975129d-01,
     1  2.99974996d-01, 2.99961526d-01, 2.99927041d-01, 2.99854809d-01,
     1  2.99712769d-01, 2.99442742d-01, 2.98942676d-01, 2.98038511d-01,
     1  2.96441259d-01, 2.93684573d-01, 2.89040478d-01, 2.81421884d-01/
      data u2f /
     1  2.69315148d-01, 2.50874185d-01, 2.24457680d-01, 1.89885662d-01,
     1  1.49894358d-01, 1.09927672d-01, 7.54041273d-02, 4.90259517d-02,
     1  3.06080023d-02, 1.85165524d-02, 1.09104125d-02, 6.27726960d-03,
     1  3.53002680d-03, 1.94049735d-03, 1.04218859d-03, 5.45964314d-04/
      data u2g / 2.77379128d-04, 1.33343739d-04, 5.32660444d-05/
      data u3a /
     1  1.86765383d-10, 1.96772458d-09, 1.19111389d-08, 5.54964761d-08,
     1  2.18340713d-07, 7.55899524d-07, 2.35604385d-06, 6.70801745d-06,
     1  1.76224112d-05, 4.30351929d-05, 9.82592148d-05, 2.10736217d-04,
     1  4.26209304d-04, 8.15657041d-04, 1.48160943d-03, 2.56186555d-03/
      data u3b /
     1  4.22851247d-03, 6.68078970d-03, 1.01317466d-02, 1.47903961d-02,
     1  2.08424987d-02, 2.84336008d-02, 3.76573037d-02, 4.85502549d-02,
     1  6.10936693d-02, 7.52198901d-02, 9.08218891d-02, 1.07763660d-01,
     1  1.25889931d-01, 1.45034247d-01, 1.65025016d-01, 1.85689556d-01/
      data u3c /
     1  2.06856371d-01, 2.28356037d-01, 2.50021072d-01, 2.71685149d-01,
     1  2.93181998d-01, 3.14344301d-01, 3.35002907d-01, 3.54986687d-01,
     1  3.74123404d-01, 3.92241969d-01, 4.09176451d-01, 4.24772089d-01,
     1  4.38893320d-01, 4.51433444d-01, 4.62324969d-01, 4.71549073d-01/
      data u3d /
     1  4.79142163d-01, 4.85197409d-01, 4.89859810d-01, 4.93314543d-01,
     1  4.95770115d-01, 4.97439231d-01, 4.98520996d-01, 4.99187563d-01,
     1  4.99576941d-01, 4.99791928d-01, 4.99903753d-01, 4.99958343d-01,
     1  4.99983239d-01, 4.99993785d-01, 4.99997902d-01, 4.99999367d-01/
      data u3e /
     1  4.99999835d-01, 4.99999965d-01, 4.99999995d-01, 5.00000000d-01,
     1  5.00000000d-01, 4.99999997d-01, 4.99999976d-01, 4.99999863d-01,
     1  4.99999315d-01, 4.99996914d-01, 4.99987300d-01, 4.99951740d-01,
     1  4.99829328d-01, 4.99435130d-01, 4.98245007d-01, 4.94883400d-01/
      data u3f /
     1  4.86081966d-01, 4.65174923d-01, 4.21856650d-01, 3.47885738d-01,
     1  2.49649938d-01, 1.51648615d-01, 7.80173239d-02, 3.47983164d-02,
     1  1.38686441d-02, 5.05765688d-03, 1.71052539d-03, 5.38966324d-04,
     1  1.57923694d-04, 4.27352191d-05, 1.05512005d-05, 2.33068621d-06/
      data u3g / 4.45404604d-07, 6.88336884d-08, 7.23875975d-09/
      data zero/0.0d0/
c
      abermx = zero
      do 10 k = 1,99
        ae1 = dabs(y(1,k) - u1(k))
        ae2 = dabs(y(2,k) - u2(k))
        ae3 = dabs(y(3,k) - u3(k))
        abermx = dmax1(abermx, ae1, ae2, ae3)
 10     continue
c
c
      return
c end of subroutine maxerr
c end of demonstration program for lsoibt.
      end
