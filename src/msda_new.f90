! --------------------------------------------------
SUBROUTINE msda(nk,nvars,sigma,delta,pf,dfmax,pmax,nlam,flmin,ulam,&
        eps,maxit,sml,verbose,nalam,theta,m,ntheta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
    INTEGER::mnl
    INTEGER::nk
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::verbose
    INTEGER::maxit
    INTEGER::m(pmax)
    INTEGER::ntheta(nlam)
    DOUBLE PRECISION::flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION::sml
    DOUBLE PRECISION::sigma(nvars,nvars)
    DOUBLE PRECISION::delta(nk,nvars)
    DOUBLE PRECISION::pf(nvars)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::theta(nk,pmax,nlam)
    DOUBLE PRECISION::alam(nlam)
    ! - - - local declarations - - -
    INTEGER::mm(nvars)
    INTEGER::k
    INTEGER::j
    INTEGER::jj
    INTEGER::l
    INTEGER::vrg
    INTEGER::ni
    INTEGER::me
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::v
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::thetanew(nk,nvars)
    DOUBLE PRECISION::thetaold(nk,nvars)
    DOUBLE PRECISION::r(nk,nvars)
    DOUBLE PRECISION::ab(nk,nvars)
    DOUBLE PRECISION::d(nk)
    DOUBLE PRECISION::theta_sum(nk)
    DOUBLE PRECISION::thetatmp(nk)
    DOUBLE PRECISION::u(nk)
    DOUBLE PRECISION::loss_diff
    DOUBLE PRECISION::penalty_diff
    DOUBLE PRECISION::dev
    DOUBLE PRECISION::dev_tmp
    DOUBLE PRECISION::tmp1_new
    DOUBLE PRECISION::tmp2_new
    DOUBLE PRECISION::dev_new
    DOUBLE PRECISION::dev1_new
    DOUBLE PRECISION::dev2_new
    DOUBLE PRECISION::dev3_new
    DOUBLE PRECISION::tmp1_old
    DOUBLE PRECISION::tmp2_old
    DOUBLE PRECISION::dev_old
    DOUBLE PRECISION::dev1_old
    DOUBLE PRECISION::dev2_old
    DOUBLE PRECISION::dev3_old
! - - - begin - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)
! - - - some initial setup - - -
    mnl = Min (mnlam, nlam)
    r = delta
    thetanew=0.0D0
    thetaold=0.0D0
    dev=0.0D0
    m=0
    mm=0
    npass=0
    ni=npass
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    DO l=1,nlam
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al=0.0D0
                DO j = 1,nvars
                    IF(pf(j)>0.0D0) THEN
                            u = delta(:,j)
                            v = sqrt(dot_product(u,u))
                            al=max(al,v/pf(j))
                    ENDIF
                END DO
                al=al*alf
            ENDIF
        ENDIF
! --------- outer loop ----------------------------
        DO
            IF(ni>0) thetaold(:,m(1:ni))=thetanew(:,m(1:ni))
! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                dev_tmp = dev
                DO k=1,nvars
                    thetatmp=thetanew(:,k)
                    u = r(:,k)/sigma(k,k) + thetatmp
                    unorm = sqrt(dot_product(u,u))
                    v = unorm-al*pf(k)/sigma(k,k)
                    IF(v > 0.0D0) THEN
                        thetanew(:,k) = v*u/unorm
                    ELSE
                        thetanew(:,k)=0.0D0
                    ENDIF
                    d=thetanew(:,k)-thetatmp
                    theta_sum=thetanew(:,k)+thetatmp
                    IF(any(d/=0.0D0)) THEN
                        dif=max(dif,maxval(abs(d)))
                        loss_diff = sum(d*(0.5*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                        penalty_diff = al*pf(k)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                        - sqrt(dot_product(thetatmp,thetatmp)))
                        dev = dev + loss_diff + penalty_diff
                        ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                        r=r-ab
                        IF(mm(k)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                                mm(k)=ni
                                m(ni)=k
                            ENDIF
                        ENDIF
                ENDDO
                IF(abs((dev-dev_tmp)/dev)<sml) EXIT
                IF(ni>pmax) EXIT
                IF(dif<eps) EXIT
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
               ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    dev_tmp = dev
                    DO j=1,ni
                        k=m(j)
                        thetatmp=thetanew(:,k)
                        u = r(:,k)/sigma(k,k) + thetatmp
                        unorm = sqrt(dot_product(u,u))
                        v = unorm-al*pf(k)/sigma(k,k)
                        IF(v > 0.0D0) THEN
                            thetanew(:,k) = v*u/unorm
                        ELSE
                            thetanew(:,k)=0.0D0
                        ENDIF
                        d=thetanew(:,k)-thetatmp
                        theta_sum=thetanew(:,k)+thetatmp
                        IF(any(d/=0.0D0)) THEN
                            dif=max(dif,maxval(abs(d)))
                            loss_diff = sum(d*(0.5*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                            penalty_diff = al*pf(k)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                            - sqrt(dot_product(thetatmp,thetatmp)))
                            dev = dev + loss_diff + penalty_diff
                            ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                            r=r-ab
                        ENDIF
                    ENDDO
                    IF(abs((dev-dev_tmp)/dev)<sml) EXIT
                    IF(dif<eps) EXIT
                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- this is the final check ------------------------
            vrg=1
            DO j=1,ni
                IF(maxval(abs(thetanew(:,m(j))-thetaold(:,m(j))))>=eps) THEN
                    vrg=0
                    EXIT
                ENDIF
            ENDDO
            IF(vrg==1) EXIT
            ! test deviance loop
            dev1_new = 0.0
            dev2_new = 0.0
            dev1_old = 0.0
            dev2_old = 0.0
            DO jj = 1,nk
                tmp1_new = dot_product(MATMUL(thetanew(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetanew(jj,m(1:ni)))
                tmp2_new = dot_product(thetanew(jj,m(1:ni)),delta(jj,m(1:ni)))
                dev1_new = dev1_new + tmp1_new
                dev2_new = dev2_new + tmp2_new
                tmp1_old = dot_product(MATMUL(thetaold(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetaold(jj,m(1:ni)))
                tmp2_old = dot_product(thetaold(jj,m(1:ni)),delta(jj,m(1:ni)))
                dev1_old = dev1_old + tmp1_old
                dev2_old = dev2_old + tmp2_old
            ENDDO
            dev3_new = al * sum(pf(m(1:ni)) * sqrt(sum(thetanew(:,m(1:ni)) * thetanew(:,m(1:ni)), DIM = 1)))
            dev3_old = al * sum(pf(m(1:ni)) * sqrt(sum(thetaold(:,m(1:ni)) * thetaold(:,m(1:ni)), DIM = 1)))
            dev_new = 0.5 * dev1_new - dev2_new + dev3_new
            dev_old = 0.5 * dev1_old - dev2_old + dev3_old
            IF(verbose==1) THEN
                CALL intpr('Current Lambda',-1,l,1)
                CALL dblepr('Obj-func Jump',-1,abs((dev_new-dev_old)/dev_new),1)
            ENDIF
            IF(abs((dev_new-dev_old)/dev_new)<sml) EXIT
            ! test deviance loop end
        ENDDO
!--- final update variable save results------------
        IF(ni>pmax) THEN
!             jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) theta(:,1:ni,l)=thetanew(:,m(1:ni))
        me = count(maxval(abs(theta(:,1:ni,l)),dim=1)/=0.0D0)
        IF(me>dfmax) EXIT
        ntheta(l)=ni
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
    ENDDO
    RETURN
END SUBROUTINE msda
