MODULE DBEM
	!**********************************************************
	!Updated version 01/2021 with carrying capacity adjusted for NPP
	!orcinus version
	!**********************************************************
		IMPLICIT NONE
		!General variables
		INTEGER,PARAMETER	::nrow=360,ncol=720,inTimeStep=24
		INTEGER,PARAMETER	:: iIniStartYear=1971, iIniEndYear=2000,iSimStartYear=1951,iSimEndYear=2100
			INTEGER,PARAMETER       :: iIniStartYear2=1951,iIniEndYear2=1970,MPAScStartYear=1951
		INTEGER,PARAMETER	::	iLarTimeStep=7
		CHARACTER(10),PARAMETER	::	cSetfilename='settings'
		INTEGER		::	SppNo,iSpp
		CHARACTER (LEN=30)	::	cResultspath
		CHARACTER (LEN=30)	::	cTaxonDatapath
		CHARACTER(LEN=60)	::	cfolderpath
		INTEGER,ALLOCATABLE	::	iTaxonList(:)
		LOGICAL	::	lerrTaxon
		CHARACTER(20)	::	cSppListFile
	
		!INTEGER, DIMENSION(ncol,nrow,2)	:: MSeg
	
		REAL,ALLOCATABLE	::	rAbdMaster(:,:)
		!Variables for dispersal and movement
		REAL, DIMENSION(nrow*ncol)	:: rAbd,rLarv
		INTEGER, ALLOCATABLE	:: iHorsegncell(:),iVersegncell(:),iHorcellref(:,:),iVercellref(:,:)
		REAL, ALLOCATABLE	::	rA(:,:),rB(:,:),rC(:,:),rD(:,:),rE(:,:),rF(:,:),rG(:,:),rH(:,:)
		REAL,PARAMETER	:: iDeltaT=1,rDScaler=2.0 !4
		REAL, PARAMETER		:: rMort=0,rLarDiffCoef=100,rLarSurv=0.15,rLarSettle=0.15
		REAL,ALLOCATABLE	::	rDiffCoef(:)
		REAL	::	rtDiffCoef(nrow*ncol)
	
		!Taxon biological and biogeographical variables
		REAL,ALLOCATABLE	::	rtaxMinD(:),rtaxMaxD(:),rIntR(:),rMaxCatch(:),rCoral(:),rUpwelling(:)
		CHARACTER(1),ALLOCATABLE	::	cDemPel(:)
	
		!Variables for habitat suitability
		INTEGER, PARAMETER	::	iTempBinNum=48
		REAL,PARAMETER	:: rTempBinSize=1
		REAL	::	rtAbdBin(iTempBinNum),rTempBin(iTempBinNum),rtLtemp,rtUtemp
		REAL,ALLOCATABLE	:: rTempPref(:,:),rAbdBin(:,:),rUTemp(:),rLTemp(:)
		REAL,ALLOCATABLE	::	rHabSuit(:,:),rIniTotAbd(:),rK(:,:),rKMax(:),rHabSuitMax(:),rKOrig(:,:)
		REAL,ALLOCATABLE	::	rHabAssocIce(:),rHabAssocSal(:,:)
	
		!Ocean variables
		CHARACTER(40)	:: cCCScenario,cSSP,cHSmapfile
		REAL, DIMENSION(nrow*ncol)	::	pwater,lon,lat,U,V,rDX0,rDY0,rDepth,rtIce,rUpwellMap
		REAL, DIMENSION(nrow*ncol)	::	rSST,rBtT, rMaxElev,rMinElev,rAvgElev,rCoralMap,rFMortyr
		REAL	::	rSalPScale(nrow*ncol,5)	!salinity scale
		REAL,ALLOCATABLE	::	rIniSST(:,:),rIniBtT(:,:),rInshore(:),rShelf(:),rOffshore(:)
		REAL,ALLOCATABLE	::	dlResCatch(:)
		REAL	::	rFfactorHS,rFfactorEEZ
	END MODULE DBEM
	! Recursive Fortran 95 quicksort routine
	! sorts real numbers into ascending numerical order
	! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
	! Based on algorithm from Cormen et al., Introduction to Algorithms,
	! 1997 printing
	
	! Made F conformant by Walt Brainerd
	
	MODULE qsort_c_module
	
	IMPLICIT NONE
	
	public :: QsortC
	private :: Partition
	
	contains
	
	recursive subroutine QsortC(A)
	  DOUBLE PRECISION, intent(in out), dimension(:) :: A
	  integer :: iq
	
	  if(size(A) > 1) then
		 call Partition(A, iq)
		 call QsortC(A(:iq-1))
		 call QsortC(A(iq:))
	  endif
	end subroutine QsortC
	
	subroutine Partition(A, marker)
	  DOUBLE PRECISION, intent(in out), dimension(:) :: A
	  integer, intent(out) :: marker
	  integer :: i, j
	  DOUBLE PRECISION :: temp
	  DOUBLE PRECISION :: x      ! pivot point
	  x = A(1)
	  i= 0
	  j= size(A) + 1
	
	  do
		 j = j-1
		 do
			if (A(j) <= x) exit
			j = j-1
		 end do
		 i = i+1
		 do
			if (A(i) >= x) exit
			i = i+1
		 end do
		 if (i < j) then
			! exchange A(i) and A(j)
			temp = A(i)
			A(i) = A(j)
			A(j) = temp
		 elseif (i == j) then
			marker = i+1
			return
		 else
			marker = i
			return
		 endif
	  end do
	
	END subroutine Partition
	
	END MODULE qsort_c_module
	
	PROGRAM ADR
		!Main program for running DBEM
		USE qsort_c_module
		USE DBEM, ONLY :	pwater,lon,lat,U,V,ncol,nrow,rDX0,rDY0,rAbd,inTimeStep,iDeltaT,iSpp
		USE DBEM, ONLY : rTempPref,rTempBin,rAbdBin,rUTemp,rLTemp,iTempBinNum,SppNo,rAbdMaster
		USE DBEM, ONLY	:	rTempBinSize,rtAbdBin,rtUTemp,rtLTemp,cfolderpath,rHabSuit,rSST,rFMortyr
		USE DBEM, ONLY	:	iIniStartYear,iIniEndYear,iSimStartYear,iSimEndYear,rtDiffCoef,rDiffCoef
		USE DBEM, ONLY  :       iIniStartYear2,iIniEndYear2,MPAScStartYear
		USE DBEM, ONLY	:	rIniTotAbd,rK,rKOrig,rDepth,rBtT,rtaxMaxD,rtaxMinD,rMaxElev,rMinElev,rAvgElev
		USE DBEM, ONLY	: rIntR,rKMax,cDemPel,rDScaler,rtIce,rHabAssocIce,iTaxonList,iLarTimeStep,rLarv
		USE DBEM,ONLY	: rIniSST,rIniBtT,rHabSuitMax,rSalPScale,rHabAssocSal
		USE DBEM,ONLY	:	rMaxCatch,rCoral,rUpwelling,cCCScenario,cResultspath,cTaxonDatapath,rCoralMap,rUpWellMap,cSSP,cHSmapfile
		USE DBEM,ONLY	:	rInshore,rShelf,rOffshore,dlResCatch,rFfactorHS,rFfactorEEZ
	
	
		IMPLICIT NONE
		CHARACTER(40)	:: fname
		INTEGER	:: i,j,itemp,islen,inTime,isppcnt,iMonStep,iTotCellNum
		REAL	::	rTempProb,rtSST(nrow*ncol),rtCR,rDeltaAbd,rUDummy(nrow*ncol),rVDummy(nrow*ncol)
		REAL	::	rSTDummy(nrow*ncol)
		CHARACTER(4)	:: cTemp
		CHARACTER(10)	::	cTemp2
		CHARACTER(80)	:: cTemp3
		REAL,ALLOCATABLE	::	rtLarvMax(:),rKMin(:),riniAbd(:,:)
		INTEGER,PARAMETER	::	iSizeBinNum=40
		REAL	::	rWinfDummy,rVBonKDummy,rWmatDummy,rLinfDummy
		INTEGER	::	iTmaxDummy
		REAL,ALLOCATABLE	::	rWinf(:),rVBonK(:),rLwA(:),rLwB(:),rLinf(:),rNMortDummy(:),rNMort(:)
		DOUBLE PRECISION,ALLOCATABLE	::	rOrigO2Btm(:),rNewO2Btm(:),rOrigHBtm(:),rOrigHSurf(:)
		DOUBLE PRECISION,ALLOCATABLE	::	rNewO2Surf(:),rOrigO2Surf(:),rrNewO2Btm(:),rNewHBtm(:),rNewHSurf(:), &
		rOrigMeanW(:,:),rMeanW(:),rprop(:),rnewprop(:),rproptemp(:)
		REAL,ALLOCATABLE	::	 rArrhCoef(:,:),rAlloScal(:,:),rHQ10Met(:),rKVar(:),rMortTemp(:)
	
		!!--------------- TRAVIS1 -------------!!
		!!------------- Set up REAL objects for pH adult mortality and larval mortality ----------!!
		REAL,ALLOCATABLE	::	 rHQ10Mort(:),rHQ10LarvMort(:)
		!!--------------- TRAVIS1 -------------!!
	
		REAL,ALLOCATABLE	::	rOrigSST(:),rOrigBtT(:),rOrigH(:)
		LOGICAL ::	lLinfExceed,lFixedBin
		REAL	::	rlencl(iSizeBinNum),rWgtcl(iSizeBinNum)
		REAL,ALLOCATABLE	::	rFMort(:,:),rFMortScen(:,:),rLarvSortTemp(:),rMCorrect(:) !new variable Travis0,
		REAL,ALLOCATABLE	::	rHSmap(:),rHSMPAmap(:)  !rHSMPAmap for high seas MPA analysis
		DOUBLE PRECISION,ALLOCATABLE	::	dTCatch(:),dTCatchTemp(:)
		DOUBLE PRECISION,ALLOCATABLE ::	nppthreshold(:)
		INTEGER	::	iqcount
		REAL,ALLOCATABLE	::	rGCoef(:)
		DOUBLE PRECISION	::	dlAnnCatch
		LOGICAL	::	lfileexist1,lfileexist2
		!OpenMP variables
		INTEGER :: NTHREADS, TID,  OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS, CHUNK
		INTEGER, PARAMETER	:: CHUNKSIZE=1
	
		REAL, PARAMETER	::	senfactor=1.0, senfactorDD=1.00
	
		lFixedBin=.FALSE.
	
		iTotCellNum=nrow*ncol
	
		!For cedar version
		cfolderpath='~/projects/rrg-wailung/Data/'
		
		!Read basic settings
		CALL ReadSettings('../scripts/run_dbem/settings.txt')
	
		!************************************************************************
		!Allocate variables dimension
		ALLOCATE(rTempPref(SppNo,iTempBinNum))
		ALLOCATE(rAbdBin(SppNo,iTempBinNum))
	
		ALLOCATE(rAbdMaster((nrow*ncol),SppNo))
		ALLOCATE(rHabSuit((nrow*ncol),SppNo))
		ALLOCATE(rK((nrow*ncol),SppNo))
		ALLOCATE(rKOrig(nrow*ncol,SppNo))
		ALLOCATE(rKMin(SppNo))
		ALLOCATE(riniAbd(nrow*ncol,SppNo))
	
		ALLOCATE(rUTemp(SppNo))
		ALLOCATE(rLTemp(SppNo))
		ALLOCATE(rIniTotAbd(SppNo))
		ALLOCATE(rtaxMinD(SppNo))
		ALLOCATE(rtaxMaxD(SppNo))
		ALLOCATE(rDiffCoef(SppNo))
		ALLOCATE(rIntR(SppNo))
		ALLOCATE(rKMax(SppNo))
		ALLOCATE(cDemPel(SppNo))
		ALLOCATE(rHabAssocIce(SppNo))
		ALLOCATE(rHabAssocSal(SppNo,5))
		ALLOCATE(iTaxonList(SppNo))
		ALLOCATE(rtLarvMax(SppNo))
		ALLOCATE(rHabSuitMax(SppNo))
		ALLOCATE(rIniSST(nrow*ncol,(iIniEndYear-iIniStartYear+1)))
		ALLOCATE(rIniBtT(nrow*ncol,(iIniEndYear-iIniStartYear+1)))
		ALLOCATE(rOrigSST(nrow*ncol))
		ALLOCATE(rOrigBtT(nrow*ncol))
		ALLOCATE(rprop(nrow*ncol))
		ALLOCATE(rproptemp(nrow*ncol))
		ALLOCATE(rnewprop(nrow*ncol))
		ALLOCATE(rLarvSortTemp(nrow*ncol))
		ALLOCATE(rGCoef(nrow*ncol))
		ALLOCATE(dTCatch(nrow*ncol))
		ALLOCATE(dTCatchTemp(nrow*ncol))
	
	
		ALLOCATE(rLinf(SppNo))
		ALLOCATE(rWinf(SppNo))
		ALLOCATE(rVBonK(SppNo))
		ALLOCATE(rLwA(SppNo))
		ALLOCATE(rLwB(SppNo))
		ALLOCATE(rKVar(SppNo))
		ALLOCATE(rNMort(iTotCellNum))
		ALLOCATE(rOrigO2Surf(iTotCellNum))
		ALLOCATE(rOrigO2Btm(iTotCellNum))
		ALLOCATE(rNewO2Surf(iTotCellNum))
		ALLOCATE(rNewO2Btm(iTotCellNum))
		ALLOCATE(rOrigHSurf(iTotCellNum))
		ALLOCATE(rNewHSurf(iTotCellNum))
		ALLOCATE(rOrigHBtm(iTotCellNum))
		ALLOCATE(rNewHBtm(iTotCellNum))
		ALLOCATE(rNMortDummy(iSizeBinNum))
		ALLOCATE(rArrhCoef(2,SppNo))
		ALLOCATE(rAlloScal(2,SppNo))
		ALLOCATE(rHQ10Met(SppNo))
		ALLOCATE(nppthreshold(SppNo))
	
		!!--------------- TRAVIS2 -------------!!
		!!--------------- Allocate dimensions to pH adult mortality and larval mortality -------------!!
		ALLOCATE(rHQ10Mort(SppNo))
		ALLOCATE(rHQ10LarvMort(SppNo))
		!!--------------- TRAVIS2 -------------!!
	
		ALLOCATE(rOrigMeanW(nrow*ncol,SppNo))
		ALLOCATE(rMeanW(nrow*ncol))
		ALLOCATE(rMaxCatch(SppNo))
		ALLOCATE(rCoral(SppNo))
		ALLOCATE(rUpwelling(SppNo))
		ALLOCATE(rInshore(SppNo))
		ALLOCATE(rOffshore(SppNo))
		ALLOCATE(rShelf(SppNo))
		ALLOCATE(dlResCatch(nrow*ncol))
		ALLOCATE(rFMort(nrow*ncol,SppNo))
		ALLOCATE(rMortTemp(nrow*ncol))
		ALLOCATE(rHSmap(nrow*ncol))
		ALLOCATE(rHSMPAmap(nrow*ncol))
		  ALLOCATE(rMCorrect(SppNo)) !Travis
		ALLOCATE(rFMortScen(nrow*ncol,SppNo)) ! Different MPA Scenario
		!************************************************************************
	
		WRITE (*,*)	'Initializing, loading files...'
		!Read species list
		CALL	ReadTaxonList(cfolderpath)
	
		WRITE(*,*)	SppNo
	
		fname='Environment/WArea.txt'
		CALL ReadCSV(fname,pwater,nrow*ncol,-9999)
	
		CALL MapSeg
	
		fname='Environment/Lat_Lon.txt'
		CALL ReadLatLon(fname)
	
		!Read grid cell distance
		fname='Environment/HorDX.txt'
		CALL ReadCSV(fname,rDX0,nrow*ncol,0)
	
		fname='Environment/VerDY.txt'
		CALL ReadCSV(fname,rDY0,nrow*ncol,0)
	
		fname='Environment/EleMin.txt'
		CALL ReadCSV(fname,rMinElev,nrow*ncol,-9999)
	
		fname='Environment/EleMax.txt'
		CALL ReadCSV(fname,rMaxElev,nrow*ncol,-9999)
	
		fname='Environment/EleAvg.txt'
		CALL ReadCSV(fname,rAvgElev,nrow*ncol,-9999)
	
		fname='Environment/Coral.txt'
		CALL ReadCSV(fname,rCoralMap,nrow*ncol,-9999)
	
		fname='Environment/Upwelling.txt'
		CALL ReadCSV(fname,rUpWellMap,nrow*ncol,-9999)
	
		fname='Environment/HSMap.txt'
		CALL ReadCSV(fname,rHSMap,nrow*ncol,-9999)
	
		!JEPA: Incorporate MPA layer for fishing mortality
		!(0 = mpa; >0 to 1< mpa proportion, 1 != mpa)
		fname='MPA_Scenarios/'//TRIM(cHSmapfile)//'.txt'
		CALL ReadCSV(fname,rHSMPAmap,nrow*ncol,-9999)
	
		!Initialization to store climate data
		rtSST=0
		rSST=0
		rBtT=0
		rtIce=0
		rKVar=30.0
		rprop=0
		rproptemp=0
		rnewprop=0
		rNMort(:)=0
		rGCoef(:)=0
	
		DO inTime=iIniStartYear,iIniEndYear
			Write(cTemp, '(i4)') inTime
	
			fname='Climate/'//TRIM(cCCScenario)//'/bot_temp_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,rtSST,nrow*ncol,-9999)
			rIniBtT(:,inTime-iIniStartYear+1)=rtSST
			rbtT=rBtT+rtSST/(iIniEndYear-iIniStartYear+1)
	
			fname='Climate/'//TRIM(cCCScenario)//'/SST_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,rtSST,nrow*ncol,-9999)
			rIniSST(:,inTime-iIniStartYear+1)=rtSST
			rSST=rSST+rtSST/(iIniEndYear-iIniStartYear+1)
	
			fname='Climate/'//TRIM(cCCScenario)//'/IceExt_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,rtSST,nrow*ncol,0)
			IF(MAXVAL(rtSST)>1) THEN
					rtSST=rtSST/100
			END IF
			rtIce=rtIce+rtSST/(iIniEndYear-iIniStartYear+1)
	
			!Read advection
			fname='Climate/'//TRIM(cCCScenario)//'/AdvectionU_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,U,nrow*ncol,0)
			U=U*60*60*24/1000*iLarTimeStep
			rUDummy=rUDummy+U/(iIniEndYear-iIniStartYear+1)
	
			fname='Climate/'//TRIM(cCCScenario)//'/AdvectionV_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,V,nrow*ncol,0)
			V=(-V*60*60*24/1000*iLarTimeStep)
			rVDummy=rVDummy+V/(iIniEndYear-iIniStartYear+1)
	
			fname='Climate/'//TRIM(cCCScenario)//'/O2_surf_'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rNewO2Surf,nrow*ncol,0)
			rOrigO2Surf=rOrigO2Surf+rNewO2Surf/(iIniEndYear-iIniStartYear+1)
	
			fname='Climate/'//TRIM(cCCScenario)//'/htotal_surf_'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rNewHSurf,nrow*ncol,0)
			rOrigHSurf=rOrigHSurf+rNewHSurf/(iIniEndYear-iIniStartYear+1)
	
			fname='Climate/'//TRIM(cCCScenario)//'/O2_btm_'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rNewO2btm,nrow*ncol,0)
			rOrigO2btm=rOrigO2btm+rNewO2btm/(iIniEndYear-iIniStartYear+1)
	
			fname='Climate/'//TRIM(cCCScenario)//'/htotal_btm_'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rNewHbtm,nrow*ncol,0)
			rOrigHbtm=rOrigHbtm+rNewHbtm/(iIniEndYear-iIniStartYear+1)
	
			fname='Climate/'//TRIM(cCCScenario)//'/totalphy2'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rnewprop,nrow*ncol,0)
			rprop=rprop+rnewprop/(iIniEndYear-iIniStartYear+1)
	
			!IF(inTime==iIniStartYear) THEN
			!	fname='Effort/'//TRIM(cSSP)//'/FMort'//TRIM(ctemp)//'.txt'
			!	CALL READFMort(fname,rFMortyr,nrow*ncol,0)
			!END IF
	
			!WRITE(*,*) SUM(rFMortyr)
		END DO
	
		!Call Salinity map based on demersal or pelagic
		IF(cDemPel(iSpp)=='D') THEN
				CALL SalinityMap(iIniStartYear,iIniEndYear,'Salinity_btm_')
		ELSE
				CALL SalinityMap(iIniStartYear,iIniEndYear,'Salinity_surf_')
		END IF
	
		rNewHbtm=0
		rNewO2btm=0
		rNewHSurf=0
		rNewO2Surf=0
	
		U=rUDummy
		V=rVDummy
	
		DO iSpp=1,SppNo
		  Write(cTemp2, '(i6)') iTaxonList(iSpp)
		  fname='/Species/Distributions/S'//TRIM(cTemp2)//'.csv' 
			INQUIRE(FILE=TRIM(cfolderpath)//TRIM(fname), EXIST=lfileexist1)
			INQUIRE(FILE=TRIM(cfolderpath)//'/Species/'//TRIM(cTaxonDatapath)//"/"//TRIM(ctemp2)//".txt", EXIST=lfileexist2)
			IF (lfileexist1 .EQV. .TRUE. .AND. lfileexist2 .EQV. .TRUE.) THEN
				Write(cTemp2, '(i6)') iTaxonList(iSpp)
				WRITE(*,*)	cTemp2
	
				ctemp3='Results/'//TRIM(cResultspath)//'/'//TRIM(cTemp2)
				
				WRITE(*,*) ctemp3
	
				!Create folder to store results
	
				CALL CreateDataFolder(ctemp3)
	
				WRITE(*,*) 'Read Parameters...'
	
				!Read in taxon specific parameters
	
				!!--------------- TRAVIS3 -------------!!
				!!--------------- Add in rHQ10Mort and rHQ10LarvMort to CALL function -------------!!
				CALL ReadTaxonPara(TRIM(cTemp2)//'.txt',rLinf(iSpp),rVBonK(iSpp),rLwA(iSpp),rLwB(iSpp), &
				rArrhCoef(:,iSpp),rAlloScal(:,iSpp),rHQ10Met(iSpp), &
				rHQ10Mort(iSpp),rHQ10LarvMort(iSpp))
				!!--------------- TRAVIS3 -------------!!
	
				!sensitivity testing
				rDiffCoef(iSpp)=rDiffCoef(iSpp)*senfactor
	
				!WRITE(*,*)	rLinf(iSpp),rVBonK(iSpp),rtaxMaxD(iSpp),rtaxMinD(iSpp)
	
				rWinf(iSpp)=rLwA(iSpp)*rLinf(iSpp)**rLwB(iSpp)
	
				!testing
				!WRITE(*,*) rWinf(iSpp)
	
				!Debugging
				!WRITE(*,*)	rLinf(iSpp),rWinf(iSpp),rVBonK(iSpp),rLwA(iSpp),rLwB(iSpp), &
				!rArrhCoef(1,iSpp),rArrhCoef(2,iSpp),rAlloScal(1,iSpp),rAlloScal(2,iSpp),rHQ10Met(iSpp)
	
				fname='Species/Distributions/S'//TRIM(cTemp2)//'.csv'
				WRITE(*,*) fname
				CALL ReadCSV(fname,rAbd,nrow*ncol,0)
	
				rAbdMaster(:,iSpp)=0
	
				!Call Salinity map based on demersal or pelagic
		  !IF(cDemPel(iSpp)=='D') THEN
				!		CALL SalinityMap(iIniStartYear,iIniEndYear,'Salinity_btm_')
				!ELSE
				!		CALL SalinityMap(iIniStartYear,iIniEndYear,'Salinity_surf_')
				!END IF
	
				!************************************************
				!Convert relative abundance to absolute abundance
				!If rMaxCatch too small (<1000), increase it by 1e4 else  10 - Travis
		  IF(rMaxCatch(iSpp)<1000) THEN
				  rMCorrect(iSpp)=10000
		  ELSE
				  rMCorrect(iSpp)=10
		  END IF
	
				rAbd=(rMaxCatch(iSpp)*rMCorrect(iSpp)*4/rIntR(iSpp))*rAbd/SUM(rAbd) !Travis
	
	!***************
				!testing
				!WRITE(*,*) 'Abd=',SUM(rAbd)
				!Fishing mortality relative to natural mortality
		  !WRITE(*,*)	rFfactorHS,rFfactorEEZ
	!***************
	!Initialization fishing mortality rates for highseas (HS) and EEZs (read in from Settings file)
		  DO i=1,iTotCellNum
					IF(rHSMap(i)==1) THEN
						rFMort(i,iSpp)=rFfactorHS*rIntR(iSpp)/2
						rFMortScen(i,iSpp)=rFfactorHS*rIntR(iSpp)/2*(rHSMPAMap(i)) !JEPA: Incorportates MPA layer 
					ELSE
						rFMort(i,iSpp)=rFfactorEEZ*rIntR(iSpp)/2
						rFMortScen(i,iSpp)=rFfactorEEZ*rIntR(iSpp)/2*(rHSMPAMap(i)) !JEPA: Incorportates MPA layer
					END IF
				END DO
	
				!calculate npp factor for carrying capacity calculation
				rproptemp(:)=rprop(:)
				DO i=1,iTotCellNum
					IF(rAbd(i)==0) THEN
						rproptemp(i)=0
					END IF
				END DO
	
				CALL QsortC(rproptemp)
	
				rproptemp(:)=rproptemp((nrow*ncol):1:-1)
				iqcount=1
				Do While (rproptemp(iqcount)>0)
					iqcount=iqcount+1
				END DO
	
				!calculate 95 percentile of NPP in distribution (will only affect
				!prediction of carrying capacity at novel range)
				iqcount=CEILING(iqcount*0.05)
				nppthreshold(ispp)=rproptemp(iqcount)
				WRITE(*,*) nppthreshold(ispp)
	
			!	WRITE(*,*) SUM(rHSMPAMap)
				!rFMort(:,iSpp)=0.0
				!WRITE(*,*) 'IntR', rIntR(iSpp)/2
				!rFMort(:,iSpp)=0
				!************************************************
				riniAbd(:,iSpp)=rAbd(:)
				rAbdMaster(:,iSpp)=rAbd
	
				!Set temperature preference
				CALL SetTempPreference
				rAbdBin(iSpp,:)=rtAbdBin
				rLTemp(iSpp)=rtLTemp
				rUTemp(iSpp)=rtUTemp
				rIniTotAbd=SUM(rAbd,MASK=rAbd>0)
				IF(cDemPel(iSpp)=='D') THEN
						CALL CalK(.TRUE.,rBtT,rLinf(iSpp))
						rSTDummy=rBtT
				ELSE
	
						CALL CalK(.TRUE.,rSST,rLinf(iSpp))
						rSTDummy=rSST
				END IF
	
				!Adjusting for NPP for theoretical carrying capacity
				DO i=1,iTotCellNum
					!IF(rnewprop(i)/rprop(i)<=2) THEN
					IF(rprop(i)>0) THEN
						!rK(i,iSpp)=rK(i,iSpp)*rnewprop(i)/rprop(i)
						IF(rprop(i)<nppthreshold(iSpp)) THEN
							rK(i,iSpp)=rK(i,iSpp)*rprop(i)/nppthreshold(iSpp)
						END IF
	
					END IF
				END DO
				!store baseline reference carrying capacity
				rKOrig(:,iSpp)=rK(:,iSpp)
	
				!WRITE(*,*) 'K=', SUM(rK(:,iSpp))
				!***********************************************************************
				!loop to create size bins
	
				IF(lFixedBin .EQV. .FALSE.) THEN
					!Create variable size bin depending on species
					rlencl(1)=(rLinf(iSpp) / (iSizeBinNum-10))/2
	
					DO itemp=2,iSizeBinNum
						rlencl(itemp) = rlencl(itemp - 1) + (rLinf(iSpp) /(iSizeBinNum-10))
					END DO
					rWgtcl=rLwA(iSpp)*rlencl**rLwB(iSpp)
				ELSE
					lLinfExceed=.FALSE.
					itemp=1
					rWgtcl(1)=1.0
					Do While ((itemp<=iSizeBinNum) .AND. (lLinfExceed .EQV. .FALSE.))
						itemp=itemp+1
						!Create a log2 series of weight bin, then convert to length bin
						rWgtcl(itemp)=2**(itemp-1)+(((2**(itemp))-(2**(itemp-1)))/2)
						rlencl(itemp)=(rWgtcl(itemp)/rLwA(iSpp))**(1/rLwB(iSpp))
						IF(rWgtcl(itemp)>=rWinf(iSpp)) THEN
							lLinfExceed=.FALSE.
						END IF
					END DO
				END IF
	
				!*********************************************************************************
	
				rOrigMeanW(:,iSpp)=0
				rMeanW(:)=0
				!********************************************************************************
				!Calculate biomass per recruit(rMeanW) per cell (and size-frequency distribution)
				!Use OpenMP for parallel computing
				CHUNK=CHUNKSize
				!$OMP	PARALLEL SHARED(rOrigMeanW) PRIVATE(i,TID)
	
				TID = OMP_GET_THREAD_NUM()
				IF (TID == 0) THEN
					NTHREADS = OMP_GET_NUM_THREADS()
				!	PRINT *, 'Number of threads =', NTHREADS
				END IF
				PRINT *, 'Thread',TID,' starting...'
	
				!$OMP DO SCHEDULE(GUIDED)
				DO i=1,iTotCellNum
	
					IF(pwater(i)>0) THEN
					IF(cDemPel(iSpp)=='D') THEN
						CALL Ecophysiology(iSizeBinNum,rWinf(iSpp),rVBonK(iSpp),rLwA(iSpp),rLwB(iSpp), &
						rArrhCoef(:,iSpp),rAlloScal(:,iSpp),rHQ10Met(iSpp),rSTDummy(i),rOrigO2Btm(i), &
						rOrigHBtm(i),rSTDummy(i),rOrigO2btm(i),rOrigHbtm(i),rKVar(iSpp),rlencl,rWgtcl, &
						rWinfDummy,rVBonKDummy,rWmatDummy,rNMortDummy, rLinfDummy,iTmaxDummy,rMeanW(i),rNMort(i), &
						!!------------ TRAVIS4 --------------!!
						!!------------ Include new variable rHQ10Mort(iSpp) ------------!!
						rHQ10Mort(iSpp))
						!!------------ TRAVIS4 --------------!!
					ELSE
						CALL Ecophysiology(iSizeBinNum,rWinf(iSpp),rVBonK(iSpp),rLwA(iSpp),rLwB(iSpp), &
						rArrhCoef(:,iSpp),rAlloScal(:,iSpp),rHQ10Met(iSpp),rSTDummy(i),rOrigO2Surf(i), &
						rOrigHSurf(i),rSTDummy(i),rOrigO2Surf(i),rOrigHSurf(i),rKVar(iSpp),rlencl,rWgtcl, &
						rWinfDummy,rVBonKDummy,rWmatDummy,rNMortDummy, rLinfDummy,iTmaxDummy,rMeanW(i),rNMort(i), &
						!!------------ TRAVIS5 --------------!!
						!!------------ Include new variable rHQ10Mort(iSpp) ------------!!
						rHQ10Mort(iSpp))
						!!------------ TRAVIS5 --------------!!
					END IF
						!Debugging
						!IF (ISNAN(rMeanW(i))) THEN
						!	WRITE(*,*)	'Winf',rWinfDummy,rLinfDummy,rVBonKDummy,rWmatDummy,rNMortDummy
						!END IF
				!		WRITE(*,100) TID,i,rOrigMeanW(i,iSpp)
				 !100    FORMAT(' Thread',I2,': C(',I5,')=',F8.2)
						!Abundance in weight into number
					END IF
	
				END DO
	
				!$OMP END DO
				!!NOWAIT
				!$OMP BARRIER
				!$OMP END PARALLEL
				!*********************************************************************************
				rOrigMeanW(:,iSpp)=rMeanW(:)
				!WRITE(*,*) SUM(rOrigMeanw(:,iSpp))
			  fname='testMeanW.csv'
				CALL WriteCSVDouble(fname,rOrigMeanW(:,iSpp),nrow*ncol)
	
				!WRITE(*,*) rOrigMeanW(101795,iSpp)
				!Calculate larval dispersal
	
				!!--------------- TRAVIS6 -------------!!
				!!--------------- Insert variables for pH effects on larval mortality ---------------!!
				CALL CalLarDisp(rHQ10LarvMort(iSpp),rOrigHSurf,rOrigHSurf)
				!!--------------- TRAVIS6 -------------!!
	
				rLarvSortTemp(:)=0
				!Store maximum larval dispersal
				!rtLarvMax(iSpp)=MAXVAL(rLarv)
	
			  !change to calculate 1% quantile
			  rLarvSortTemp(:)=rLarv(:)
			  !CALL QsortC(rLarvSortTemp)
			  !rLarvSortTemp(:)=rLarvSortTemp((nrow*ncol):1:-1)
			  !iqcount=1
			  !Do While (rLarvSortTemp(iqcount)>0)
				!	iqcount=iqcount+1
			  !END DO
	
				!rtLarvMax(iSpp)=rLarvSortTemp(floor(iqcount*0.99))
				!WRITE(*,*) rtLarvMax(iSpp)
	
				!Use a fixed threshold
				rtLarvMax(iSpp)=0.0000001
	
				!testing
				fname='testrk1.csv'
				CALL WriteCSV(fname,rK(:,iSpp),nrow*ncol)
	
				DO i=1,nrow*ncol
					IF(rLarv(i)<rtLarvMax(iSpp)) THEN
						rLarv(i)=0.0
					END IF
				END DO
	
				!WRITE(*,*)	rtLarvMax(iSpp)
			  rKMin(iSpp)=MAXVAL(rK(:,iSpp))
			  DO i=1,iTotCellNum
					IF(rAbd(i)>0 .AND. rK(i,iSpp)>0) THEN
						IF(rK(i,iSpp)<rKMin(iSpp)) THEN
							rKMin(iSpp)=rK(i,iSpp)
						END IF
					END IF
			 END DO
	
				!testing
				fname='testLarv.csv'
				CALL WriteCSV(fname,rLarv,nrow*ncol)
	
				fname='testrk2.csv'
				CALL WriteCSV(fname,rK(:,iSpp),nrow*ncol)
				!WRITE(*,*)	SUM(rK(:,iSpp))
				fname='testini.csv'
				CALL WriteCSV(fname,rAbd,nrow*ncol)
	
			END IF
		END DO
	
		DEALLOCATE(rIniSST)
		DEALLOCATE(rInIBtT)
	
		rOrigSST=rSST
		rOrigBtT=rBtT
		rnewprop=rprop
		!Spin-up
		DO j=1,4	!number of cycle to spin-up the distribution
			DO inTime=iIniStartYear2,iIniEndYear2
				WRITE(*,*) 'Time step = ', inTime
	
				Write(cTemp, '(i4)') inTime
			!	fname='Climate/'//TRIM(cCCScenario)//'/SST_'//TRIM(ctemp)//'.txt'
			!	CALL ReadCSV(fname,rSST,nrow*ncol,-9999)
			rSST=rOrigSST
			!	fname='Climate/'//TRIM(cCCScenario)//'/botTemp_'//TRIM(ctemp)//'.txt'
			!	CALL ReadCSV(fname,rBtT,nrow*ncol,-9999)
			rBtT=rOrigBtT
			!Read advection
			!fname='Climate/'//TRIM(cCCScenario)//'/AdvectionU_'//TRIM(ctemp)//'.txt'
			!CALL ReadCSV(fname,U,nrow*ncol,0)
			!U=U*60*60*24/1000*iLarTimeStep
	
			!fname='Climate/'//TRIM(cCCScenario)//'/AdvectionV_'//TRIM(ctemp)//'.txt'
			!CALL ReadCSV(fname,V,nrow*ncol,0)
			!V=(-V*60*60*24/1000*iLarTimeStep)
	
			!Ice extent
			!fname='Climate/'//TRIM(cCCScenario)//'/IceExt_'//TRIM(ctemp)//'.txt'
			!CALL ReadCSV(fname,rtIce,nrow*ncol,0)
			!IF(MAXVAL(rtIce)>1) THEN
				!	rtIce(:)=rtIce(:)/100
			!END IF
			!	fname='Climate/'//TRIM(cCCScenario)//'/O2_surf_'//TRIM(ctemp)//'.txt'
			!	CALL ReadCSVDouble(fname,rNewO2Surf,nrow*ncol,0)
			rNewO2Surf=rOrigO2Surf
			!	fname='Climate/'//TRIM(cCCScenario)//'/htotal_surf_'//TRIM(ctemp)//'.txt'
			!	CALL ReadCSVDouble(fname,rNewHSurf,nrow*ncol,0)
			rNewHSurf=rOrigHSurf
			!	fname='Climate/'//TRIM(cCCScenario)//'/O2_btm_'//TRIM(ctemp)//'.txt'
			!	CALL ReadCSVDouble(fname,rNewO2btm,nrow*ncol,0)
			rNewO2btm=rOrigO2btm
			!	fname='Climate/'//TRIM(cCCScenario)//'/htotal_btm_'//TRIM(ctemp)//'.txt'
			!	CALL ReadCSVDouble(fname,rNewHbtm,nrow*ncol,0)
			rNewHbtm=rOrigHbtm
			!	fname='Climate/'//TRIM(cCCScenario)//'/totalphy2'//TRIM(ctemp)//'.txt'
			!	CALL ReadCSVDouble(fname,rnewprop,nrow*ncol,0)
			rnewprop=rprop
	
				DO iSpp=1,SppNo
					dTCatch(:)=0.0
					!Check to see if base distribution file exists and taxon information exists
					Write(cTemp2, '(i6)') iTaxonList(iSpp)
					fname='Species/Distributions/S'//TRIM(cTemp2)//'.csv' ! changed to Mel species distribution folder
					INQUIRE(FILE=TRIM(cfolderpath)//TRIM(fname), EXIST=lfileexist1)
					INQUIRE(FILE=TRIM(cfolderpath)//'/Species/'//TRIM(cTaxonDatapath)//"/"//TRIM(ctemp2)//".txt", EXIST=lfileexist2)
					IF (lfileexist1 .EQV. .TRUE. .AND. lfileexist2 .EQV. .TRUE.) THEN
						rAbd=rAbdMaster(:,iSpp)
	
						!WRITE(*,*) SUM(rAbd)
	
						!testing
						!WRITE(ctemp2,'(i6)') iTaxonlist(iSpp)
	
						!WRITE(*,*) ctemp2
	
						!CALL ReadCatch(inTime,ctemp2,dlAnnCatch)
	
						!WRITE(*,*) dlAnnCatch
	
						!CALL AllocCatch(rAbd,dlAnnCatch)
	
						IF(cDemPel(iSpp)=='D') THEN
							CALL CalK(.FALSE.,rBtT,rLinf(iSpp))
						ELSE
							CALL CalK(.FALSE.,rSST,rLinf(iSpp))
						END IF
	
						!Adjust for NPP in calculating theoretical carrying capacity
						DO i=1,iTotCellNum
								!IF(rnewprop(i)/rprop(i)<=2) THEN
								IF(rnewprop(i)>0) THEN
										!rK(i,iSpp)=rK(i,iSpp)*rnewprop(i)/rprop(i)
										IF(rnewprop(i)<nppthreshold(iSpp)) THEN
														rK(i,iSpp)=rK(i,iSpp)*rnewprop(i)/nppthreshold(iSpp)
										END IF
	
									END IF
						END DO
	
						DO i=1,iTotCellNum
							IF(riniAbd(i,iSpp)==0) THEN
								IF((rK(i,iSpp)-rKOrig(i,iSpp))>(rKMin(iSpp)*0.5)) THEN
									rK(i,iSpp)=rK(i,iSpp)-rKOrig(i,iSpp)
								ELSE
									rK(i,iSpp)=0
								END IF
							END IF
						END DO
	
						!WRITE(*,*)	'K1',SUM(rK(:,iSpp))
						!***********************************************************************
						!loop to create size bins
	
						IF(lFixedBin .EQV. .FALSE.) THEN
							!Create variable size bin depending on species
							rlencl(1)=(rLinf(iSpp) / (iSizeBinNum-10))/2
	
							DO itemp=2,iSizeBinNum
								rlencl(itemp) = rlencl(itemp - 1) + (rLinf(iSpp) /(iSizeBinNum-10))
							END DO
							rWgtcl=rLwA(iSpp)*rlencl**rLwB(iSpp)
						ELSE
							lLinfExceed=.FALSE.
							itemp=1
							rWgtcl(1)=1.0
							Do While ((itemp<=iSizeBinNum) .AND. (lLinfExceed .EQV. .FALSE.))
								itemp=itemp+1
								!Create a log2 series of weight bin, then convert to length bin
								rWgtcl(itemp)=2**(itemp-1)+(((2**(itemp))-(2**(itemp-1)))/2)
								rlencl(itemp)=(rWgtcl(itemp)/rLwA(iSpp))**(1/rLwB(iSpp))
								IF(rWgtcl(itemp)>=rWinf(iSpp)) THEN
									lLinfExceed=.FALSE.
								END IF
							END DO
						END IF
	
						!*********************************************************************************
	
						rMeanW(:)=0
						iqcount=0
						!********************************************************************************
						!Calculate biomass per recruit(rMeanW) per cell (and size-frequency distribution)
						!Use OpenMP for parallel computing
						!WRITE(*,*)	'k2',SUM(rK(:,iSpp))
						CHUNK=CHUNKSize
						!$OMP	PARALLEL SHARED(rOrigMeanW) PRIVATE(i,TID)
	!					!$OMP	PARALLEL SHARED(rMeanW,rWinf,rVBonK,rLwA,rLwB, &
	!					!$OMP&		rArrhCoef,rAlloScal,rHQ10Met, &
	!					!$OMP&		rOrigHBtm,rSTDummy,rOrigO2btm,rKVar,rlencl,rWgtcl, &
	!					!$OMP&		pwater,cDemPel,NTHREADS,CHUNK,rNMort)  PRIVATE(i,rWinfDummy,rVBonKDummy, &
	!					!$OMP& 		rWmatDummy,rNMortDummy, rLinfDummy,iTmaxDummy,TID)
	
						TID = OMP_GET_THREAD_NUM()
						IF (TID == 0) THEN
							NTHREADS = OMP_GET_NUM_THREADS()
						!	PRINT *, 'Number of threads =', NTHREADS
						END IF
						!PRINT *, 'Thread',TID,' starting...'
	
						!$OMP DO SCHEDULE(GUIDED)
						DO i=1,iTotCellNum
	
							IF(rK(i,iSpp)>0) THEN
								IF(cDemPel(iSpp)=='D') THEN
	
									CALL Ecophysiology(iSizeBinNum,rWinf(iSpp),rVBonK(iSpp),rLwA(iSpp),rLwB(iSpp), &
									rArrhCoef(:,iSpp),rAlloScal(:,iSpp),rHQ10Met(iSpp),rOrigBtT(i),rOrigO2Btm(i),&
									rOrigHBtm(i),rBtT(i), rNewO2Btm(i),rNewHBtm(i),rKVar(iSpp),rlencl,rWgtcl, &
									rWinfDummy,rVBonKDummy,rWmatDummy,rNMortDummy,rLinfDummy,iTmaxDummy,rMeanW(i),rNMort(i), &
									!!------------ TRAVIS7 --------------!!
									!!------------ Include new variable rHQ10Mort(iSpp) ------------!!
									rHQ10Mort(iSpp))
									   !!------------ TRAVIS7 --------------!!
								ELSE
	
									CALL Ecophysiology(iSizeBinNum,rWinf(iSpp),rVBonK(iSpp),rLwA(iSpp),rLwB(iSpp), &
									rArrhCoef(:,iSpp),rAlloScal(:,iSpp),rHQ10Met(iSpp),rOrigSST(i),rOrigO2Surf(i), &
									rOrigHSurf(i),rSST(i),rNewO2Surf(i),rNewHSurf(i),rKVar(iSpp),rlencl,rWgtcl, &
									rWinfDummy,rVBonKDummy,rWmatDummy,rNMortDummy,rLinfDummy,iTmaxDummy,rMeanW(i),rNMort(i), &
									!!------------ TRAVIS8 --------------!!
									!!------------ Include new variable rHQ10Mort(iSpp) ------------!!
									rHQ10Mort(iSpp))
									!!------------ TRAVIS8 --------------!!
								END IF
	
								!Debugging
								!WRITE(*,*)	'Winf',rWinfDummy,rLinfDummy,rVBonKDummy,rWmatDummy,rNMortDummy
						!		WRITE(*,100) TID,i,rMeanW(i,iSpp)
					 !100    	FORMAT(' Thread',I2,': C(',I5,')=',F8.2)
								!Abundance in weight into number
								!rNMort(i)=rNMortDummy(1)/inTimeStep
								!IF(ISNAN(rNMort(i))) THEN
								!	WRITE(*,*) i, rNMortDummy,inTimeStep,rK(i,iSpp),rBtT(i),rMeanW(i)
								!	STOP
								!END IF
							END IF
	
						END DO
	
						!$OMP END DO
						!!NOWAIT
						!$OMP BARRIER
						!$OMP END PARALLEL
						rNMort(:)=rNMort(:)/inTimeStep
						!Delete and put in BPR factor in habitat suitability
						!WRITE(*,*)	'K3',SUM(rK(:,iSpp))
						!WRITE(*,*) SUM(rOrigMeanW(:,iSpp)),SUM(rMeanW(:))
	
	
						DO i=1,iTotCellNum
							IF (rOrigMeanW(i,iSpp)>0) THEN
								rK(i,iSpp)=rK(i,iSpp)*rMeanW(i)/rOrigMeanW(i,iSpp)
							ELSE
								rOrigMeanW(i,iSpp)=rMeanW(i)
							END IF
						END DO
						!WRITE(*,*)	'K4',SUM(rK(:,iSpp))
	
	
						!*********************************************************************************
	
	
						!debugging
						!fname='testrk.csv'
						!CALL WriteCSV(fname,rK(:,iSpp),nrow*ncol)
	
						!Calculate larval dispersal
						!!--------------- TRAVIS9 -------------!!
						!!--------------- Insert variables for pH effects on larval mortality ---------------!!
						CALL CalLarDisp(rHQ10LarvMort(iSpp),rOrigHSurf,rNewHSurf)
						!!--------------- TRAVIS9 -------------!!
						rNMort(:)=0.0
						rGCoef(:)=1.0
	
						!DO i=1,nrow*ncol
						!	IF(rLarv(i)>rtLarvMax(iSpp)) THEN
								!rK(i,iSpp)=0.0
								!rLarv(i)=0.0
						!		rNMort(i)=0.0
						!		rGCoef(i)=1.0
						!	ELSE
						!		rNMort(i)=1.0
						!		rGCoef(i)=1.0
						!	END IF
						!END DO
	
						!recreating map segments by accounting for rK
						CALL MapSeg2
						!WRITE(*,*)	'K5',SUM(rK(:,iSpp))
						!rtDiffCoef=rDiffCoef
						U=0
						V=0
						dlAnnCatch=0.0
						DO iMonStep=1,inTimeStep
							dTCatchTemp(:)=0.0
							DO i=1,nrow*ncol
								rtCR=0
								IF(rK(i,iSpp)>0) THEN
									rtCR=MAX(0.0,(1-rAbd(i)/rK(i,iSpp)))
									rtCR=Min(1.0,rtCR)
									!Add BPR factor to determine change in diffcoef
									!rtDiffCoef(i)=1/(1+rDScaler*rHabSuit(i,iSpp)*rtCR*rMeanW(i)/rOrigMeanW(i,iSpp))
									!rtDiffCoef(i)=1/(1+rDScaler*rHabSuit(i,iSpp)*rtCR)
									rtDiffCoef(i)=2/(1+EXP(rDScaler*rHabSuit(i,iSpp)*rtCR)) !1/(1+rDScaler*rHabSuit(i,iSpp)*rtCR)
									rtDiffCoef(i)=rDiffCoef(iSpp)*rtDiffCoef(i)
									!rtDiffCoef(i)=rDiffCoef(iSpp)*rtDiffCoef(i)+1
								ELSE
									rtDiffCoef(i)=rDiffCoef(iSpp)
								END IF
	
							END DO
							!recreating map segments by accounting for rK
	
							!testing
							!WRITE(*,*)	iMonStep,rAbd(79084),rtDiffCoef(79084),rHabSuit(79084,iSpp)
	
							CALL setabcdef
	
							CALL caltridag
	
							!WRITE(*,*)	SUM(rAbd),SUM(rk(:,iSpp))
							!WRITE(*,*)	'653...',SUM(rAbd),SUM(rk(:,iSpp))
	
							!WRITE(*,*)	'Abd=', SUM(rAbd(:))
							!CALL WriteCSV('rabdtest.csv',rAbd,259200)
							!rNMort(:)=0.0
							dTCatchTemp(:)=0.0
	
							!testing
							!WRITE(*,*)	rAbd(36999)/rK(36999,iSpp),rHabSuit(36999,iSpp),rtDiffCoef(36999)
	
							!testing
							!WRITE(*,*)	iMonStep,rAbd(79084),rFMort(79084,iSpp),rNMort(79084)
							CALL CalDeltaAbd(rFMort(:,iSpp),rNMort,rGCoef,dTCatchTemp)
							!WRITE(*,*)	'Abd2=',SUM(rAbd(:)) !,rFMort(iSpp)
	
							!test
							!WRITE(*,*)	iMonStep,rAbd(72930)
							dTCatch(:)=dTCatch(:)+dTCatchTemp(:)
						End DO
	
						WRITE(*,*)	'Catch=', SUM(dTCatch(:))
						DO i= 1,iTotCellNum
							IF(rAbd(i)==0 .AND. rLarv(i)>0 .AND. rK(i,iSpp)>0) THEN
								rAbd(i)=1
							END IF
							!testing
							!IF(i==72930) THEN
							!	WRITE(*,*) rAbd(i),rk(i,iSpp),rtDiffCoef(i),rMeanW(i)/rOrigMeanW(i,iSpp),rnewprop(i)/rprop(i)
							!END IF
						END DO
	
						!WRITE(*,*) SUM(rAbd(:))
						rAbdMaster(:,iSpp)=rAbd
						Write(cTemp2, '(i6)') iTaxonList(iSpp)
						fname=TRIM(cTemp2)//'Abd'//TRIM(ctemp)//'.txt'
						CALL WriteCSV(fname,rAbd,nrow*ncol)
						WRITE(*,*)	'Abd=',SUM(rAbd(:))
						!fname='MeanW'//TRIM(ctemp)//'.csv'
						!CALL WriteCSVDouble(fname,rMeanW,nrow*ncol)
	
						!fname='OrigMeanW'//TRIM(ctemp)//'.csv'
						!CALL WriteCSVDouble(fname,rOrigMeanW(:,iSpp),nrow*ncol)
	
						rAbd=0
						rOrigMeanW(:,iSpp)=rMeanW(:)
					END IF
				END DO
			END DO
		END DO
		!************************************************
		!Convert relative abundance to absolute abundance
	!turning this off for testing
	!	Do iSpp=1,SppNo
	!		Write(cTemp2, '(i6)') iTaxonList(iSpp)
	!	    fname='Distributions/S'//TRIM(cTemp2)//'.csv'! Changed to Mel spp distribution folder
	!		INQUIRE(FILE=TRIM(cfolderpath)//TRIM(fname), EXIST=lfileexist1)
	!		INQUIRE(FILE=TRIM(cfolderpath)//"TaxonDataNE/"//TRIM(ctemp2)//".txt", EXIST=lfileexist2)
	!		IF (lfileexist1 .EQV. .TRUE. .AND. lfileexist2 .EQV. .TRUE.) THEN
	!			rAbdMaster(:,iSpp)=(rMaxCatch(iSpp)*4/rIntR(iSpp)*1000)*rAbdMaster(:,iSpp)/SUM(rAbdMaster(:,iSpp))
	!			!testing
	!			!WRITE(*,*) SUM(rAbdMaster(:,iSpp)
	!		ELSE
	!			rAbdMaster(:,iSpp)=0
	!		END IF
	!	END DO
	!
		!************************************************
		DO inTime=iSimStartYear,iSimEndYear
			WRITE(*,*) 'Time step = ', inTime
	
			Write(cTemp, '(i4)') inTime
			fname='Climate/'//TRIM(cCCScenario)//'/SST_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,rSST,nrow*ncol,-9999)
	
			fname='Climate/'//TRIM(cCCScenario)//'/bot_temp_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,rBtT,nrow*ncol,-9999)
	!rSST=rOrigSST
	!rBtT=rOrigBtT
			!Read advection
			fname='Climate/'//TRIM(cCCScenario)//'/AdvectionU_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,U,nrow*ncol,0)
			U=U*60*60*24/1000*iLarTimeStep
	
			fname='Climate/'//TRIM(cCCScenario)//'/AdvectionV_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,V,nrow*ncol,0)
			V=(-V*60*60*24/1000*iLarTimeStep)
	
			!Ice extent
			fname='Climate/'//TRIM(cCCScenario)//'/IceExt_'//TRIM(ctemp)//'.txt'
			CALL ReadCSV(fname,rtIce,nrow*ncol,0)
	
			fname='Climate/'//TRIM(cCCScenario)//'/O2_surf_'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rNewO2Surf,nrow*ncol,0)
	
			fname='Climate/'//TRIM(cCCScenario)//'/htotal_surf_'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rNewHSurf,nrow*ncol,0)
	
			fname='Climate/'//TRIM(cCCScenario)//'/O2_btm_'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rNewO2btm,nrow*ncol,0)
	
			fname='Climate/'//TRIM(cCCScenario)//'/htotal_btm_'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rNewHbtm,nrow*ncol,0)
	
			fname='Climate/'//TRIM(cCCScenario)//'/totalphy2'//TRIM(ctemp)//'.txt'
			CALL ReadCSVDouble(fname,rnewprop,nrow*ncol,0)
	
			IF(MAXVAL(rtIce)>1) THEN
				rtIce=rtIce/100
			END IF
	
			!Call Salinity map based on demersal or pelagic
		IF(cDemPel(iSpp)=='D') THEN
					CALL SalinityMap(iIniStartYear,iIniEndYear,'Salinity_btm_')
			ELSE
			CALL SalinityMap(iIniStartYear,iIniEndYear,'Salinity_surf_')
			END IF
	
			!fname='Effort/'//TRIM(cSSP)//'/FMort'//TRIM(ctemp)//'.txt'
			!CALL READFMort(fname,rFMortyr,nrow*ncol,0)
	
	
			DO iSpp=1,SppNo
				dTCatch(:)=0.0
				Write(cTemp2, '(i6)') iTaxonList(iSpp)
				fname='Species/Distributions/S'//TRIM(cTemp2)//'.csv'
				INQUIRE(FILE=TRIM(cfolderpath)//TRIM(fname), EXIST=lfileexist1)
				INQUIRE(FILE=TRIM(cfolderpath)//'/Species/'//TRIM(cTaxonDatapath)//"/"//TRIM(ctemp2)//".txt", EXIST=lfileexist2)
				IF (lfileexist1 .EQV. .TRUE. .AND. lfileexist2 .EQV. .TRUE.) THEN
					rAbd=rAbdMaster(:,iSpp)
	
					!Read in catch
					WRITE(ctemp2,'(i6)') iTaxonlist(iSpp)
	
	
					!! ----------- JEPA_Scenarios------------ !!
					! Set conservation rules based on time frame
					
					IF(inTime < MPAScStartYear) THEN ! if time is before MPA start year
						rFMort = rFMort ! No MPAs
					!ELSE IF(inTime >= 2020 .AND. inTime <= 2025) THEN ! For time in between 2020 and 2025
					!	rFMort=rFMortRef ! Reference scenario kicks in (Current MPA protection)
					!ELSE IF(inTime > 2025 .AND. inTime <= 2030) THEN ! Between 2025 and 2030
					!		rFMort = rFMortScenHalf  ! 1/2 of the B30 scenario kicks in
					ELSE IF (inTime >= MPAScStartYear) THEN ! After MPA start year
						rFMort = rFMortScen ! Normal B30 scenario
					END IF		
				
					!! ------------------------------------ !!
	
					!rFMort(:,iSpp)=0
					!IF (inTime<=2006) THEN
						!CALL ReadCatch(inTime,ctemp2,dlAnnCatch)
	
						!WRITE(*,*) 'Annual catch',dlAnnCatch
	
						!Turn off catch
					!	dlAnnCatch=0.0
					!	CALL AllocCatch(rAbd,dlAnnCatch)
	
					!	rFMort(iSpp)=dlAnnCatch/SUM(rAbd(:))
											!test
					!Read fishing mortality
					!DO i=1,iTotCellNum
					!	rFMort(i,iSpp)=rIntR(iSpp)/2
					!END DO
	
					!END IF
	
					!WRITE(*,*) 'FMort',rFMort(iSpp)
	
					IF(cDemPel(iSpp) =='D') THEN
						CALL CalK(.FALSE.,rBtT,rLinf(iSpp))
					ELSE
						CALL CalK(.FALSE.,rSST,rLinf(iSpp))
					END IF
	
					DO i=1,iTotCellNum
							!IF(rnewprop(i)/rprop(i)<=2) THEN
							IF(rnewprop(i)>0) THEN
								IF(rnewprop(i)<nppthreshold(iSpp)) THEN
													rK(i,iSpp)=rK(i,iSpp)*rnewprop(i)/nppthreshold(iSpp)
								END IF
							END IF
					END DO
	
					DO i=1,iTotCellNum
						IF(riniAbd(i,iSpp)==0) THEN
							IF((rK(i,iSpp)-rKOrig(i,iSpp))>(rKMin(iSpp)*0.5)) THEN
								rK(i,iSpp)=rK(i,iSpp)-rKOrig(i,iSpp)
							ELSE
								rK(i,iSpp)=0
							END IF
						END IF
					END DO
	
					!Calculate larval dispersal
					!!--------------- TRAVIS10 -------------!!
					!!--------------- Insert variables for pH effects on larval survival ---------------!!
					CALL CalLarDisp(rHQ10LarvMort(iSpp),rOrigHSurf,rNewHSurf)
					!!--------------- TRAVIS10 -------------!!
	
					!***********************************************************************
					!loop to create size bins
	
					IF(lFixedBin .EQV. .FALSE.) THEN
						!Create variable size bin depending on species
						rlencl(1)=(rLinf(iSpp) / (iSizeBinNum-10))/2
	
						DO itemp=2,iSizeBinNum
							rlencl(itemp) = rlencl(itemp - 1) + (rLinf(iSpp) /(iSizeBinNum-10))
						END DO
						rWgtcl=rLwA(iSpp)*rlencl**rLwB(iSpp)
					ELSE
						lLinfExceed=.FALSE.
						itemp=1
						rWgtcl(1)=1.0
						Do While ((itemp<=iSizeBinNum) .AND. (lLinfExceed .EQV. .FALSE.))
							itemp=itemp+1
							!Create a log2 series of weight bin, then convert to length bin
							rWgtcl(itemp)=2**(itemp-1)+(((2**(itemp))-(2**(itemp-1)))/2)
							rlencl(itemp)=(rWgtcl(itemp)/rLwA(iSpp))**(1/rLwB(iSpp))
							IF(rWgtcl(itemp)>=rWinf(iSpp)) THEN
								lLinfExceed=.FALSE.
							END IF
						END DO
					END IF
	
					!*********************************************************************************
	
					rMeanW(:)=0
	
					!********************************************************************************
					!Calculate biomass per recruit(rMeanW) per cell (and size-frequency distribution)
					!Use OpenMP for parallel computing
					CHUNK=CHUNKSize
					!$OMP	PARALLEL SHARED(rOrigMeanW) PRIVATE(i,TID)
	!				!$OMP	PARALLEL SHARED(rMeanW,rWinf,rVBonK,rLwA,rLwB, &
	!				!$OMP&		rArrhCoef,rAlloScal,rHQ10Met, &
	!				!$OMP&		rOrigHBtm,rSTDummy,rOrigO2btm,rKVar,rlencl,rWgtcl, &
	!				!$OMP&		pwater,cDemPel,NTHREADS,CHUNK,rNMort)  PRIVATE(i,rWinfDummy,rVBonKDummy, &
	!				!$OMP& 		rWmatDummy,rNMortDummy, rLinfDummy,iTmaxDummy,TID)
					TID = OMP_GET_THREAD_NUM()
					IF (TID == 0) THEN
						NTHREADS = OMP_GET_NUM_THREADS()
					!	PRINT *, 'Number of threads =', NTHREADS
					END IF
					!PRINT *, 'Thread',TID,' starting...'
					!$OMP DO SCHEDULE(GUIDED)
	
					DO i=1,iTotCellNum
						IF(rK(i,iSpp)>0) THEN
							IF(cDemPel(iSpp)=='D') THEN
	
							CALL Ecophysiology(iSizeBinNum,rWinf(iSpp),rVBonK(iSpp),rLwA(iSpp),rLwB(iSpp), &
							rArrhCoef(:,iSpp),rAlloScal(:,iSpp),rHQ10Met(iSpp),rOrigBtT(i),rOrigO2Btm(i),&
							rOrigHBtm(i),rBtT(i), rNewO2Btm(i),rNewHBtm(i),rKVar(iSpp),rlencl,rWgtcl, &
							rWinfDummy,rVBonKDummy,rWmatDummy,rNMortDummy,rLinfDummy,iTmaxDummy,rMeanW(i),rNMort(i), &
							!!------------ TRAVIS11 --------------!!
							!!------------ Include new variable rHQ10Mort(iSpp) ------------!!
							rHQ10Mort(iSpp))
							!!------------ TRAVIS11 --------------!!
	
							ELSE
	
							CALL Ecophysiology(iSizeBinNum,rWinf(iSpp),rVBonK(iSpp),rLwA(iSpp),rLwB(iSpp), &
							rArrhCoef(:,iSpp),rAlloScal(:,iSpp),rHQ10Met(iSpp),rOrigSST(i),rOrigO2Surf(i), &
							rOrigHSurf(i),rSST(i),rNewO2Surf(i),rNewHSurf(i),rKVar(iSpp),rlencl,rWgtcl, &
							rWinfDummy,rVBonKDummy,rWmatDummy,rNMortDummy,rLinfDummy,iTmaxDummy,rMeanW(i),rNMort(i), &
							!!------------ TRAVIS12 --------------!!
							!!------------ Include new variable rHQ10Mort(iSpp) ------------!!
							rHQ10Mort(iSpp))
							!!------------ TRAVIS12 --------------!!
	
							END IF
							!rNMort(i)=rNMortDummy(1)/inTimeStep
						END IF
	
					END DO
	
					!$OMP END DO
					!!NOWAIT
					!$OMP BARRIER
					!$OMP END PARALLEL
	
					rNMort(:)=rNMort(:)/inTimeStep
	
	
					!DO i=1,iTotCellNum
						!Fishing mortality rate
					!	rFMort(i,iSpp)=rFMortyr(i)*rIntR(iSpp)/2*(rHSMPAMap(i))
						!If the cell is MPA, set fishing mortality to zero
					!	IF(inTime>=MPAScStartYear)THEN
					!		IF(rHSMPAMap(i)==0) THEN !changed value to 0 if it is an MPA
					!			rFMort(i,iSpp)=0.0
					  !		END IF
					!	END IF
					!END DO
	
					!WRITE(*,*) 'KOrig= ',SUM(rK(:,iSpp))
					DO i=1,iTotCellNum
						IF (rOrigMeanW(i,iSpp)>0) THEN
							rK(i,iSpp)=rK(i,iSpp)*rMeanW(i)/rOrigMeanW(i,iSpp)
	
						ELSE
							rOrigMeanW(i,iSpp)=rMeanW(i)
						END IF
					END DO
	
					!*********************************************************************************
					!WRITE(*,*) SUM(rOrigMeanW(:,iSpp)), SUM(rMeanW(:))
	
	
					!WRITE(*,*)	'K=', SUM(rK(:,iSpp))
	
					IF (ISNAN(SUM(rK(:,iSpp)))) THEN
							fname='testrk2error.csv'
							CALL WriteCSV(fname,rK(:,iSpp),nrow*ncol)
	
							fname='testrWerror.csv'
							CALL WriteCSVDouble(fname,rMeanW(:),nrow*ncol)
	
							fname='testrNerror.csv'
							CALL WriteCSV(fname,rNMort(:),nrow*ncol)
	
					END IF
	
					WRITE(*,*) 'Abd=',SUM(rAbd)
	
					U=0
					V=0
	
					!recreating map segments by accounting for rK
					CALL MapSeg2
	
					rNMort(:)=0.0
					rGCoef(:)=1.0
					!Store maximum larval dis
					!DO i=1,nrow*ncol
					!	IF(rLarv(i)>rtLarvMax(iSpp)) THEN
					!		rGCoef(i)=1.0
					!	ELSE
					!		rGCoef(i)=0.0
					!	END IF
					!END DO
	
	
					DO iMonStep=1,inTimeStep
						dTCatchTemp(:)=0
						DO i=1,nrow*ncol
							rtCR=0
							IF(rK(i,iSpp)>0) THEN
								rtCR=MAX(0.0,(1-rAbd(i)/rK(i,iSpp)))
								rtCR=Min(1.0,rtCR)
								!Add BPR factor to determine change in diffcoef
								!rtDiffCoef(i)=1/(1+rDScaler*rHabSuit(i,iSpp)*rtCR*rMeanW(i)/rOrigMeanW(i,iSpp))
								!rtDiffCoef(i)=1/(1+rDScaler*rHabSuit(i,iSpp)*rtCR)
								rtDiffCoef(i)=2/(1+EXP(rDScaler*rHabSuit(i,iSpp)*rtCR))
								rtDiffCoef(i)=rDiffCoef(iSpp)*rtDiffCoef(i)
							ELSE
								rtDiffCoef(i)=rDiffCoef(iSpp)
							END IF
						END DO
						!testing
						!WRITE(*,*)	rAbd(19837)/rK(19837,iSpp),rHabSuit(19837,iSpp),rtDiffCoef(19837)
	
						CALL setabcdef
	
						CALL caltridag
						!rNMort(:)=0.0
						dTCatchTemp(:)=0.0
						CALL CalDeltaAbd(rFMort(:,iSpp),rNMort,rGCoef,dTCatchTemp)
						dTCatch(:)=dTCatchTemp(:)
	
						IF (SUM(rAbd)==0) THEN
							EXIT
						END IF
					End DO
					WRITE(*,*)	'Catch=', SUM(dTCatch(:))
					DO i= 1,iTotCellNum
						IF(rAbd(i)==0 .AND. rLarv(i)>0 .AND. rK(i,iSpp)>0) THEN
							rAbd(i)=1
						END IF
					END DO
					!fname='testDCf.txt'
					!CALL WriteCSV(fname,rtDiffCoef,nrow*ncol)
	
					rAbdMaster(:,iSpp)=rAbd
					Write(cTemp2, '(i6)') iTaxonList(iSpp)
					fname=TRIM(cTemp2)//'Abd'//TRIM(ctemp)//'.txt'
					CALL WriteCSV(fname,rAbd/rMCorrect(iSpp),nrow*ncol) !Travis
	
					Write(cTemp2, '(i6)') iTaxonList(iSpp)
					fname=TRIM(cTemp2)//'Catch'//TRIM(ctemp)//'.txt'
					CALL WriteCSVDouble(fname,dTCatch/rMCorrect(iSpp),nrow*ncol)
	
					!fname='MeanW'//TRIM(ctemp)//'.csv'
					!CALL WriteCSVDouble(fname,rMeanW,nrow*ncol)
	
					rAbd=0
					rOrigMeanW(:,iSpp)=rMeanW(:)
				END IF
			END DO
		END DO
	
		WRITE(*,*)	'Completed!'
	
	END PROGRAM ADR
	
	SUBROUTINE MapSeg
	
		USE DBEM, ONLY :	pwater, iHorsegncell,iVersegncell,iHorcellref,iVercellref, &
		Lat,nrow,ncol,rA,rB,rC,rD,rE,rF,rG,rH
		IMPLICIT NONE
		INTEGER	:: icolcnt,irowcnt,icolcnt2,i,isegcnt,iseq,icellcnt,iMaxSegLen
		LOGICAL	:: lsegstart
	
		!Check if variables are allocated, if yes, deallocate first.
		IF(ALLOCATED(rA) .EQV. .TRUE.) THEN
	
			DEALLOCATE (iHorcellref)
			DEALLOCATE (iHorsegncell)
			DEALLOCATE (rA)
			DEALLOCATE (rB)
			DEALLOCATE (rC)
			DEALLOCATE (rH)
			DEALLOCATE (iVercellref)
			DEALLOCATE (iVersegncell)
			DEALLOCATE (rD)
			DEALLOCATE (rE)
			DEALLOCATE (rF)
			DEALLOCATE (rG)
		END IF
	
		!This loop to create horizontal segment
		iseq=0	!set sequence number counter to 0
		isegcnt=0	!set segment counter to 0
		icellcnt=0
		lsegstart=.FALSE. !set start segment pointer as false
		iMaxSegLen=0
		DO irowcnt=1,nrow
			DO icolcnt=1,ncol !loop through all column
				iseq=(irowcnt-1)*ncol+icolcnt !calculate sequence number
				IF (pwater(iseq)>0 ) THEN
					IF(lsegstart .EQV. .FALSE.) THEN
						lsegstart = .TRUE.
						isegcnt=isegcnt+1
						icellcnt=1
					ELSE
						icellcnt=icellcnt+1
					END IF
				END IF
				IF(icolcnt/=ncol .AND. pwater((irowcnt-1)*ncol+icolcnt+1)<=0) THEN
					lsegstart = .FALSE.
					iMaxSegLen=MAX(iMaxSegLen,icellcnt)
				END IF
				!Spherical horizontal segment by ending the eastern most segment of the map
				IF (icolcnt==ncol .AND. pwater(iseq)>0) THEN
					icolcnt2=0
					DO WHILE (pwater((irowcnt-1)*ncol+icolcnt+1)>0 .AND. icolcnt2<=100)
						icolcnt2=icolcnt2+1
						icellcnt=icellcnt+1
					END DO
				END IF
			END DO
			lsegstart = .FALSE.
			IF (pwater(iseq)>0) THEN
				iMaxSegLen=MAX(iMaxSegLen,icellcnt)
			END IF
		END DO
		!Re-dimension segment pointers
		ALLOCATE (iHorcellref(isegcnt,iMaxSegLen))
		ALLOCATE (iHorsegncell(isegcnt))
		ALLOCATE (rA(isegcnt,iMaxSegLen))
		ALLOCATE (rB(isegcnt,iMaxSegLen))
		ALLOCATE (rC(isegcnt,iMaxSegLen))
		ALLOCATE (rH(isegcnt,iMaxSegLen))
	
		!This loop to create vertical segment
		iseq=0	!set sequence number counter to 0
		isegcnt=0	!set segment counter to 0
		icellcnt=0
		lsegstart=.FALSE. !set start segment pointer as false
		iMaxSegLen=0
		DO icolcnt=1,ncol
			DO irowcnt=1,nrow
				iseq=(irowcnt-1)*ncol+icolcnt
				IF (pwater(iseq)>0 ) THEN
					IF(lsegstart .EQV. .FALSE.) THEN
						lsegstart = .TRUE.
						isegcnt=isegcnt+1
						icellcnt=1
					ELSE
						icellcnt=icellcnt+1
					END IF
				END IF
				IF(irowcnt/=nrow.AND. pwater((irowcnt-1+1)*ncol+icolcnt)<=0) THEN
					lsegstart = .FALSE.
					iMaxSegLen=MAX(iMaxSegLen,icellcnt)
				END IF
			END DO
			lsegstart = .FALSE.
			IF (pwater(iseq)>0) THEN
				iMaxSegLen=MAX(iMaxSegLen,icellcnt)
			END IF
		END DO
	
		!Re-dimension segment pointers
		ALLOCATE (iVercellref(isegcnt,iMaxSegLen))
		ALLOCATE (iVersegncell(isegcnt))
		ALLOCATE (rD(isegcnt,iMaxSegLen))
		ALLOCATE (rE(isegcnt,iMaxSegLen))
		ALLOCATE (rF(isegcnt,iMaxSegLen))
		ALLOCATE (rG(isegcnt,iMaxSegLen))
	
		!WRITE(*,*) isegcnt,iMaxSegLen !for testing
	
		!Repeat the horizontal and vertical loops, but store value this time.
		!This loop to create horizontal segment
		iseq=0	!set sequence number counter to 0
		isegcnt=0	!set segment counter to 0
		icellcnt=0
		lsegstart=.FALSE. !set start segment pointer as false
		DO irowcnt=1,nrow
			DO icolcnt=1,ncol !loop through all column
				iseq=(irowcnt-1)*ncol+icolcnt !calculate sequence number
				IF (pwater(iseq)>0 ) THEN
					IF(lsegstart .EQV. .FALSE.) THEN
						lsegstart = .TRUE.
						isegcnt=isegcnt+1
						icellcnt=1
					ELSE
						icellcnt=icellcnt+1
					END IF
					iHorCellref(isegcnt,icellcnt)=iseq
				END IF
				IF(icolcnt/=ncol .AND. pwater((irowcnt-1)*ncol+icolcnt+1)<=0) THEN
					lsegstart = .FALSE.
					iHorsegncell(isegcnt)=icellcnt
				END IF
				!Spherical horizontal segment by ending the eastern most segment of the map
				IF (icolcnt==ncol .AND. pwater(iseq)>0) THEN
					icolcnt2=1
					DO WHILE (pwater((irowcnt-1)*ncol+icolcnt+icolcnt2)>0 .AND. icolcnt2<=100)
						icellcnt=icellcnt+1
						iHorCellref(isegcnt,icellcnt)=(irowcnt-1)*ncol+icolcnt2
						icolcnt2=icolcnt2+1
					END DO
					iHorsegncell(isegcnt)=icellcnt
				END IF
			END DO
			lsegstart = .FALSE.
		END DO
	
	
	
		!This loop to create vertical segment
		iseq=0	!set sequence number counter to 0
		isegcnt=0	!set segment counter to 0
		icellcnt=0
		lsegstart=.FALSE. !set start segment pointer as false
		DO icolcnt=1,ncol
			DO irowcnt=1,nrow
				iseq=(irowcnt-1)*ncol+icolcnt
				IF (pwater(iseq)>0 ) THEN
					IF(lsegstart .EQV. .FALSE.) THEN
						lsegstart = .TRUE.
						isegcnt=isegcnt+1
						icellcnt=1
					ELSE
						icellcnt=icellcnt+1
					END IF
					iVerCellref(isegcnt,icellcnt)=iseq
				END IF
				IF(irowcnt/=nrow.AND. pwater((irowcnt-1+1)*ncol+icolcnt)<=0) THEN
					lsegstart = .FALSE.
					iVersegncell(isegcnt)=icellcnt
				END IF
			END DO
			if (pwater(iseq)>0) THEN
				iVersegncell(isegcnt)=icellcnt
			END IF
			lsegstart = .FALSE.
		END DO
	
		!WRITE(*,*)	SIZE(iVerCellref,1), SIZE(iVerCellref,2) !testing
		!WRITE(*,*)	SIZE(iHorCellref,1), SIZE(iHorCellref,2) !testing
		!WRITE(*,*)	SIZE(iVersegncell,1) !testing
		!WRITE(*,*)	SIZE(iHorsegncell,1) !testing
	END SUBROUTINE MapSeg
	
	SUBROUTINE MapSeg2
	
		USE DBEM, ONLY :	pwater, iHorsegncell,iVersegncell,iHorcellref,iVercellref, &
		Lat,nrow,ncol,rA,rB,rC,rD,rE,rF,rG,rH,rK,iSpp
		IMPLICIT NONE
		INTEGER	:: icolcnt,irowcnt,icolcnt2,i,isegcnt,iseq,icellcnt,iMaxSegLen
		LOGICAL	:: lsegstart
	
		DEALLOCATE (iHorcellref)
		DEALLOCATE (iHorsegncell)
		DEALLOCATE (rA)
		DEALLOCATE (rB)
		DEALLOCATE (rC)
		DEALLOCATE (rH)
		DEALLOCATE (iVercellref)
		DEALLOCATE (iVersegncell)
		DEALLOCATE (rD)
		DEALLOCATE (rE)
		DEALLOCATE (rF)
		DEALLOCATE (rG)
	
		!This loop to create horizontal segment
		iseq=0	!set sequence number counter to 0
		isegcnt=0	!set segment counter to 0
		icellcnt=0
		lsegstart=.FALSE. !set start segment pointer as false
		iMaxSegLen=0
		DO irowcnt=1,nrow
			DO icolcnt=1,ncol !loop through all column
				iseq=(irowcnt-1)*ncol+icolcnt !calculate sequence number
				IF (pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
					IF(lsegstart .EQV. .FALSE.) THEN
						lsegstart = .TRUE.
						isegcnt=isegcnt+1
						icellcnt=1
					ELSE
						icellcnt=icellcnt+1
					END IF
				END IF
				IF(icolcnt/=ncol .AND. (pwater((irowcnt-1)*ncol+icolcnt+1)<=0 &
				.OR. rK((irowcnt-1)*ncol+icolcnt+1,iSpp)==0)) THEN
					lsegstart = .FALSE.
					iMaxSegLen=MAX(iMaxSegLen,icellcnt)
				END IF
				!Spherical horizontal segment by ending the eastern most segment of the map
				IF (icolcnt==ncol .AND. pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
					icolcnt2=0
					DO WHILE (pwater((irowcnt-1)*ncol+icolcnt+1)>0 .AND. icolcnt2<=100 &
					.AND. rK((irowcnt-1)*ncol+icolcnt+1,iSpp)>0)
						icolcnt2=icolcnt2+1
						icellcnt=icellcnt+1
					END DO
				END IF
			END DO
			lsegstart = .FALSE.
			IF (pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
				iMaxSegLen=MAX(iMaxSegLen,icellcnt)
			END IF
		END DO
		!Re-dimension segment pointers
		ALLOCATE (iHorcellref(isegcnt,iMaxSegLen))
		ALLOCATE (iHorsegncell(isegcnt))
		ALLOCATE (rA(isegcnt,iMaxSegLen))
		ALLOCATE (rB(isegcnt,iMaxSegLen))
		ALLOCATE (rC(isegcnt,iMaxSegLen))
		ALLOCATE (rH(isegcnt,iMaxSegLen))
	
	
		!This loop to create vertical segment
		iseq=0	!set sequence number counter to 0
		isegcnt=0	!set segment counter to 0
		icellcnt=0
		lsegstart=.FALSE. !set start segment pointer as false
		iMaxSegLen=0
		DO icolcnt=1,ncol
			DO irowcnt=1,nrow
				iseq=(irowcnt-1)*ncol+icolcnt
				IF (pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
					IF(lsegstart .EQV. .FALSE.) THEN
						lsegstart = .TRUE.
						isegcnt=isegcnt+1
						icellcnt=1
					ELSE
						icellcnt=icellcnt+1
					END IF
				END IF
				IF(irowcnt/=nrow .AND. (lsegstart .EQV. .TRUE.) .AND. (pwater((irowcnt-1+1)*ncol+icolcnt)<=0 &
				.OR. rK((irowcnt-1+1)*ncol+icolcnt,iSpp)==0)) THEN
					lsegstart = .FALSE.
					iMaxSegLen=MAX(iMaxSegLen,icellcnt)
				END IF
			END DO
			lsegstart = .FALSE.
			IF (pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
				iMaxSegLen=MAX(iMaxSegLen,icellcnt)
			END IF
		END DO
	
		!Re-dimension segment pointers
		ALLOCATE (iVercellref(isegcnt,iMaxSegLen))
		ALLOCATE (iVersegncell(isegcnt))
		ALLOCATE (rD(isegcnt,iMaxSegLen))
		ALLOCATE (rE(isegcnt,iMaxSegLen))
		ALLOCATE (rF(isegcnt,iMaxSegLen))
		ALLOCATE (rG(isegcnt,iMaxSegLen))
	
	
		!WRITE(*,*) isegcnt,iMaxSegLen !for testing
	
		!Repeat the horizontal and vertical loops, but store value this time.
		!This loop to create horizontal segment
		iseq=0	!set sequence number counter to 0
		isegcnt=0	!set segment counter to 0
		icellcnt=0
		lsegstart=.FALSE. !set start segment pointer as false
		DO irowcnt=1,nrow
			DO icolcnt=1,ncol !loop through all column
				iseq=(irowcnt-1)*ncol+icolcnt !calculate sequence number
				IF (pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
					IF(lsegstart .EQV. .FALSE.) THEN
						lsegstart = .TRUE.
						isegcnt=isegcnt+1
						icellcnt=1
					ELSE
						icellcnt=icellcnt+1
					END IF
					iHorCellref(isegcnt,icellcnt)=iseq
				END IF
				IF(icolcnt/=ncol .AND. (lsegstart .EQV. .TRUE.) .AND. (pwater((irowcnt-1)*ncol+icolcnt+1)<=0 &
				.OR. rK((irowcnt-1)*ncol+icolcnt+1,iSpp)==0)) THEN
					lsegstart = .FALSE.
					iHorsegncell(isegcnt)=icellcnt
				END IF
				!Spherical horizontal segment by ending the eastern most segment of the map
				IF (icolcnt==ncol .AND. pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
					icolcnt2=1
					DO WHILE (pwater((irowcnt-1)*ncol+icolcnt+icolcnt2)>0 .AND. icolcnt2<=100 &
					.AND. rK((irowcnt-1)*ncol+icolcnt+icolcnt2,iSpp)>0)
						icellcnt=icellcnt+1
						iHorCellref(isegcnt,icellcnt)=(irowcnt-1)*ncol+icolcnt2
						icolcnt2=icolcnt2+1
					END DO
					iHorsegncell(isegcnt)=icellcnt
				END IF
			END DO
			lsegstart = .FALSE.
		END DO
	
	
		!This loop to create vertical segment
		iseq=0	!set sequence number counter to 0
		isegcnt=0	!set segment counter to 0
		icellcnt=0
		lsegstart=.FALSE. !set start segment pointer as false
		DO icolcnt=1,ncol
			DO irowcnt=1,nrow
				iseq=(irowcnt-1)*ncol+icolcnt
				IF (pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
					IF(lsegstart .EQV. .FALSE.) THEN
						lsegstart = .TRUE.
						isegcnt=isegcnt+1
						icellcnt=1
					ELSE
						icellcnt=icellcnt+1
					END IF
					iVerCellref(isegcnt,icellcnt)=iseq
				END IF
	
				IF(irowcnt/=nrow .AND. (lsegstart .EQV. .TRUE.) .AND. (pwater((irowcnt-1+1)*ncol+icolcnt)<=0 &
				.OR. rK((irowcnt-1+1)*ncol+icolcnt,iSpp)==0)) THEN
					lsegstart = .FALSE.
					iVersegncell(isegcnt)=icellcnt
				END IF
			END DO
			if (pwater(iseq)>0 .AND. rK(iseq,iSpp)>0) THEN
				iVersegncell(isegcnt)=icellcnt
			END IF
			lsegstart = .FALSE.
		END DO
	
		!WRITE(*,*)	SIZE(iVerCellref,1), SIZE(iVerCellref,2) !testing
		!WRITE(*,*)	SIZE(iHorCellref,1), SIZE(iHorCellref,2) !testing
		!WRITE(*,*)	SIZE(iVersegncell,1) !testing
		!WRITE(*,*)	SIZE(iHorsegncell,1) !testing
	END SUBROUTINE MapSeg2
	
	SUBROUTINE setabcdef
	
	!*************************************************************************************************************
	!Subroutine to set the coefficients (a, b, c, d, e, f) of the advection-diffusion equation solution by Sibert
	!Fortran version
	!Created by: William Cheung
	!Date: 23 Dec 2013
	!*************************************************************************************************************
	
		USE DBEM, ONLY : U,V,iHorsegncell,iVersegncell,iHorcellref,iVercellref, &
		nrow,ncol,rtDiffCoef,iDeltaT,rDX0,rDY0,rMort,rA,rB,rC,rD,rE,rF
		IMPLICIT NONE
		INTEGER	:: icnt, isegcnt,islen,itempseq,i,itempseqm1,itempseqp1
		REAL	:: rUm1,rU0,rUp1,rtDX0,rtDXm1,rtDXp1,rVm1,rV0,rVp1,rtDY0,rtDYm1,rtDYp1
	
		!Horizontal segment
		DO isegcnt=1,SIZE(iHorsegncell)
			islen=iHorsegncell(isegcnt) !no of cells in the segment
	
			DO icnt=1,islen
				itempseq=iHorcellref(isegcnt,icnt) !sequence number of the segment cell
	
	
	
				!Extract Current U value
				rU0=U(itempseq)
	
				!Extract distance
				rtDX0=rDX0(itempseq)
	
				IF (icnt>1) THEN
					rUm1=U(iHorcellref(isegcnt,icnt-1))
					rtDXm1=rDX0(iHorcellref(isegcnt,icnt-1))
					itempseqm1=iHorcellref(isegcnt,icnt-1) !sequence number of the segment cell (-1)
				END IF
				IF (icnt<islen) THEN
					rUp1=U(iHorcellref(isegcnt,icnt+1))
					rtDXp1=rDX0(iHorcellref(isegcnt,icnt+1))
					itempseqp1=iHorcellref(isegcnt,icnt+1) !sequence number of the segment cell (+1)
				END IF
	
				IF (icnt==1 .AND. islen>1) THEN !if i=1
					rA(isegcnt,icnt)=0
					rB(isegcnt,icnt)=2/iDeltaT+rtDiffCoef(itempseq)/(rtDX0**2)+rMort
					IF (icnt/=islen) THEN !if i is not end of segment
						IF(rUp1>=0) THEN
							rC(isegcnt,icnt)=(-(rtDiffCoef(itempseqp1)/(rtDXp1**2)))
						ELSE
							rC(isegcnt,icnt)=-((rtDiffCoef(itempseqp1)/(rtDXp1**2))-rUp1/rtDXp1)
						END IF
					END IF
				ELSE IF (icnt==1 .AND. islen==1) THEN
					rA(isegcnt,icnt)=0
					rB(isegcnt,icnt)=rMort
					rC(isegcnt,icnt)=0
				ELSE IF (icnt>1 .AND. icnt/=islen) THEN !if i>1 and < end of segment
					IF (rUm1>0) THEN
						rA(isegcnt,icnt)=(-(rtDiffCoef(itempseqm1)/(rtDXm1**2)+rUm1/rtDXm1))
					ELSE
						rA(isegcnt,icnt)=(-rtDiffCoef(itempseqm1)/(rtDXm1**2))
					END IF
					IF (rU0>0) THEN
						rB(isegcnt,icnt)=2/iDeltaT+2*rtDiffCoef(itempseq)/(rtDX0**2)+rU0/rtDx0+rMort
					ELSE
						rB(isegcnt,icnt)=2/iDeltaT+2*rtDiffCoef(itempseq)/(rtDx0**2)-rU0/rtDX0+rMort
					END IF
					IF (rUp1>=0) THEN
						rC(isegcnt,icnt)=(-(rtDiffCoef(itempseqp1)/(rtDXp1**2)))
					ELSE
						rC(isegcnt,icnt)=(-(rtDiffCoef(itempseqp1)/(rtDXp1**2)-rUp1/rtDXp1))
					END IF
	
				ELSE IF (icnt==islen) THEN
					IF (rUm1>0) THEN
						rA(isegcnt,icnt)=(-(rtDiffCoef(itempseqm1)/(rtDXm1**2)+rUm1/rtDXm1))
					ELSE
						rA(isegcnt,icnt)=(-rtDiffCoef(itempseqm1)/(rtDXm1**2))
					END IF
					rB(isegcnt,icnt)=2/iDeltaT+rtDiffCoef(itempseq)/(rtDX0**2)+rMort
					rC(isegcnt,icnt)=0
	
				END IF
			!	IF (itempseq==74766) THEN
			!		WRITE (*,*) rA(isegcnt,icnt),rB(isegcnt,icnt),rC(isegcnt,icnt), itempseq
			!	END IF
			END DO
		END DO
	
		!WRITE (*,*) SIZE(rA,1),SIZE(rA,2),SIZE(iHorsegncell),iHorsegncell(1) !testing
		!Vertical segment
		DO isegcnt=1,SIZE(iVersegncell)
			islen=iVersegncell(isegcnt)
			DO icnt=1,islen
				itempseq=iVercellref(isegcnt,icnt) !sequence number of the segment cell
	
				!Extract Current V value
				rV0=V(itempseq)
	
				!Extract distance
				rtDY0=rDY0(itempseq)
	
				IF (icnt>1) THEN
					rVm1=V(iVercellref(isegcnt,icnt-1))
					rtDYm1=rDY0(iVercellref(isegcnt,icnt-1))
					itempseqm1=iVercellref(isegcnt,icnt-1) !sequence number of the segment cell
	
				END IF
				IF(icnt<islen) THEN
					rVp1=V(iVercellref(isegcnt,icnt+1))
					rtDYp1=rDY0(iVercellref(isegcnt,icnt+1))
					itempseqp1=iVercellref(isegcnt,icnt+1) !sequence number of the segment cell
	
				END IF
	
				IF (icnt==1 .AND. islen>1) THEN
					rD(isegcnt,icnt)=0
					rE(isegcnt,icnt)=2/iDeltaT+rtDiffCoef(itempseq)/(rtDY0**2)
					IF (islen>0) THEN
						IF (rVp1>=0) THEN
							rF(isegcnt,icnt)=(-rtDiffCoef(itempseqp1)/(rtDYp1**2))
						ELSE
							rF(isegcnt,icnt)=(-(rtDiffCoef(itempseqp1)/(rtDYp1**2)-rVp1/rtDYp1))
						END IF
					END IF
				ELSE IF (icnt==1 .AND. islen==1) THEN
					rD(isegcnt,icnt)=0
					rE(isegcnt,icnt)=1
					rF(isegcnt,icnt)=0
				ELSE IF (icnt>1 .AND. icnt<islen) THEN
					IF (rVm1>=0) THEN
						rD(isegcnt,icnt)=(-(rtDiffCoef(itempseqm1)/(rtDYm1**2)+rVm1/rtDYm1))
					ELSE
						rD(isegcnt,icnt)=(-(rtDiffCoef(itempseqm1)/(rtDYm1**2)))
					END IF
					IF (rV0>=0) THEN
						rE(isegcnt,icnt)=2/iDeltaT+2*rtDiffCoef(itempseq)/(rtDY0**2)+rV0/rtDY0
					ELSE
						rE(isegcnt,icnt)=2/iDeltaT+2*rtDiffCoef(itempseq)/(rtDY0**2)-rV0/rtDY0
					END IF
					IF (rVp1>=0) THEN
						rF(isegcnt,icnt)=(-rtDiffCoef(itempseqp1)/(rtDYp1**2))
					ELSE
						rF(isegcnt,icnt)=(-(rtDiffCoef(itempseqp1)/(rtDYp1**2)-rVp1/rtDYp1))
					END IF
				ELSE IF (icnt==islen .AND. islen>1) THEN
					IF (rVm1>=0) THEN
						rD(isegcnt,icnt)=(-(rtDiffCoef(itempseqm1)/(rtDYm1**2)+rVm1/rtDYm1))
					ELSE
						rD(isegcnt,icnt)=(-(rtDiffCoef(itempseqm1)/(rtDYm1**2)))
					END IF
					rE(isegcnt,icnt)=2/iDeltaT+2*rtDiffCoef(itempseq)/(rtDY0**2)
					rF(isegcnt,icnt)=0
				END IF
	
			END DO
		END DO
	
	END SUBROUTINE setabcdef
	
	SUBROUTINE setg
		!Subroutine to calculate the g parameter
		USE DBEM, ONLY : rD, rE, rF,rG,iVersegncell,iVercellref,iDeltaT,rAbd,rtDiffCoef
		IMPLICIT NONE
		INTEGER	:: icnt,islen,isegcnt,itempseq0,itempseqM1,itempseqP1
		REAL	:: rTime
		rTime=REAL(iDeltaT)/2
	
		DO isegcnt=1,SIZE(iVersegncell)
			islen=iVersegncell(isegcnt)
			DO icnt=1,islen
				!sequence number of the segment cell
				itempseq0=iVercellref(isegcnt,icnt)
				IF (icnt>1) THEN
					itempseqM1=iVercellref(isegcnt,icnt-1)
				END IF
				IF (icnt<islen) THEN
					itempseqP1=iVercellref(isegcnt,icnt+1)
				END IF
				!calculate g parameter for tridiagonal matrix
				!condition on different position on the segment
				IF (icnt==1) THEN
					IF (islen==1) THEN
						rG(isegcnt,icnt)=rAbd(itempseq0)
					ELSE IF (icnt<islen) THEN
						rG(isegcnt,icnt)=(2/rTime-rE(isegcnt,icnt))*rAbd(itempseq0)
						rG(isegcnt,icnt)=rG(isegcnt,icnt)-rF(isegcnt,icnt)*rAbd(itempseqP1)
					END IF
				ELSE IF (icnt>1 .AND. icnt<islen) THEN
					rG(isegcnt,icnt)=(-rD(isegcnt,icnt)*rAbd(itempseqM1)+(2/rTime-rE(isegcnt,icnt))*rAbd(itempseq0))
					rG(isegcnt,icnt)=rG(isegcnt,icnt)-rF(isegcnt,icnt)*rAbd(itempseqP1)
				ELSE IF (icnt==islen .AND. islen>1) THEN
					rG(isegcnt,icnt)=-rD(isegcnt,icnt)*rAbd(itempseqM1)+(2/rTime-rE(isegcnt,icnt))*rAbd(itempseq0)
				END IF
			END DO
		END DO
	
	END SUBROUTINE setg
	
	SUBROUTINE seth
		!calculate h parameter
		USE DBEM, ONLY : rA,rB,rC,rH,iHorsegncell,iHorcellref,iDeltaT,rAbd,rtDiffCoef
		IMPLICIT NONE
		INTEGER	:: icnt,islen,isegcnt,itempseq0,itempseqM1,itempseqP1
		REAL	:: rTime
		rTime=REAL(iDeltaT)/2
		DO isegcnt=1,SIZE(iHorsegncell)
			islen=iHorsegncell(isegcnt)
			DO icnt=1,islen
				!sequence number of the segment cell
				itempseq0=iHorcellref(isegcnt,icnt)
				IF (icnt>1) THEN
					itempseqM1=iHorcellref(isegcnt,icnt-1)
				END IF
				IF (icnt<islen) THEN
					itempseqP1=iHorcellref(isegcnt,icnt+1)
				END IF
				!calculate h parameter for tridiagonal matrix
				!condition on different position on the segment
				IF (icnt==1) THEN
					IF (islen==1) THEN
						rH(isegcnt,icnt)=rAbd(itempseq0)
					ELSE
						rH(isegcnt,icnt)=(2/rTime-rB(isegcnt,icnt))*rAbd(itempseq0)-rC(isegcnt,icnt)*rAbd(itempseqP1)
					END IF
				ELSE IF (icnt>1 .AND. icnt<islen) THEN
					rH(isegcnt,icnt)=(-rA(isegcnt,icnt)*rAbd(itempseqM1)+(2/rTime-rB(isegcnt,icnt))*rAbd(itempseq0))
					rH(isegcnt,icnt)=rH(isegcnt,icnt)-rC(isegcnt,icnt)*rAbd(itempseqP1)
				ELSE IF (icnt==islen .AND. islen>1) THEN
					rH(isegcnt,icnt)=-rA(isegcnt,icnt)*rAbd(itempseqM1)+(2/rTime-rB(isegcnt,icnt))*rAbd(itempseq0)
				END IF
	
			END DO
		END DO
		RETURN
	END SUBROUTINE seth
	
	
	SUBROUTINE caltridag
		!Subroutine to calculate dispersal using the tridiagonal matrix
		!To loop through all segments and then call tridag
		USE DBEM, ONLY	: rA,rB,rC,rD,rE,rF,rG,rH,iHorsegncell,iHorcellref,iVersegncell, &
			iVercellref,rAbd,nrow,ncol
		IMPLICIT NONE
		INTEGER	::	icnt,islen,isegcnt,itempseq
		INTEGER,ALLOCATABLE	:: itseq(:)
		REAL,ALLOCATABLE	::	rtA(:),rtB(:),rtC(:),rtD(:),rtE(:),rtF(:),rtG(:),rtH(:),rtempAbd(:)
		LOGICAL	:: ltemp
		CHARACTER(20)	:: fname
		REAL,DIMENSION(nrow*ncol)	:: rdummyAbd
		REAL	:: rtotaltemp1,rtotaltemp2
	
		!Record total abundance
		rtotaltemp1=0
		DO icnt=1,nrow*ncol
			rtotaltemp1=rtotaltemp1+rAbd(icnt)
		END DO
		!Solve vertical segments
		!testing
		!OPEN (UNIT=1,file='/Users/wwlcheung/Fortran/Data/testvarT.csv',ACTION="write")
		!OPEN (UNIT=2,file='/Users/wwlcheung/Fortran/Data/testseg1.csv',ACTION="write")
		CALL Setg
	
		DO isegcnt=1,SIZE(iVersegncell)
			islen=iVersegncell(isegcnt)
			ALLOCATE (itseq(islen))
			ALLOCATE (rtD(islen))
			ALLOCATE (rtE(islen))
			ALLOCATE (rtF(islen))
			ALLOCATE (rtG(islen))
			ALLOCATE (rtempAbd(islen))
			DO icnt=1,islen
				!get sequence number
				itempseq=iVercellref(isegcnt,icnt)
				itseq(icnt)=itempseq
				rtD(icnt)=rD(isegcnt,icnt)
				rtE(icnt)=rE(isegcnt,icnt)
				rtF(icnt)=rF(isegcnt,icnt)
				rtG(icnt)=rG(isegcnt,icnt)
				rtempAbd(icnt)=rAbd(itempseq)
				!testseq(itempseq)=icnt
			END DO
			ltemp=tridag(rtD,rtE,rtF,rtG,islen,rtempAbd,itseq)
	
			!Write in abundance values
			DO icnt=1,islen
				itempseq=iVercellref(isegcnt,icnt)
				rAbd(itempseq)=rtempAbd(icnt)
	
				!IF(Lat(itempseq)==27.75 .AND. Lon(itempseq)==124.25) THEN
				!	WRITE(*,*)	rD(isegcnt,icnt),rE(isegcnt,icnt),rF(isegcnt,icnt),rG(isegcnt,icnt),rAbd(itempseq),isegcnt,icnt, &
				!	Lat(itempseq),Lon(itempseq)
				!END IF
			END DO
	
			DEALLOCATE (itseq)
			DEALLOCATE (rtD)
			DEALLOCATE (rtE)
			DEALLOCATE (rtF)
			DEALLOCATE (rtG)
			DEALLOCATE (rtempAbd)
		END DO
		!CLOSE(1)
		!testing
		!fname='test1.csv'
		!CALL WriteCSV(fname,rAbd,259200)
	
		!Solve Horizontal segments
		rdummyAbd=0
		Call seth
		!OPEN(UNIT=1,FILE='/Users/wwlcheung/Fortran/Data/testseg.csv')
		DO isegcnt=1,SIZE(iHorsegncell)
			islen=iHorsegncell(isegcnt)
			ALLOCATE (itseq(islen))
			ALLOCATE (rtA(islen))
			ALLOCATE (rtB(islen))
			ALLOCATE (rtC(islen))
			ALLOCATE (rtH(islen))
			ALLOCATE (rtempAbd(islen)) !temp variable for Abd
			DO icnt=1,islen
				!get sequence number
				itempseq=iHorcellref(isegcnt,icnt)
				itseq(icnt)=itempseq
				rtA(icnt)=rA(isegcnt,icnt)
				rtB(icnt)=rB(isegcnt,icnt)
				rtC(icnt)=rC(isegcnt,icnt)
				rtH(icnt)=rH(isegcnt,icnt)
				rtempAbd(icnt)=rAbd(itempseq)
		!		IF(Lat(itempseq)==27.75 .AND. Lon(itempseq)==124.25) THEN
		!			WRITE(*,*)	rAbd(itempseq),isegcnt,icnt, &
		!			Lat(itempseq),Lon(itempseq)
		!		END IF
			END DO
			ltemp=tridag(rtA,rtB,rtC,rtH,islen,rtempAbd,itseq)
	
			DO icnt=1,islen
				itempseq=iHorcellref(isegcnt,icnt)
				rdummyAbd(itempseq)=rtempabd(icnt) !storing Abd in dummy variable to update spherical segment
				!IF (isegcnt==997) THEN
				!	WRITE(1,*)	rA(isegcnt,icnt),rB(isegcnt,icnt),rC(isegcnt,icnt),rH(isegcnt,icnt),rAbd(itempseq), &
				!	rdummyAbd(itempseq),isegcnt,icnt, &
				!	Lat(itempseq),Lon(itempseq)
				!END IF
			END DO
	
			DEALLOCATE (itseq)
			DEALLOCATE (rtA)
			DEALLOCATE (rtB)
			DEALLOCATE (rtC)
			DEALLOCATE (rtH)
			DEALLOCATE (rtempAbd)
		END DO
	
		!CLOSE(1)
		!update final Abd
		rtotaltemp2=0
		DO icnt=1,nrow*ncol
			IF(rdummyabd(icnt)<1) THEN
				rdummyabd(icnt)=0
			END IF
			rtotaltemp2=rtotaltemp2+rdummyabd(icnt)
		END DO
		rAbd=rdummyAbd/rtotaltemp2*rtotaltemp1
		!testing
		!fname='test2.csv'
		!CALL WriteCSV(fname,rAbd,259200)
		CONTAINS
	
		LOGICAL FUNCTION tridag(rt1,rt2,rt3,rtR,itn,TriU,itseq)
			!solution to tridiagonal matrix
			!Based on Presse et al. 1988 Numerical Recipes for Fortran
	
			IMPLICIT NONE
			REAL,DIMENSION(itn),INTENT(IN)	:: rt1,rt2,rt3,rtR
			INTEGER,DIMENSION(*),INTENT(IN)	::	itseq
			INTEGER,INTENT(IN)	::	itn
			REAL,DIMENSION(*),INTENT(OUT)	:: TriU
			INTEGER	:: icnt
			REAL	:: bet,gam(itn)
	
			!if segment only have one cell then
			IF (itn==1) THEN
				TriU(1)=rt2(1)*rtR(1)
			ELSE IF (itn>1) THEN
				bet=rt2(1)
				TriU(1)=rtR(1)/bet
				DO icnt=2,itn
					gam(icnt)=rt3(icnt-1)/bet
					bet=rt2(icnt)-rt1(icnt)*gam(icnt)
					IF(bet==0) THEN
						WRITE(*,*)	'Error in tridag: Division by 0'
					END IF
					TriU(icnt)=(rtR(icnt)-rt1(icnt)*TriU(icnt-1))/bet
				END DO
				DO icnt=(itn-1),1,-1
					TriU(icnt)=TriU(icnt)-gam(icnt+1)*TriU(icnt+1)
				END DO
			END IF
			tridag=.TRUE.
			RETURN
	
		END FUNCTION tridag
	
	END SUBROUTINE caltridag
	
	SUBROUTINE ReadCSV(fname,rVar,ncnt,iNAValue)
		!Generic subroutine to read csv file
		USE DBEM,ONLY	:	cfolderpath
		IMPLICIT NONE
		CHARACTER(40), INTENT(IN)	:: fname
		INTEGER, INTENT(IN)	::	ncnt,iNAValue
		REAL, INTENT(OUT)	:: rVar(ncnt)
		INTEGER	::	i
	
		OPEN (UNIT=1,file=TRIM(cfolderpath)//fname)
		DO i=1,259200
			READ(1,*)	rVar(i)
			IF(rVar(i)==(-9999)) THEN
				rVar(i)=iNAValue
			END IF
		END DO
		CLOSE(1)
		RETURN
	END SUBROUTINE ReadCSV
	
	SUBROUTINE ReadFMort(fname,rVar,ncnt,iNAValue)
		!Generic subroutine to read csv file
		USE DBEM,ONLY	:	cfolderpath
		IMPLICIT NONE
		CHARACTER(40), INTENT(IN)	:: fname
		INTEGER, INTENT(IN)	::	ncnt,iNAValue
		REAL, INTENT(OUT)	:: rVar(ncnt)
		INTEGER	::	i
	
		OPEN (UNIT=1,file=TRIM(cfolderpath)//fname)
		DO i=1,259200
			READ(1,'(F6.3)')	rVar(i)
			IF(rVar(i)==(-9999)) THEN
				rVar(i)=iNAValue
			END IF
		END DO
		CLOSE(1)
		RETURN
	END SUBROUTINE ReadFMort
	
	SUBROUTINE ReadCSVDouble(fname,dVar,ncnt,iNAValue)
		!Generic subroutine to read csv file
		USE DBEM,ONLY	:	cfolderpath
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN)	:: fname
		INTEGER, INTENT(IN)	::	ncnt,iNAValue
		DOUBLE PRECISION, INTENT(OUT)	:: DVar(ncnt)
		INTEGER	::	i
		OPEN (UNIT=1,file=TRIM(cfolderpath)//fname)
		DO i=1,259200
			READ(1,*)	DVar(i)
			IF(DVar(i)==(-9999)) THEN
				DVar(i)=iNAValue
			END IF
		END DO
		CLOSE(1)
		RETURN
	END SUBROUTINE ReadCSVDouble
	
	SUBROUTINE ReadLatLon(fname)
		!Subroutine to read latitude and longitude
		USE DBEM, ONLY : lon, lat,cfolderpath
		IMPLICIT NONE
		REAL	:: rtemp
		CHARACTER(*), INTENT(IN)	:: fname
		INTEGER	::	i,n=SIZE(lon)
		OPEN (UNIT=1,file=TRIM(cfolderpath)//fname)
		DO i=1,n
			READ(1,*)	rtemp,lon(i),lat(i)
		END DO
		CLOSE(1)
		RETURN
	END SUBROUTINE ReadLatLon
	
	SUBROUTINE WriteCSV(fname,rVar,itn)
		!Subroutine to write to CSV
		USE DBEM,ONLY	:	cfolderpath,cResultspath,iSpp,iTaxonList
		IMPLICIT NONE
		CHARACTER(*),INTENT(IN)		:: fname
		CHARACTER(6)	::	cTemp
		REAL,INTENT(IN)	:: rVar(*)
		INTEGER,INTENT(IN)	:: itn
		INTEGER	::	i
	
		Write(cTemp, '(i6)') iTaxonList(iSpp)
		OPEN (UNIT=1,file='/scratch/jepa/'//'Results/'//TRIM(cResultspath)//'/'//TRIM(cTemp)//'/'//fname) !JEPA: Change user name according to computecanada
		!for Cedar
		!OPEN (UNIT=1,file='/home/wailung/projects/rrg-wailung/CMIP6/DBEM_outputs/'//'txt/'//TRIM(cResultspath)//'/'//TRIM(cTemp)//'/'//fname) !JEPA: change user name according to computecanada
	
		!OPEN (UNIT=1,file='/Users/wwlcheung/Fortran/Data/'//'Results/'//TRIM(cResultspath)//'/'//TRIM(cTemp)//'/'//fname)
		DO i=1,itn
			IF(rVar(i)>0 .OR. ISNAN(rVar(i))) THEN
				WRITE (1,*)	i,rVar(i)
			END IF
		END DO
		CLOSE(1)
	END SUBROUTINE WriteCSV
	
	SUBROUTINE WriteCSVDouble(fname,rVar,itn)
		!Subroutine to write to CSV
		USE DBEM,ONLY	:	cfolderpath,iSpp,iTaxonList,cResultspath
		IMPLICIT NONE
		CHARACTER(*),INTENT(IN)		:: fname
		CHARACTER(6)	::	cTemp
		DOUBLE PRECISION,INTENT(IN)	:: rVar(*)
		INTEGER,INTENT(IN)	:: itn
		INTEGER	::	i
	
		Write(cTemp, '(i6)') iTaxonList(iSpp)
		OPEN (UNIT=1,file='/scratch/jepa/'//'Results/'//TRIM(cResultspath)//'/'//TRIM(cTemp)//'/'//fname) !JEPA: change user name according to computecanada
		!for cedar
		!OPEN (UNIT=1,file='/home/wailung/projects/rrg-wailung/CMIP6/DBEM_outputs/'//'txt/'//TRIM(cResultspath)//'/'//TRIM(cTemp)//'/'//fname) !JEPA change user name according to computecanada
	
		!OPEN (UNIT=1,file='/Users/wwlcheung/Fortran/Data/'//'Results/'//TRIM(cResultspath)//'/'//TRIM(cTemp)//'/'//fname)
		DO i=1,itn
			IF(rVar(i)>0 .OR. ISNAN(rVar(i))) THEN
				WRITE (1,*)	i,rVar(i)
			END IF
		END DO
		CLOSE(1)
	END SUBROUTINE WriteCSVDouble
	
	!!-------------- TRAVIS13 -----------------!!
	!!---------------------------------- Include rTHQ10Mort and rTHQ10LarvMort in Subroutine line (T=Temporary) -------------!!
	SUBROUTINE ReadTaxonPara(fname,rTLinf,rTVBonK,rTLwA,rTLwB,rTArrhCoef,rTAlloScal,rTHQ10Met,rTHQ10Mort,rTHQ10LarvMort)
		USE DBEM, ONLY	:	iSpp,cfolderpath,cTaxonDatapath,rtaxMinD,rtaxMaxD,rDiffCoef,lerrTaxon,rIntR,cDemPel,rHabAssocIce
		USE DBEM, ONLY	:	rHabAssocSal,rMaxCatch,rCoral,rUpwelling,rInshore,rOffshore,rShelf
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN)	:: fname
	
		!!-------------- Set up temporary variables for pH adult and larval mortality  -----------------!!
		REAL,INTENT(OUT)	::	rTLinf,rTVBonK,rTLwA,rTLwB,rTArrhCoef(2),rTAlloScal(2),rTHQ10Met,rTHQ10Mort,rTHQ10LarvMort
		!!-------------- TRAVIS13 -----------------!!
	
		CHARACTER(20)	:: ctemp1,ctemp2
		REAL		::	rtemp,rtemp2(5)
		INTEGER	::	iIO,icnt,j
		lerrTaxon=.FALSE.
		rtemp2=0
		OPEN (UNIT=1,file=TRIM(cfolderpath)//'/Species/'//TRIM(cTaxonDatapath)//'/'//fname)
			DO icnt=1,26
				SELECT CASE (icnt)
					CASE (1)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (2)	!Taxon name
						READ(1,*,IOSTAT=iIO)	ctemp1,ctemp2
						IF (iIO==0) THEN
							WRITE(*,*)	'Running Taxon '//ctemp2
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (3)	!Demersal/Pelagic
						READ(1,*,IOSTAT=iIO)	ctemp1,ctemp2
						IF (iIO==0) THEN
							cDemPel(iSpp)=ctemp2
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (4)	!Maximum depth
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rtaxMaxD(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (5)	!Minimum depth
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rtaxMinD(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (6)	!Diffcoef
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rDiffCoef(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (7)	!Intrinsic R
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rIntR(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (8)	!Linf
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTLinf=rtemp
						ElSE
							lerrTaxon=.TRUE.
							WRITE(*,*)	'Winf error'
							EXIT
						END IF
					CASE (9)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTVBonK=rtemp
						ElSE
							lerrTaxon=.TRUE.
							WRITE(*,*)	'rTVBonK error'
							EXIT
						END IF
					CASE (10)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTLwA=rtemp
						ElSE
							lerrTaxon=.TRUE.
							WRITE(*,*)	'rTLwA error'
							EXIT
						END IF
					CASE (11)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTLwB=rtemp
						ElSE
							lerrTaxon=.TRUE.
							WRITE(*,*)	'rTLwB error'
							EXIT
						END IF
					CASE (12)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTArrhCoef(1)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (13)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTArrhCoef(2)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (14)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTAlloScal(1)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (15)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTAlloScal(2)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (16)	!Ice habitat
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rHabAssocIce(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (17)	!Salinity habitat
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp2(1),rtemp2(2),rtemp2(3),rtemp2(4),rtemp2(5)
						IF (iIO==0) THEN
							!WRITE(*,*)	rtemp2(1)
							rHabAssocSal(iSpp,1)=rtemp2(1)
							rHabAssocSal(iSpp,2)=rtemp2(2)
							rHabAssocSal(iSpp,3)=rtemp2(3)
							rHabAssocSal(iSpp,4)=rtemp2(4)
							rHabAssocSal(iSpp,5)=rtemp2(5)
							!WRITE(*,*) 'Habitat',	rHabAssocSal(1,1)
						ElSE
							lerrTaxon=.TRUE.
							WRITE(*,*)	'Habitat error'
							EXIT
						END IF
	
					CASE (18)	!Coral
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rCoral(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (19) !Upwelling
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rUpwelling(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (20) !Inshore
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rInshore(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (21) !Shelf
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rShelf(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (22) !Offshore
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rOffshore(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (23)	!pH sensitivity (metabolism)
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTHQ10Met=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (24)	!Max catch
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rMaxCatch(iSpp)=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					!!----------- TRAVIS14 --------------!!
					!!-------------------------ADD CASE (25) pH adult mortality and CASE (26) pH larval mortality------------------------!!
					CASE (25)	!pH adult mortality
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTHQ10Mort=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (26)	!pH larval mortality
						READ(1,*,IOSTAT=iIO)	ctemp1,rtemp
						IF (iIO==0) THEN
							rTHQ10LarvMort=rtemp
						ElSE
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					!!----------- TRAVIS14 --------------!!
	
				END SELECT
			END DO
		CLOSE(1)
		RETURN
	END SUBROUTINE
	
	SUBROUTINE ReadSettings(fname)
		!Subroutine to read in list of species to run
		USE DBEM, ONLY	:	SppNo,cfolderpath,cCCScenario,cSppListFile,cResultspath,cSetfilename,cTaxonDatapath,cSSP,cHSmapfile
		USE DBEM, ONLY	: 	rFfactorHS,rFfactorEEZ
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN)	:: fname
		INTEGER*1	:: itemp,ifilenumber
		CHARACTER(20)	::	ctemp,ctemp2
		LOGICAL	::	lerrTaxon
		INTEGER	::	iIO,icnt,ifile
	
		lerrTaxon=.FALSE.
		icnt=0
		!OPEN (UNIT=1,file=TRIM(cfolderpath)//TRIM(cSetfilename)//'.txt',ACTION="READ")
		OPEN (UNIT=1,file='../../scripts/run_dbem/settings.txt',ACTION="READ")
			DO icnt=1,10
				SELECT CASE (icnt)
					CASE (1) !Number of Species
	
						READ(1,*,IOSTAT=iIO)	ctemp,itemp
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
						SppNo=itemp
					CASE (2) !Climate change scenario folder name
	
						READ(1,*,IOSTAT=iIO)	ctemp,cCCScenario
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (3) !Shared Socio-economic Pathways
	
						READ(1,*,IOSTAT=iIO)	ctemp,cSSP
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (4) !Species list
	
						READ(1,*,IOSTAT=iIO)	ctemp,cSppListFile
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (5) !Result directory
	
						READ(1,*,IOSTAT=iIO)	ctemp,cResultspath
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (6) !Taxon metadata directory
	
						READ(1,*,IOSTAT=iIO)	ctemp,cTaxonDatapath
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (7) !filenumber
	
						READ(1,*,IOSTAT=iIO)	ctemp,ifilenumber
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (8)	!fishing mortality in HS
						READ(1,*,IOSTAT=iIO)	ctemp,rFfactorHS
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
	
					CASE (9)	!fishing mortality in EEZ
						READ(1,*,IOSTAT=iIO)	ctemp,rFfactorEEZ
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
					CASE (10)	!Filepath for HS MPA Scenario
						READ(1,*,IOSTAT=iIO)	ctemp,cHSmapfile
						IF (iIO/=0) THEN
							lerrTaxon=.TRUE.
							EXIT
						END IF
				END SELECT
			END DO
		CLOSE(1)
		WRITE(*,*) 'Species list: ',ifilenumber
		WRITE(ctemp2,'(i2)')	ifilenumber
		ctemp=cSppListFile
		cSppListFile=TRIM(cSppListFile)//TRIM(ctemp2)//'.txt'
		WRITE(*,*)	cSppListFile
		WRITE(ctemp2,'(i2)')	ifilenumber+1
	
		!OPEN (UNIT=1,file=TRIM(cfolderpath)//TRIM(cSetfilename)//'.txt', &
		OPEN (UNIT=1,file='../../scripts/run_dbem/settings.txt', &
		!OPEN (UNIT=1,file='../../scripts/run_dbem/settings.txt',ACTION="READ")
		ACTION="readwrite")
			WRITE(1,*) 'SppNo ',SppNo
			WRITE(1,*) 'CCSc ',cCCScenario
			WRITE(1,*) 'SSP ', cSSP
			WRITE(1,*) 'rsfile ', ctemp
			WRITE(1,*) 'rpath ', cResultspath
			WRITE(1,*) 'tpath ', cTaxonDatapath
			WRITE(1,*)	'ifile ',ctemp2
			WRITE(1,'(A5,F7.2)')	' FHS '	,rFfactorHS
			WRITE(1,'(A6,F7.2)')	' FEEZ '	,rFfactorEEZ
			WRITE(1,*) 'MPApath ',	cHSmapfile
		CLOSE(1)
	END SUBROUTINE
	
	SUBROUTINE ReadTaxonList(fname)
		!Subroutine to read in list of species to run
		USE DBEM, ONLY	:	SppNo,cfolderpath,iTaxonList,cSppListFile
		IMPLICIT NONE
		CHARACTER(*), INTENT(IN)	:: fname
		INTEGER	:: itemp
		!CHARACTER(6)	::	ctemp
		LOGICAL	::	lerrTaxon
		INTEGER	::	iIO,icnt
	
		lerrTaxon=.FALSE.
		icnt=0
		OPEN (UNIT=1,file=TRIM(cfolderpath)//'Species/TaxonList/'//TRIM(cSppListFile),ACTION="READ")
		!OPEN (UNIT=1,file='/Data/RunSppList.csv',ACTION="READ")
			DO icnt=1,SppNo
				READ(1,'(i6)',IOSTAT=iIO)	itemp
				IF (iIO/=0) THEN
					lerrTaxon=.TRUE.
					EXIT
				END IF
				iTaxonList(icnt)=itemp
			END DO
		CLOSE(1)
		!WRITE(*,*) itemp
	
	END SUBROUTINE
	
	SUBROUTINE SetTempPreference
		!Subroutine to set the temperature preference (unit = relative abundance/area vs temperature)
		USE DBEM,ONLY	: iIniStartYear,iIniEndYear,rAbd,pwater,iTempBinNum,rTempBinSize,rTempBin, &
		ncol,nrow,rtAbdBin, iSpp,rtLTemp,rtUTemp,rIniSST,rIniBtT,cDemPel
		IMPLICIT NONE
		INTEGER	::	iYr,icnt,icnt2,i
		REAL,DIMENSION(ncol*nrow)	:: rST, rCellTemp,rDummyAbd,rtArea
		CHARACTER(4)	:: ctemp
		REAL	::	rTotalAbd
		REAL,DIMENSION(iTempBinNum)	::	rtAreaBin
	
		rtArea=0
		rtAreaBin=0
		!loop through all years of historical data
		DO icnt=1, iTempBinNum
			IF (icnt==1) THEN
				rTempBin(icnt)=(-5-rTempBinSize/2)
			ELSE
				rTempBin(icnt)=rTempBin(icnt-1)+rTempBinSize
			END IF
		END DO
		rtAbdBin=0
		rtUTemp=-9999
		rtLTemp=9999
		DO	iYr=iIniStartYear,iIniEndYear
			rST=0
			icnt2=0
			rTotalAbd=0
			rDummyAbd=0
			IF(cDemPel(iSpp)=='P') THEN
				rST=rIniSST(:,(iYr-iIniStartYear+1))
			ELSE
				rST=rIniBtT(:,iYr-(iIniStartYear+1))
			END IF
	
			DO icnt=1,ncol*nrow
				IF (rAbd(icnt)>0 .AND. rST(icnt)/=-9999) THEN
					icnt2=icnt2+1
					rCellTemp(icnt2)=rST(icnt)
					rDummyAbd(icnt2)=rAbd(icnt)
					rtUTemp=MAX(rtUTemp,rST(icnt))
					rtLTemp=MIN(rtLTemp,rST(icnt))
					rtArea(icnt2)=pwater(icnt)
				END IF
			END DO
			DO icnt=1,icnt2
				DO i=1, iTempBinNum
					IF (rCellTemp(icnt)>=(rTempBin(i)-rTempBinSize/2) .AND. (rCellTemp(icnt)<(rTempBin(i)+rTempBinSize/2))) THEN
						rtAbdBin(i)=rtAbdBin(i)+rDummyAbd(icnt)
						rtAreaBin(i)=rtAreaBin(i)+rtArea(icnt)
					END IF
				END DO
			END DO
		END DO
	
		!Smooth temperature preference
		CALL FindSmoothPara
	
		rTotalAbd=0
		DO icnt=1,iTempBinNum
			If(rtAreaBin(icnt)>0) THEN
				rTotalAbd=rTotalAbd+rtAbdBin(icnt)/rtAreaBin(icnt)
				rtAbdBin(icnt)=rtAbdBin(icnt)/rtAreaBin(icnt)
			END IF
		END DO
		!Use absolute per area relative abundance
		rtAbdBin=rtAbdBin !/rTotalAbd
	
		!DO	icnt=1,iTempBinNum
		!	WRITE(*,*)	rTempBin(icnt),rtAbdBin(icnt)
		!END DO
	
		!WRITE(*,*)	rtLTemp,rtUTemp
	END SUBROUTINE SetTempPreference
	
	SUBROUTINE FindSmoothPara
		USE DBEM,ONLY	:	rtAbdBin,iTempBinNum,rTempBin,iSpp
		IMPLICIT NONE
		INTEGER	::	icnt,iPeakSum,icntKernel,ismoothcnt
		REAL ::	rtempsum1,rtempsum2,rtAbd(iTempBinNum),rAbdDiff,rtAbd2(iTempBinNum)
	
		iPeakSum=0
		ismoothcnt=1
		rtAbd=0
		rtAbd2=0
		rtAbd2=rtAbdBin
		DO WHILE (iPeakSum==0)
			!Smoothing function
			DO icnt=1,iTempBinNum
				rtempsum1=0
				rtempsum2=0
				IF (ismoothcnt==0) THEN
					rtAbd(icnt)=rtAbdBin(icnt)
				ELSE
					DO	icntKernel=INT(icnt-ismoothcnt/2),INT(icnt+ismoothcnt/2)
						IF (icntKernel<1 .OR. icntKernel>iTempBinNum) THEN
							rtempsum2=rtempsum2+1
						ELSE
							rtempsum1=rtempsum1+rtAbdBin(icntKernel)
							rtempsum2=rtempsum2+1
						END IF
					END DO
					rtAbd(icnt)=rtempsum1/rtempsum2
				END IF
			END DO
			DO icnt=1,(iTempBinNum-2)
				rAbdDiff = (rtAbd(icnt + 1) - rtAbd(icnt - 1)) * (rtAbd(icnt + 2) - rtAbd(icnt))
				IF (rAbdDiff < 0) Then
					iPeakSum = iPeakSum + 1
				End If
			END DO
			ismoothcnt=ismoothcnt+2
		END DO
		rtAbdBin=rtAbd
	
	END SUBROUTINE
	
	SUBROUTINE CalK(lInitiation,rtST,rtLinf)
		!Subroutine to calculate habitat suitability in each cell
		USE	DBEM,ONLY	:	rAbdBin,rTempBin,rLTemp,rUTemp,rTempBinSize,nrow,ncol
		USE	DBEM,ONLY	:	iSpp,rHabSuit,pwater,rK,rAbd,rKMax,rtIce,rHabAssocIce,rHabSuitMax
		USE DBEM,ONLY	:	rHabAssocSal,rSalPScale,cDemPel,rCoral,rCoralMap,rUpWellMap,rUpwelling
		USE DBEM,ONLY	:	rtaxMinD,rtaxMaxD,rShelf,rInshore,rOffshore
		IMPLICIT NONE
		LOGICAL,INTENT(IN)	:: lInitiation
		REAL,INTENT(IN)	::	rtST(nrow*ncol)
		REAL,INTENT(IN)	::	rtLinf
		INTEGER	::	icnt,j
		REAL	::	rTempProb,rtMax,rDepthProb,rIceProb,rSalProb
		REAL	::	rtempHabSuit(nrow*ncol)
		REAL	:: rTempCoral,rTempUpwell
	
		rtMax=0
		!Calculate habitat suitability base on temperature (rAbdBin)
		DO icnt=1,nrow*ncol
			IF (pwater(icnt)>0 .AND. pwater(icnt)/=-9999) THEN
				CALL EstEnvtSuit(rtST(icnt),rAbdBin,rTempBin,rLTemp,rUTemp,rTempBinSize,rTempProb)
			Else
				rTempProb=0
			END IF
			rtempHabSuit(icnt)=rTempProb
		END DO
		!WRITE(*,*) rLTemp,rUTemp
		!DO icnt=1,SIZE(rAbdBin(iSpp,:))
		!	WRITE(*,*) rTempBin(icnt),rAbdBin(iSpp,icnt)
		!END DO
		!*****************************
		!Calculate carrying capacity
		!******************************
		!Calculate maximum temperature preference value
		!rtMax=MAXVAL(rAbdBin)
	
		!if first time step, initialize variables
		IF (lInitiation .EQV. .TRUE.) THEN
			!Initialize variable
			rHabSuit(:,iSpp)=0
			rK(:,iSpp)=0
	
			DO icnt=1,nrow*ncol
				IF(pwater(icnt)>0) THEN
					IF(cDemPel(iSpp)=='D') THEN
						rDepthProb=getDepthValue(icnt,rtaxMinD(iSpp),rtaxMaxD(iSpp))
					ELSE
						rDepthProb=1
						IF(rOffshore(iSpp)<=0.25) THEN
							IF(rInshore(iSpp)>=0.25 .OR. rShelf(iSpp)>=0.25) THEN
								rDepthProb=getDepthValue(icnt,1.0,300.0)
							ELSE
								rDepthProb=1
							END IF
						END IF
					END IF
	
					rIceProb=getIceValue(icnt)
					rSalProb=0
					DO j=1,5
						rSalProb=rHabAssocSal(iSpp,j)*rSalPScale(icnt,j)+rSalProb
					END DO
					rSalProb=MAX(0.0,rSalProb)
					!WRITE(*,*)	(rHabAssocSal(iSpp,j),j=1,5)
					!WRITE(*,*)	(rSalPScale(icnt,j),j=1,5)
					!WRITE(*,*)	rSalProb
	
					!Coral habitat suirtability
					IF(rCoral(iSpp)/=0.0) THEN
						rTempCoral=getCoralValue(icnt)
						rTempCoral=MIN(1.0,rTempCoral)
						rTempCoral=1+rTempCoral*rCoral(iSpp)
					ELSE
						rTempCoral=1
					END IF
	
	
					!Upwelling
					!IF(rUpwelling(iSpp)>0.5) THEN
					!	rTempUpwell=rUpWellMap(icnt)
					!ELSE
						rTempUpwell=1
					!END IF
	
					IF(rK(icnt,iSpp)==0.0 .AND. pwater(icnt)>0.0) THEN
						rK(icnt,iSpp)=rtempHabSuit(icnt)*pwater(icnt)*rDepthProb*rIceProb*rSalProb*rTempUpwell
					END IF
					!Store habitat suitability for next time step
					rHabSuit(icnt,iSpp)=rtempHabSuit(icnt)*rDepthProb*rIceProb*rSalProb*rTempCoral
				END IF
			END DO
			rKMax(iSpp)=MAXVAL(rK(:,iSpp))
			rHabSuitMax(iSpp)=MAXVAL(rHabSuit(:,iSpp))
			DO icnt=1,nrow*ncol
				IF(rK(icnt,iSpp)<(rKMax(iSpp)*0.01)) THEN
					rK(icnt,iSpp)=0
				END IF
			END DO
				rHabSuit(:,iSpp)=rHabSuit(:,iSpp)/rHabSuitMax(iSpp)
		ELSE
			!if K>0, then calculate new K based on change in habitat suitability
			!Otherwise, base on average abundance in habitat suitability bin
			DO icnt=1,nrow*ncol
				IF(pwater(icnt)>0) THEN
	
					!Depth
					IF(cDemPel(iSpp)=='D') THEN
						rDepthProb=getDepthValue(icnt,rtaxMinD(iSpp),rtaxMaxD(iSpp))
					ELSE
						rDepthProb=1
						IF(rOffshore(iSpp)<=0.25) THEN
							IF(rInshore(iSpp)>=0.25 .OR. rShelf(iSpp)>=0.25) THEN
								rDepthProb=getDepthValue(icnt,1.0,300.0)
							ELSE
								rDepthProb=1
							END IF
						END IF
					END IF
	
					!Salinity
					rSalProb=0
					DO j=1,5
						rSalProb=rHabAssocSal(iSpp,j)*rSalPScale(icnt,j)+rSalProb
					END DO
					rSalProb=MAX(0.0,rSalProb)
					rIceProb=getIceValue(icnt)
					IF(rCoral(iSpp)/=0.0) THEN
						rTempCoral=getCoralValue(icnt)
						rTempCoral=MIN(1.0,rTempCoral)
						rTempCoral=1+rTempCoral*rCoral(iSpp)
					ELSE
						rTempCoral=1
					END IF
	
					!Upwelling
					!IF(rUpwelling(iSpp)>0.5) THEN
					!	rTempUpwell=rUpWellMap(icnt)
					!ELSE
						rTempUpwell=1
					!END IF
	
					rK(icnt,iSpp)=rtempHabSuit(icnt)*pwater(icnt)*rDepthProb*rIceProb*rSalProb
					IF(rK(icnt,iSpp)<(rKMax(iSpp)*0.01)) THEN
						rK(icnt,iSpp)=0
					END IF
					rHabSuit(icnt,iSpp)=rtempHabSuit(icnt)*rDepthProb*rIceProb*rSalProb*rTempCoral*rTempUpWell
					rHabSuit(icnt,iSpp)=rHabSuit(icnt,iSpp)/rHabSuitMax(iSpp)
	
					!testing
					!IF(icnt==36999) THEN
					!	WRITE(*,*) rtempHabSuit(icnt),rDepthProb,rIceProb,rSalProb,rTempCoral,rTempUpWell
					!END IF
				END IF
			END DO
		END IF
	
		RETURN
	
		CONTAINS
	
		REAL FUNCTION getDepthValue(itempseq,rtMinD,rtMaxD)
			!Function to calculate depth suitability
			USE DBEM,ONLY	: rMaxElev,rMinElev,rAvgElev
			REAL	::	rDepthVal,rMaxD,rMinD,rEleMax,rEleMin,rEleAvg
			INTEGER,INTENT(IN)	::	itempseq
			REAL,INTENT(IN)	::	rtMinD,rtMaxD
			rMaxD = -1*rtMaxD
			rMinD = -1*rtMinD
			rEleMax = rMaxElev(itempseq)
			rEleMin = rMinElev(itempseq)
			rEleAvg = rAvgElev(itempseq)
	
			IF ((rMaxD <= rEleMax .AND. rMaxD >= rEleMin) .OR. (rMinD <= rEleMax .AND. rMaxD >= rEleMin) &
			.OR. (rMinD >= rEleMin .AND. rMaxD <= rEleMax) .OR. (rMinD >= rEleMin .AND. rMinD <= rEleMax)) THEN
				rDepthVal = 1
			ELSE
				rDepthVal = 0
			End IF
	
			getDepthValue=rDepthVal
	
		END FUNCTION getDepthValue
	
		REAL FUNCTION getIceValue(itempseq)
			!Function to calculate ice-related suitability
			USE DBEM,ONLY	:rtIce,rHabAssocIce
			INTEGER,INTENT(IN)	:: itempseq
			REAL	::	rIceVal
	
			!If ice association is = 0, then ice habitat is negatively and directly proportional
			!to ice extent(%)
			!IF(rtIce(itempseq)>0.5) THEN
			!	rtIce(itempseq)=1
			!END IF
			!testing
			!IF(itempseq==36999) THEN
			!	WRITE(*,*) rHabAssocIce(iSpp),rtIce(36999)
			!END IF
			IF(rHabAssocIce(iSpp)==0) THEN
				rIceVal=1-rtIce(itempseq)
			ELSE
				rIceVal=MAX(rtIce(itempseq)*rHabAssocIce(iSpp)+(1-rHabAssocIce(iSpp)),0.1)
			END IF
			getIceValue=rIceVal
		END FUNCTION getIceValue
		REAL FUNCTION getCoralValue(itempseq)
			!Function to calculate coral-related suitability
			!Assume a non-linear relationship with = 7*Coral/(1+6*Coral)
			USE DBEM, ONLY : rCoralMap,pwater
			INTEGER,INTENT(IN)	::	itempseq
			REAL	::	rcoraltemp
	
			rcoraltemp=rCoralMap(itempseq)/pwater(itempseq)
	
			getCoralValue=rcoraltemp*7/(1+6*rcoraltemp)
	
		END FUNCTION getCoralValue
	END SUBROUTINE CalK
	
	SUBROUTINE EstEnvtSuit(rVar,rPDF,rPDFGroup,rLTT,rUTT,rtBinSize,rProb)
		!Function to calculate the probability of occupancy given either temperature or depth
		!Re-coded by William on 27 Dec 2013
		USE DBEM, ONLY	:	iSpp,SppNo,iTempBinNum
		IMPLICIT NONE
		REAL,INTENT(IN)	::	rVar,rPDF(SppNo,iTempBinNum),rPDFGroup(iTempBinNum),rLTT(SppNo),rUTT(SppNo),rtBinSize
		REAL,INTENT(OUT)	::	rProb
		INTEGER	::	icnt
	
		DO icnt=1,iTempBinNum
			IF (icnt==1) THEN
				IF (rVar<rPDFGroup(icnt)) THEN
					IF (rVar>=rLTT(iSpp) .AND. rVar<rPDFGroup(icnt)) THEN
						rProb=(rVar-(rPDFGroup(icnt)-rtBinSize/2))*rPDF(iSpp,icnt)
						rProb=rProb/(rPDFGroup(icnt)-(rPDFGroup(icnt)-rtBinSize/2))
						rProb=MAX(0.0,rProb)
					ELSEIF	(rVar<rLTT(iSpp)) THEN
						rProb=0
					END IF
				ELSE IF (rVar>=rPDFGroup(icnt) .AND. rVar<rPDFGroup(icnt)+rtBinSize/2) THEN
					rProb = rPDF(iSpp,icnt)
					rProb=rProb+(rVar-rPDFGroup(icnt))*(rPDF(iSpp,icnt+1)-rPDF(iSpp,icnt))/(rPDFGroup(icnt + 1)-rPDFGroup(icnt))
					rProb=MAX(0.0,rProb)
				END IF
			ELSE
				IF (rVar<rPDFGroup(icnt) .AND. rVar>=rPDFGroup(icnt)-rtBinSize/2) THEN
					rProb = rPDF(iSpp,icnt - 1)
					rProb=rProb+(rVar-rPDFGroup(icnt - 1))*(rPDF(iSpp,icnt)-rPDF(iSpp,icnt-1))/(rPDFGroup(icnt)-rPDFGroup(icnt-1))
					rProb=MAX(0.0,rProb)
				ELSE IF (rVar>=rPDFGroup(icnt) .AND. rVar<rPDFGroup(icnt)+rtBinSize/2) THEN
					IF (icnt/=iTempBinNum) THEN
						rProb=rPDF(iSpp,icnt)
						rProb=rProb+(rVar-rPDFGroup(icnt))*(rPDF(iSpp,icnt+1)-rPDF(iSpp,icnt))/(rPDFGroup(icnt+1)-rPDFGroup(icnt))
						rProb=MAX(0.0,rProb)
					ELSE
						rProb=(rVar - rPDFGroup(icnt))*rPDF(iSpp,icnt)/(rPDFGroup(icnt)+rtBinSize/2-rPDFGroup(icnt))
						rProb = rPDF(iSpp,icnt)+rProb
						rProb=MAX(0.0,rProb)
					END	IF
				END IF
			END IF
		END DO
		IF (rVar>rUTT(iSpp)) THEN
			rProb=0
		END IF
		IF (rVar<rLTT(iSpp)) THEN
			rProb=0
		END IF
		RETURN
	END SUBROUTINE EstEnvtSuit
	
	SUBROUTINE	 CalDeltaAbd(dlTempCatch,rMTemp,rGCoef,dTCatch)
		!Subroutine to calculate change in abundance based on logistic growth model and rk4
		!numerical integration
		USE DBEM,ONLY	:	rK,rIntR,rAbd,inTimeStep,iSpp,nrow,ncol
	
		IMPLICIT NONE
		REAL	::	rDeltaAbd
		REAL,PARAMETER	::	rhpara=0.1
		DOUBLE PRECISION	::	rXTemp(0:CEILING(1/rhpara)),rtk1,rtk2,rtk3,rtk4
		INTEGER	::	icnt,itempseq,icntrk
		REAL, INTENT(IN)	::	rGCoef(259200),rMTemp(259200)	!0 = no recruitment
		REAL, INTENT(IN)	::	dlTempCatch(nrow*ncol)
		DOUBLE PRECISION, INTENT(OUT)	::	dTCatch(259200)
		DOUBLE PRECISION	:: rTCR,rCTemp(0:CEILING(1/rhpara)),dTCatch1,dTCatch2,dTCatch3,dTCatch4
	
	
		!WRITE(*,*) dlTempCatch
		icntrk=CEILING(1/rhpara)
		DO itempseq=1,nrow*ncol
			!testing
			!IF(dlTempCatch(itempseq)>0) THEN
			!	WRITE(*,*) 'FishMort', dlTempCatch(itempseq)
			!END IF
			IF (rK(itempseq,iSpp)>0) THEN
				rXTemp(0)=rAbd(itempseq)
				rCTemp(0)=0.0
				DO icnt=1,icntrk
					rTCR=rXTemp(icnt-1)/rK(itempseq,iSpp)
	
					!IF(rTCR>2.0) THEN
					!	rTCR=2.0
					!END IF
	
	
					rtk1=rXTemp(icnt-1)*rIntR(iSpp)*rGCoef(itempseq)*(1/inTimeStep)*(1-rTCR)
					rtk1=rtk1-(1-EXP(-dlTempCatch(itempseq)*(1/inTimeStep)))*rxTemp(icnt-1) !(itempseq)
					rtk1=rhpara*rtk1
					dtCatch1=(1-EXP(-dlTempCatch(itempseq)*(1/inTimeStep)))*rxTemp(icnt-1)
					dTCatch1=rhpara*dTCatch1
					!IF(isnan(rtk1*rtk1)) THEN
						!rtk1=rtk1-rtest*rXTemp(icnt-1)
					!	WRITE(*,*) 'error',icnt,rAbd(itempseq),rK(itempseq,iSpp),rXTemp(icnt-1),rtest*rXTemp(icnt-1)
					!END IF
	
	
					rTCR=(rXTemp(icnt-1)+rtk1/2)/rK(itempseq,iSpp)
	
					!IF(rTCR>2.0) THEN
					!	rTCR=2.0
					!END IF
					rtk2=(rXTemp(icnt-1)+rtk1/2)*rIntR(iSpp)*rGCoef(itempseq)*(1/inTimeStep+rhpara/2)*(1-rTCR)
					rtk2=rtk2-(1-EXP(-dlTempCatch(itempseq)*(1/inTimeStep+rhpara/2)))*(rXTemp(icnt-1)+rtk1/2)	!(itempseq)
					rtk2=rhpara*rtk2
	
					dTCatch2=(1-EXP(-dlTempCatch(itempseq)*(1/inTimeStep+rhpara/2)))*(rXTemp(icnt-1)+rtk1/2)
					dTCatch2=rhpara*dTCatch2
					rTCR=(rXTemp(icnt-1)+rtk2/2)/rK(itempseq,iSpp)
	
					!IF(rTCR>2.0) THEN
					!	rTCR=2.0
					!END IF
	
					rtk3=(rXTemp(icnt-1)+rtk2/2)*rIntR(iSpp)*rGCoef(itempseq)*(1/inTimeStep+rhpara/2)*(1-rTCR)
					rtk3=rtk3-(1-EXP(-dlTempCatch(itempseq)*(1/inTimeStep+rhpara/2)))*(rXTemp(icnt-1)+rtk2/2)	!(itempseq)
					rtk3=rhpara*rtk3
	
					dTCatch3=(1-EXP(-dlTempCatch(itempseq)*(1/inTimeStep+rhpara/2)))*(rXTemp(icnt-1)+rtk2/2)
					dTCatch3=rhpara*dTCatch3
					rTCR=(rXTemp(icnt-1)+rtk3)/rK(itempseq,iSpp)
	
					!IF(rTCR>2.0) THEN
					!	rTCR=2.0
					!END IF
	
					rtk4=(rXTemp(icnt-1)+rtk3)*rIntR(iSpp)*rGCoef(itempseq)*(1/inTimeStep+rhpara)*(1-rTCR)
					rtk4=rtk4-(1-EXP(-dlTempCatch(itempseq)*(1/inTimeStep+rhpara)))*(rXTemp(icnt-1)+rtk3)	!(itempseq)
					rtk4=rhpara*rtk4
	
					dTCatch4=(1-EXP(-dlTempCatch(itempseq)*(1/inTimeStep+rhpara)))*(rXTemp(icnt-1)+rtk3)
					dTCatch4=rhpara*dTCatch4
	
					rXTemp(icnt)=rXTemp(icnt-1)+(rtk1+2*rtk2+2*rtk3+rtk4)/6
					rCTemp(icnt)=rCTemp(icnt-1)+(dTCatch1+2*dTCatch2+2*dTCatch3+dTCatch4)/6
					!if (rtest*rXTemp(icnt)/=0 .OR. rXTemp(icnt)>100000000.00 .OR. rXTemp(icnt)<(-100000000.00)) THEN
					!	WRITE(*,*) 'Big no',icnt,rxTemp(icnt),rK(itempseq,iSpp)
					!	WRITE(*,*) 'rk para',rhpara,rXTemp(icnt-1),rIntR(iSpp),inTimeStep
					!	WRITE(*,*) 'rk para2',rtk1,rtk2,rtk3,rtk4
					!END IF
					!WRITE(*,*)	'rGCoef',rGCoef(itempseq),'rMTemp',rMTemp(itempseq)
				END DO
				rDeltaAbd=rXTemp(icntrk)-rXTemp(0)
				dTCatch(itempseq)=rCTemp(icntrk)
			ELSE
				rDeltaAbd=-rAbd(itempseq)
				dTCatch(itempseq)=0
			END IF
	
			!testing
			!IF(ISNAN(dTCatch(itempseq))) THEN
			!	WRITE(*,*)	itempseq, rIntR(iSpp),rDeltaAbd,rAbd,rK(itempseq,iSpp)
	
			!END IF
	
			!IF(ISNAN(rAbd(itempseq)) .OR. ISNAN(rDeltaAbd)) THEN
			!	WRITE(*,*) itempseq,rAbd(itempseq),rDeltaAbd,'rGCoef',rGCoef(itempseq),'rMTemp',rMTemp(itempseq)
			!	STOP
			!END IF
	
			rAbd(itempseq)=rAbd(itempseq)+rDeltaAbd
			IF(rAbd(itempseq)<0) THEN
				rAbd(itempseq)=0
			END IF
			IF(dTCatch(itempseq)<0) THEN
				dTCatch(itempseq)=0
			END IF
	
		END DO
	
	END SUBROUTINE CalDeltaAbd
	
	!!---------------- TRAVIS15 -------------!!
	!!---------------- H+ EFFECTS ON LARVAL MORTALITY HERE -------------!!
	SUBROUTINE CalLarDisp(rTempHQ10LarvMort,rTempOrigHSurf,rTempNewHSurf)
	!!---------------- TRAVIS15 -------------!!
		!Subroutine to calculate larval dispersal using ADR model, with U and V fields, and PLD
		USE DBEM, ONLY	:	U,V,nrow,ncol,rLarDiffCoef,rAbd,iLarTimeStep,rLarv,rLarSurv,rLarSettle
		USE DBEM,ONLY	:	rtDiffCoef
		IMPLICIT NONE
		INTEGER	::	iLarDur,icnt
		REAL	::	rdummyAbd(nrow*ncol),rTempLarv(nrow*ncol)
	
		!!---------------- TRAVIS16 -------------!!
		   !! ------------------------ Set variable for storing pH effects on larval mortality ------------------!!
		DOUBLE PRECISION,INTENT(IN)	:: rTempOrigHSurf,rTempNewHSurf   !! Read in as an array from original program (don't need "(i)")
		!! -------------- Set variables rTempHQ10LarvMort ---------------!!
		REAL,INTENT(IN) ::	rTempHQ10LarvMort
		!!---------------- pH Mortality multiplying factor; Temp variable for maintaining original LarSurv -------------!!
		REAL	::	rTHLarvMortfactor(nrow*ncol)
		REAL	::	rTLarSurv(nrow*ncol)
		!!---------------- TRAVIS16 -------------!!
	
		CHARACTER(12)	::	fname
		!Store the abundance of species i temporarily (for reading before exiting the subroutine)
		rdummyAbd=rAbd
	
	
		!Egg production (assuming an average GSI of 30%)
		rAbd=rAbd*0.3
	
		rLarv=0
		rTempLarv=0
	
		!Generate map segment
		CALL MapSeg
	
		!Read in larval diffcoef coefficient
		rtDiffCoef=rLarDiffCoef*60*60*24/(1000*1000)*iLarTimeStep
	
		!Calculate pelagic larval duration
		iLarDur=PLD(iLarTimeStep)
		CALL setabcdef
	
		!!---------------- TRAVIS17 -------------!!
		!!---------------- Put in pH effects on larval mortality. Algorithm here ---------!!
		!!---------------- Change rLarSurv to a temporary variable rTLarSurv ---------!!
		rTHLarvMortfactor=MAX(rTempHQ10LarvMort*rTempNewHSurf/rTempOrigHSurf+1-rTempHQ10LarvMort,1.0)
		rTLarSurv=MIN(rLarSurv*rTHLarvMortfactor,1.0)
	
		DO icnt=1,iLarDur
	
			CALL caltridag
	
			rTempLarv=rAbd*(1-rTLarSurv)*rLarSettle
	
			rAbd=rAbd*(1-rTLarSurv)-rTempLarv
	
			rLarv=rLarv+rTempLarv
	
			rTempLarv=0
		END DO
	
		rAbd=rdummyAbd
	
		!!---------------- TRAVIS17 -------------!!
	
	
		CONTAINS
		INTEGER FUNCTION PLD(itLarTimeStep)
			!A function to calculate pelagic larval duration
			USE DBEM,ONLY	:	rAbd,rSST,iSpp,pwater,nrow,ncol
			IMPLICIT NONE
			INTEGER,INTENT(IN)	::	itLarTimeStep
			REAL	:: rTempPLD(nrow*ncol),rPLD,rBo,rMeanInSST,rtemp,rtarea,rtempSST(nrow*ncol)
			REAL	::	rnum
			INTEGER	::	icnt
			CALL init_random_seed()
	
			CALL RANDOM_NUMBER(rnum)
	
			rtemp=0
			rtarea=0
			DO	icnt=1,nrow*ncol
				IF(rAbd(icnt)>0) THEN
					IF(rSST(icnt)>1) THEN
						rtemp=rtemp+LOG(rSST(icnt))*pwater(icnt)
					END IF
					rtarea=rtarea+pwater(icnt)
				END IF
				IF(rSST(icnt)>1) THEN
					rtempSST(icnt)=LOG(rSST(icnt)/15)
				ELSE
					rtempSST(icnt)=LOG((1.0/15.0))
				END IF
	
			END DO
			rMeanInSST=rtemp/rtarea
	
			rBo=0.739 + 0.739 * rMeanInSST
	
			rTempPLD=rBo-1.368*rtempSST-0.283*(rtempSST)**2
	
			rTempPLD=EXP(rTempPLD)*pwater
	
			rPLD=0
			DO icnt=1,nrow*ncol
				IF(rAbd(icnt)>0) THEN
					rPLD=rTempPLD(icnt)+rPLD
				END IF
			END DO
			rPLD=rPLD/rtarea
			IF(rnum<(1/5)) THEN
				rPLD=rPLD*2
			END IF
			rPLD=MIN(rPLD,365.0)
			rPLD=MAX(1.0,rPLD)
			PLD=CEILING(rPLD/itLarTimeStep)
			RETURN
		END FUNCTION PLD
	END SUBROUTINE CalLarDisp
	
	SUBROUTINE CreateDataFolder(cDatafname)
		!Subroutine to create folder in data directory
		USE DBEM,ONLY	:	cfolderpath
		IMPLICIT NONE
		CHARACTER(*),INTENT(IN)	::	cDatafname
	
		call system('mkdir -p out '//'/scratch/jepa/'//TRIM(cDatafname)) !JEPA change user name according to computecanada.ca
		!for cedar
		!call system('mkdir -p out '//'/home/wailung/projects/rrg-wailung/CMIP6/'//TRIM(cDatafname)) !JEPA change user name according to computecanada.ca
	
		!call system('mkdir -p out '//'/Users/wwlcheung/Fortran/Data/'//TRIM(cDatafname))
		RETURN
	END SUBROUTINE
	
	SUBROUTINE SalinityMap(istyear,iendyear,fname)
		!Subroutine to read salinity map, then convert ppm into salinity scale
		USE DBEM, ONLY	:	nrow,ncol,rSalPScale,pwater,cCCScenario
		IMPLICIT NONE
		INTEGER,INTENT(IN)	::	istyear,iendyear
		CHARACTER(*),INTENT(IN)	::	fname
		CHARACTER(LEN=4)	::	ctemp
		CHARACTER(LEN=40)	::	ctfname
		INTEGER	::	icnt
		REAL	::	rSalinity(nrow*ncol),rtSal(nrow*ncol)
	
		rtSal=0
		rSalinity=0
		DO icnt=istyear,iendyear
			Write(cTemp, '(i4)') icnt
			ctfname='Climate/'//TRIM(cCCScenario)//'/'//fname//TRIM(ctemp)//'.txt'
			CALL ReadCSV(ctfname,rtSal,nrow*ncol,0)
			rSalinity=rSalinity+rtSal/(iendyear-istyear+1)
		END DO
	
		rSalPScale=0
		DO icnt=1,nrow*ncol
			IF(pwater(icnt)>0) THEN
				IF (rSalinity(icnt)>=40) THEN	!Hyperhaline and metahaline
					rSalPScale(icnt,1)=1
				ELSEIF (rSalinity(icnt)>=29 .AND. rSalinity(icnt)<40) THEN	!mexoeuhaline
					rSalPScale(icnt,2)=1
				ELSEIF (rSalinity(icnt)>=18 .AND. rSalinity(icnt)<30) THEN	!polyhaline
					rSalPScale(icnt,3)=1
				ELSEIF (rSalinity(icnt)>=5 .AND. rSalinity(icnt)<18) THEN	!mesohaline
					rSalPScale(icnt,4)=1
				ELSEIF (rSalinity(icnt)<5) THEN	!oligohaline
					rSalPScale(icnt,5)=1
				END IF
			END IF
		END DO
		RETURN
	END SUBROUTINE SalinityMap
	
	!!---------------- TRAVIS15 ------------!!
	!!---------------- Add in variables rTempHQ10Mort to Subroutine line (Temp = Temporary) ----------------!!
	SUBROUTINE Ecophysiology(itSizeBinNum,rTWinf,rTVBonK,rTLwA,rTLwB,rTArrhCoef,rTAlloSCal,rTHQ10,rTempOrigST, &
	rTempOrigO2,rTempOrigH,rTempNewST,rTempNewO2,rTempNewH,rtKvar,rtlencl,rtWgcl,rTNewWinf,rTNewVBonK,rTNewWmat,rTNewNMort,rTNewLinf, &
	iTNewTmax,rTMeanW,rNMortTemp, &
	!!--------------- Added here ----------!!
	rTempHQ10Mort)
		!Subroutine to calculate changes in VBGF and natural mortality rate based on
		!Cheung et al. (2011) ICES JMS
		USE DBEM,ONLY	:	cfolderpath
		IMPLICIT NONE
		INTEGER,INTENT(IN)	::	itSizeBinNum
		REAL,INTENT(IN)	::	rTWinf,rTVBonK,rTLwA,rTLwB,rTArrhCoef(2),rTAlloScal(2),rtKvar
		REAL,INTENT(IN)	::	rtlencl(itSizeBinNum),rtWgcl(itSizeBinNum)
		REAL,INTENT(IN)	::	rTempOrigST,rTempNewST,rTHQ10
	
		!! -------------- Set variables rTempHQ10Mort ---------------!!
		REAL,INTENT(IN) ::	rTempHQ10Mort
		!!---------------- TRAVIS15 ------------!!
	
		DOUBLE PRECISION,INTENT(IN)	:: rTempOrigO2,rTempOrigH,rTempNewO2,rTempNewH
		REAL,INTENT(OUT)	::	rTNewWinf,rTNewVBonK,rTNewWmat,rTNewNMort(itSizeBinNum),rTNewLinf
		REAL,INTENT(OUT)	::	rNMortTemp
		DOUBLE PRECISION,INTENT(OUT)	::	rTMeanW
		INTEGER,INTENT(OUT)	::	iTNewTmax
	
		!!---------------- TRAVIS16 ------------!!
		!!---------------- create variable for pH effect on adult mortality factor ------------!!
		REAL	::	rTO2factor,rTHfactor,rTAnaCatCoef(2),rTHMfactor
		!!---------------- TRAVIS16 ------------!!
	
		DOUBLE PRECISION	::	rktemp1,rktemp2
		INTEGER	::	icnt,j
		REAL ::	rBlen(itSizeBinNum),rSSBlen(itSizeBinNum),rLen(itSizeBinNum),rLenTM(itSizeBinNum,itSizeBinNum)
		DOUBLE PRECISION	::	rLenFreq(0:40,itSizeBinNum)
		!Calculate initial parameters for growth model
		rktemp1=rtVBonK/(1-rTAlloscal(1))
	
		rTAnaCatCoef(2)=rktemp1/(EXP(-rTArrhCoef(2)/(273.15+rTempOrigST)*1000))
	
		rTAnaCatCoef(1)=rktemp1*rTWinf**(1-rTAlloscal(1))
		rTAnaCatCoef(1)=rTAnaCatCoef(1)/EXP(-rTArrhCoef(1)/(273.15+rTempOrigST)*1000)
	
		!Calculate coefficient for O2 and H+
		!!---------------- TRAVIS17 ------------!!
		!!---------------- Calculate coefficient for H+ mortality ------------!!
		rTO2factor=MIN(1.0,rTempNewO2/rTempOrigO2)
		rTHfactor=MAX(rTHQ10*rTempNewH/rTempOrigH + 1 - rTHQ10, 1.0)
		rTHMfactor=MAX(rTempHQ10Mort*rTempNewH/rTempOrigH + 1 - rTempHQ10Mort, 1.0)
		!!---------------- TRAVIS17 ------------!!
	
		!Calculate new Winf
		rTNewWinf=rTO2factor*rTAnaCatCoef(1)*EXP(-rTArrhCoef(1)/(273.15+rTempNewST)*1000)
	
	
		rktemp2=(rTAnaCatCoef(2)*rtHfactor*EXP(-rTArrhCoef(2)/(273.15+rTempNewST)*1000))
		rTNewWinf=(rTNewWinf/rktemp2)**(1/(1-rTAlloScal(1)))
	
		!Calculate new K
		rTNewVBonK=rTHfactor*rTAnaCatCoef(2)
		rktemp2=EXP(-rTArrhCoef(2)/(273.15+rTempNewSt)*1000)
		rktemp2=rktemp2*(1-rTAlloScal(1))
		rTNewVBonK=rTNewVBonK*rktemp2
	
		!Calculate new Wmat
		rTNewWmat=rTNewWinf*(0.714)**(1/(1-rTAlloScal(1)))
	
		!Calculate new Linf
		rTNewLinf=(rTNewWinf/rTLwA)**(1/rTLwB)
	
		!Calculate natural mortality rate
		DO icnt=1,itSizeBinNum
	
		!!---------------- TRAVIS18 ------------!!
			rTNewNMort(icnt)=EXP(4.355-0.085*LOG(rTNewWinf)+6.39*rTNewWinf/(rTNewWinf**3)+0.62*LOG(rTNewVBonK)-1190.43/(rTempNewSt+273.15))
			!rTNewNMort(icnt)=MIN(10.0,rTNewNMort(icnt))
	
			!!---------------- Add in pH effects on adult mortality ------------!!
			!!---------------- Replace above line with...: ------------!!
			rTNewNMort(icnt)=MIN(10.0,rTNewNMort(icnt)*rTHMfactor)
	
			!!---------------- TRAVIS18 ------------!!
	
			!IF(ISNAN(rTNewNMort(icnt)) .OR. rTNewNMort(icnt)>10000) THEN
			!	WRITE(*,*) rTNewWinf,rTempNewO2,rTempOrigO2,rTO2factor, rktemp2,rTNewVBonK,rTempNewSt
			!	STOP
			!END IF
		END DO
		rNMortTemp=rTNewNMort(1)
		!Calculate Tmax
		iTNewTmax=CEILING(3.0/rTNewVBonK)
		iTNewTMax=MIN(iTNewTMax,40)
	
		CALL LenTransMatrix(itSizeBinNum,rTNewLinf,rTNewVBonK,rtKVar,rtlencl,rLenTM)
	
		!Calculate per recruitment matrix
		CALL LMProjection(itSizeBinNum,iTNewTmax,rTNewWmat,rTNewNMort,rTLwA, &
						rTLwB,rLenTM,rtlencl,rLenFreq,rBlen,rSSBlen,rLen)
	
		!WRITE(*,*)	(rBlen(icnt),icnt=1,40)
	
		!Calculate Biomass per recruit
		rTMeanW=0.0
		rktemp2=0.0
		rktemp1=0.0
		DO icnt=1,itSizeBinNum
			!rktemp1=0.0
			!DO j=1,40
			!	rktemp1=rktemp1+rLenFreq(j,icnt)
			!	rTMeanW=rTMeanW+rLenFreq(j,icnt)
			!END DO
			!rktemp2=rktemp1*rtWgcl(icnt)+rktemp2
	
			!WRITE(*,*)	icnt
			!WRITE(*,*) rktemp1,rtWgcl(icnt),rktemp2
	
			rktemp1=rktemp1+rlen(icnt)
			rktemp2=rktemp2+rBlen(icnt)
		END DO
		IF(rktemp1>0) THEN
			rTMeanW=rktemp2/rktemp1 !/rTMeanW
		ELSE
			rTMeanW=0
		END IF
		!IF(ISNAN(rTMeanW)) THEN
		!	WRITE(*,*) 'rk',rktemp2,rktemp1
		!END IF
		!WRITE(*,*) itSizeBinNum,rTMeanW
	
					!Debugging
					!OPEN(UNIT=1,FILE=cfolderpath//'LenTM.csv')
					!	DO	icnt=1,itSizeBinNum
					!		WRITE(1,*)	(rLenTM(icnt,j),j=1,itSizeBinNum)
					!	END DO
					!CLOSE(1)
					!	WRITE(*,*)	(rNMortDummy(j),j=1,iSizeBinNum)
					!CALL LMProjection(iSizeBinNum,iTmaxDummy,rWmatDummy,rNMortDummy,rLwA(iSpp), &
					!	rLwB(iSpp),rLenTM,rlencl,rLenFreq,rBlen,rSSBlen,rLen)
	
					!Debugging
					!OPEN(UNIT=1,FILE=cfolderpath//'rlenfreq.csv')
					!DO	icnt=1,itSizeBinNum
					!WRITE(1,*)	(rBlen(icnt),icnt=0,40)
					!		WRITE(1,*)	(rLenFreq(j,icnt),j=1,40)
					!END DO
					!CLOSE(1)
	
					!WRITE(*,*)	'Mean weight = ', rTMeanW
	
		RETURN
	END SUBROUTINE Ecophysiology
	
	SUBROUTINE LenTransMatrix(itSizeBinNum,rTLinf,rTVBonK,rTKVar,rTlencl,rTLenTM)
		!Subroutine to set up the population length transition matrix
		IMPLICIT NONE
		INTEGER,INTENT(IN)	::	itSizeBinNum
		REAL,INTENT(IN)	::	rTLinf,rTVBonK,rTKVar,rTlencl(itSizeBinNum)
		REAL,INTENT(OUT)	::	rTLenTM(itSizeBinNum,itSizeBinNum)
		INTEGER	::	itrow,itcol
		REAL	::	rtemp1,rtemp2,rTpara1,rTpara2(itSizeBinNum,itSizeBinNum)
	
	
		rtemp1=0
		rTpara1=rTKvar/100*(rTLinf/5)
		Do itrow=1,itSizeBinNum
			DO itcol=1,itSizeBinNum
				rTpara2(itrow,itcol)=(rTlencl(itcol)-(rTLinf*(1-EXP(-rTVBonK))+rTlencl(itrow)*EXP(-rTVBonK)))**2
				rTpara2(itrow,itcol)=EXP(-1*rTpara2(itrow,itcol)/(2*(rTpara1**2)))
			END DO
		END DO
	
		DO itcol=1,itSizeBinNum
			rtemp1=0
			DO itrow=1,itSizeBinNum
				rtemp1=rTpara2(itrow,itcol)+rtemp1
			END DO
			rTLenTM(:,itcol)=rTpara2(:,itcol)/rtemp1
		END DO
		RETURN
	
	END SUBROUTINE LenTransMatrix
	
	SUBROUTINE LMProjection(itSizeBinSize,iTTmax,rTWmat,rTNMort,rTLwA,rTLwB,rTLenTM,rTlencl, &
	 rTLenFreq,rTBlen,rTSSBlen,rTLen)
		!Subroutine to calculate projection of abundance vs length per recruit
		IMPLICIT NONE
		INTEGER,INTENT(IN)	::	iTTmax,itSizeBinSize
		REAL,INTENT(IN)	::	rTWmat,rTNMort(itSizeBinSize),rTLwA,rTLwB
		REAL,INTENT(IN)	::	rTLenTM(itSizeBinSize,itSizeBinSize),rTlencl(itSizeBinSize)
		DOUBLE PRECISION,INTENT(OUT)	::	rTLenFreq(0:40,itSizeBinSize)
		REAL,INTENT(OUT)	::	rTBlen(itSizeBinSize),rTSSBlen(itSizeBinSize),rTLen(itSizeBinSize)
	
		INTEGER	::	iAge,itrow,itcol,icnt
		REAL	::	rNNext,rNNext2,rSurv,rtLm
	
		rNNext=0
		rTLenFreq=0
		rTLenFreq(0,1)=10.0
		DO iAge=1,iTTmax
			DO itcol=1,itSizeBinSize
				rNNext=0
				DO itrow=1,itSizeBinSize
					IF(rTLenTM(itrow,itcol)>0) THEN
						rSurv=EXP(-rTNMort(itrow))
						rNNext=rNNext+rTLenFreq(iAge-1,itrow)*rSurv*rTLenTM(itrow,itcol)
					END IF
				END DO
				rTLenFreq(iAge,itcol)=rNNext
			END DO
		END DO
	
		rTBlen=0
		rTSSBlen=0
		rTLen=0
	
		DO iAge=0,iTTmax
			DO	itcol=1,itSizeBinSize
				rTLen(itcol)=rTLen(itcol)+rTLenFreq(iAge,itcol)
				rTBlen(itcol)=rTBlen(itcol)+rTLenFreq(iAge,itcol)*rTLwA*rTlencl(itcol)**rTLwB
				IF (rTlencl(itcol)>=rTLm) THEN
					rTSSBlen(itcol)=rTSSBlen(itcol)+rTLenFreq(iAge,itcol)*rTLwA*rTlencl(itcol)**rTLwB
				END IF
			END DO
		END DO
		RETURN
	
	END SUBROUTINE LMProjection
	
	SUBROUTINE AllocCatch(rtempAbd,dltcatch)
		!Subroutine to allocate catch to cell that is proportional to the relative abundance in a cell
		USE DBEM,ONLY	:	nrow,ncol,dlResCatch
	
		IMPLICIT NONE
	
		REAL,INTENT(IN)	:: rtempAbd(nrow*ncol)
		DOUBLE PRECISION,INTENT(IN)	::	dltcatch
		INTEGER	::	icnt
		REAL	::	rTSum
		REAL	:: dlTempC(nrow*ncol)
		rTSum=SUM(rtempAbd)
	
		DO icnt=1,nrow*ncol
			dlTempC(icnt)=dltcatch*rtempAbd(icnt)/rTSum
		END DO
		dlResCatch=dlTempC
	END SUBROUTINE
	
	SUBROUTINE ReadCatch(itempyr,cTtaxon,dltACatch)
		!routine to read catch from a specific year
	
		USE DBEM,ONLY	:	cfolderpath
	
		IMPLICIT NONE
		INTEGER,INTENT(IN)	::	itempyr
		CHARACTER(*),INTENT(IN)	::	cTtaxon
		DOUBLE PRECISION,INTENT(OUT)	::	dltACatch
		CHARACTER(80)	::	ctemp3
		CHARACTER(6)	::	ctemp4
		INTEGER	::	icnt,icntyr
		DOUBLE PRECISION	::	dltemp
	
		icntyr=itempyr-1949
		ctemp4='Catch/'
		!WRITE(*,*) TRIM(cfolderpath)//TRIM(ctemp4)//TRIM(cTtaxon)//'.txt'
		OPEN(UNIT=4,file=TRIM(cfolderpath)//TRIM(ctemp4)//TRIM(cTtaxon)//'.txt')
			DO icnt=1,icntyr
				READ(4,*)	dltemp
			END DO
		CLOSE(4)
	
		WRITE(*,*) dltemp
	
		dltACatch=dltemp
	
		RETURN
	END SUBROUTINE
	
	subroutine init_random_seed()
		!Subroutine to randomize seed to generate random number
				USE IFPORT
				implicit none
				integer, allocatable :: iseed(:)
				integer :: i, n, un, istat, dt(8), t(2), s
				integer(4)	:: pid
				integer(8) :: count, tms
	
				call random_seed(size = n)
				allocate(iseed(n))
				! First try if the OS provides a random number generator
				open(newunit=un, file="/dev/urandom", access="stream", &
					 form="unformatted", action="read", status="old", iostat=istat)
				if (istat == 0) then
				   read(un) iseed
				   close(un)
				else
				   ! Fallback to XOR:ing the current time and pid. The PID is
				   ! useful in case one launches multiple instances of the same
				   ! program in parallel.
				   call system_clock(count)
				   if (count /= 0) then
					  t = transfer(count, t)
				   else
					  call date_and_time(values=dt)
					  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
						   + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
						   + dt(3) * 24 * 60 * 60 * 60 * 1000 &
						   + dt(5) * 60 * 60 * 1000 &
						   + dt(6) * 60 * 1000 + dt(7) * 1000 &
						   + dt(8)
					  t = transfer(tms, t)
				   end if
				   s = ieor(t(1), t(2))
				   pid = GETPID()
				   pid=pid + 1099279 ! Add a prime
				   s = ieor(s, pid)
				   if (n >= 3) then
					  iseed(1) = t(1) + 36269
					  iseed(2) = t(2) + 72551
					  iseed(3) = pid
					  if (n > 3) then
						 iseed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
					  end if
				   else
					  iseed = s + 37 * (/ (i, i = 0, n - 1 ) /)
				   end if
				end if
				call random_seed(put=iseed)
	end subroutine init_random_seed
	