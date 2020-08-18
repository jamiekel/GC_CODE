!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: PAH_MOD
!
! !DESCRIPTION: Module pah_mod contains variables and routines for the 
!  GEOS-Chem "online" PAH simulation. (PDI) 
!\\   
!\\   
! !INTERFACE: 
!
MODULE PAH_MOD
! 
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: PAHCHEM
!
! !DEFINED PARAMETERS:
!
  REAL(fp), PARAMETER   :: SMALLNUM   = 1e-20_fp
!
! !LOCAL VARIABLES:
!
  INTEGER               :: AS

  !Scalers
  REAL(fp), ALLOCATABLE :: BCCONV(:,:,:)
  REAL(fp), ALLOCATABLE :: OCCONV(:,:,:)
	  
  ! Species ID flags
  INTEGER    :: id_ACE,   id_ACEBCPO,   id_ACEBCPI,   id_ACEOCPO,   id_ACEOCPI
  INTEGER    :: id_ACY,   id_ACYBCPO,   id_ACYBCPI,   id_ACYOCPO,   id_ACYOCPI
  INTEGER    :: id_ANT,   id_ANTBCPO,   id_ANTBCPI,   id_ANTOCPO,   id_ANTOCPI
  INTEGER    :: id_BAA,   id_BAABCPO,   id_BAABCPI,   id_BAAOCPO,   id_BAAOCPI
  INTEGER    :: id_BAP,   id_BAPBCPO,   id_BAPBCPI,   id_BAPOCPO,   id_BAPOCPI
  INTEGER    :: id_BAPB,  id_BAPBBCPO,  id_BAPBBCPI,  id_BAPBOCPO,  id_BAPBOCPI
  INTEGER    :: id_BAPC,  id_BAPCBCPO,  id_BAPCBCPI,  id_BAPCOCPO,  id_BAPCOCPI
  INTEGER    :: id_BBF,   id_BBFBCPO,   id_BBFBCPI,   id_BBFOCPO,   id_BBFOCPI
  INTEGER    :: id_BHP,   id_BHPBCPO,   id_BHPBCPI,   id_BHPOCPO,   id_BHPOCPI
  INTEGER    :: id_BKF,   id_BKFBCPO,   id_BKFBCPI,   id_BKFOCPO,   id_BKFOCPI
  INTEGER    :: id_CHR,   id_CHRBCPO,   id_CHRBCPI,   id_CHROCPO,   id_CHROCPI
  INTEGER    :: id_DAHA,  id_DAHABCPO,  id_DAHABCPI,  id_DAHAOCPO,  id_DAHAOCPI
  INTEGER    :: id_FLA,   id_FLABCPO,   id_FLABCPI,   id_FLAOCPO,   id_FLAOCPI
  INTEGER    :: id_FLO,   id_FLOBCPO,   id_FLOBCPI,   id_FLOOCPO,   id_FLOOCPI
  INTEGER    :: id_ICDP,  id_ICDPBCPO,  id_ICDPBCPI,  id_ICDPOCPO,  id_ICDPOCPI
  INTEGER    :: id_NAP,   id_NAPBCPO,   id_NAPBCPI,   id_NAPOCPO,   id_NAPOCPI
  INTEGER    :: id_PHE,   id_PHEBCPO,   id_PHEBCPI,   id_PHEOCPO,   id_PHEOCPI
  INTEGER    :: id_PYR,   id_PYRBCPO,   id_PYRBCPI,   id_PYROCPO,   id_PYROCPI
  INTEGER    :: id_BCPI,  id_BCPO,      id_OCPI,      id_OCPO
  
CONTAINS
!EOP
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PAHCHEM
!
! !DESCRIPTION: Subroutine PAHCHEM handles the partitioning of PAH and the
!  hydrophobic to hydrophlic conversion of the PAH on particulate tracers!.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PAHCHEM( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD 
    USE ErrCode_Mod 
    USE ERROR_MOD,          ONLY : DEBUG_MSG, ERROR_STOP, ALLOC_ERR
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
	USE State_Chm_Mod,      ONLY : Ind_
    USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!	
    LOGICAL, SAVE      :: FIRSTCHEM = .TRUE.
    LOGICAL            :: prtDebug
    LOGICAL            :: LPRT
	LOGICAL            :: PPFER_PARTITIONING

    ! Pointers
    REAL(fp), POINTER  :: Spc(:,:,:,:)  ! Species [Kg]
	REAL(fp), POINTER  :: T(:,:,:)      ! Temperature [K]
	REAL(fp), POINTER ::  AIRVOL(:,:,:) ! Grid box volume [m3] (dry air)

    ! For getting fields from HEMCO 
    LOGICAL            :: aIR
	
	!Define species IDs
	id_ACE       = IND_('ACE' )
	id_ACEBCPI   = IND_('ACEBCPI' )
    id_ACEBCPO   = IND_('ACEBCPO' )
	id_ACEOCPI   = IND_('ACEOCPI' )
    id_ACEOCPO   = IND_('ACEOCPO' )
	id_ACY       = IND_('ACY' )
	id_ACYBCPI   = IND_('ACYBCPI' )
    id_ACYBCPO   = IND_('ACYBCPO' )
	id_ACYOCPI   = IND_('ACYOCPI' )
    id_ACYOCPO   = IND_('ACYOCPO' )
	id_ANT       = IND_('ANT' )
	id_ANTBCPI   = IND_('ANTBCPI' )
    id_ANTBCPO   = IND_('ANTBCPO' )
	id_ANTOCPI   = IND_('ANTOCPI' )
    id_ANTOCPO   = IND_('ANTOCPO' )
	id_BAA       = IND_('BAA' )
	id_BAABCPI   = IND_('BAABCPI' )
    id_BAABCPO   = IND_('BAABCPO' )
	id_BAAOCPI   = IND_('BAAOCPI' )
    id_BAAOCPO   = IND_('BAAOCPO' )
	id_BAP       = IND_('BAP' )
	id_BAPBCPI   = IND_('BAPBCPI' )
    id_BAPBCPO   = IND_('BAPBCPO' )
	id_BAPOCPI   = IND_('BAPOCPI' )
    id_BAPOCPO   = IND_('BAPOCPO' )
    id_BAPB       = IND_('BAPB' )
    id_BAPBBCPI   = IND_('BAPBBCPI' )
    id_BAPBBCPO   = IND_('BAPBBCPO' )
    id_BAPBOCPI   = IND_('BAPBOCPI' )
    id_BAPBOCPO   = IND_('BAPBOCPO' )
    id_BAPC       = IND_('BAPC' )
    id_BAPCBCPI   = IND_('BAPCBCPI' )
    id_BAPCBCPO   = IND_('BAPCBCPO' )
    id_BAPCOCPI   = IND_('BAPCOCPI' )
    id_BAPCOCPO   = IND_('BAPCOCPO' )
	id_BBF       = IND_('BBF' )
	id_BBFBCPI   = IND_('BBFBCPI' )
    id_BBFBCPO   = IND_('BBFBCPO' )
	id_BBFOCPI   = IND_('BBFOCPI' )
    id_BBFOCPO   = IND_('BBFOCPO' )
	id_BHP       = IND_('BHP' )
	id_BHPBCPI   = IND_('BHPBCPI' )
    id_BHPBCPO   = IND_('BHPBCPO' )
	id_BHPOCPI   = IND_('BHPOCPI' )
    id_BHPOCPO   = IND_('BHPOCPO' )
	id_BKF       = IND_('BKF' )
	id_BKFBCPI   = IND_('BKFBCPI' )
    id_BKFBCPO   = IND_('BKFBCPO' )
	id_BKFOCPI   = IND_('BKFOCPI' )
    id_BKFOCPO   = IND_('BKFOCPO' )
	id_CHR       = IND_('CHR' )
	id_CHRBCPI   = IND_('CHRBCPI' )
    id_CHRBCPO   = IND_('CHRBCPO' )
	id_CHROCPI   = IND_('CHROCPI' )
    id_CHROCPO   = IND_('CHROCPO' )
	id_DAHA      = IND_('DAHA' )
	id_DAHABCPI  = IND_('DAHABCPI' )
    id_DAHABCPO  = IND_('DAHABCPO' )
	id_DAHAOCPI  = IND_('DAHAOCPI' )
    id_DAHAOCPO  = IND_('DAHAOCPO' )
	id_FLA       = IND_('FLA' )
	id_FLABCPI   = IND_('FLABCPI' )
    id_FLABCPO   = IND_('FLABCPO' )
	id_FLAOCPI   = IND_('FLAOCPI' )
    id_FLAOCPO   = IND_('FLAOCPO' )
	id_FLO       = IND_('FLO' )
	id_FLOBCPI   = IND_('FLOBCPI' )
    id_FLOBCPO   = IND_('FLOBCPO' )
	id_FLOOCPI   = IND_('FLOOCPI' )
    id_FLOOCPO   = IND_('FLOOCPO' )
	id_ICDP      = IND_('ICDP' )
	id_ICDPBCPI  = IND_('ICDPBCPI' )
    id_ICDPBCPO  = IND_('ICDPBCPO' )
	id_ICDPOCPI  = IND_('ICDPOCPI' )
    id_ICDPOCPO  = IND_('ICDPOCPO' )
	id_NAP       = IND_('NAP' )
	id_NAPBCPI   = IND_('NAPBCPI' )
    id_NAPBCPO   = IND_('NAPBCPO' )
	id_NAPOCPI   = IND_('NAPOCPI' )
    id_NAPOCPO   = IND_('NAPOCPO' )
	id_PHE       = IND_('PHE' )
	id_PHEBCPI   = IND_('PHEBCPI' )
    id_PHEBCPO   = IND_('PHEBCPO' )
	id_PHEOCPI   = IND_('PHEOCPI' )
    id_PHEOCPO   = IND_('PHEOCPO' )
	id_PYR       = IND_('PYR' )
	id_PYRBCPI   = IND_('PYRBCPI' )
    id_PYRBCPO   = IND_('PYRBCPO' )
	id_PYROCPI   = IND_('PYROCPI' )
    id_PYROCPO   = IND_('PYROCPO' )
    id_BCPO      = IND_('BCPO' )
    id_BCPO      = IND_('BCPI' )
    id_OCPO      = IND_('OCPO' )
    id_OCPI      = IND_('OCPI' )

    WRITE(6,*) ('JMK - Checking new PAHs have correct indexes')
    WRITE(6,*) ('id_BAP =', id_BAP)
    WRITE(6,*) ('id_BAPBCPI =', id_BAPBCPI)
    WRITE(6,*) ('id_BAPBCPO =', id_BAPBCPO)
    WRITE(6,*) ('id_BAPOCPI =', id_BAPOCPI)
    WRITE(6,*) ('id_BAPOCPO =', id_BAPOCPO)
    WRITE(6,*) ('id_BAPB =', id_BAPB)
    WRITE(6,*) ('id_BAPBBCPI =', id_BAPBBCPI)
    WRITE(6,*) ('id_BAPBBCPO =', id_BAPBBCPO)
    WRITE(6,*) ('id_BAPBOCPI =', id_BAPBOCPI)
    WRITE(6,*) ('id_BAPBOCPO =', id_BAPBOCPO)
    WRITE(6,*) ('id_BAPC =', id_BAPC)
    WRITE(6,*) ('id_BAPCBCPI =', id_BAPCBCPI)
    WRITE(6,*) ('id_BAPBCPO  =', id_BAPCBCPO)
    WRITE(6,*) ('id_BAPCOCPI =', id_BAPCOCPI)
    WRITE(6,*) ('id_BAPCOCPO =', id_BAPCOCPO)   

	ALLOCATE( BCCONV( IIPAR, JJPAR, LLPAR ))
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCCONV' )
    BCCONV = 0e+0_fp
	
	ALLOCATE( OCCONV( IIPAR, JJPAR, LLPAR ))
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCCONV' )
    OCCONV = 0e+0_fp

    !=================================================================
    ! CHEMBAPCARBON begins here!
    !=================================================================

    ! Assume success
    RC                   = GC_SUCCESS

    ! am I root? 
    aIR                  = am_I_Root

    ! Copy fields from INPUT_OPT to local variables for use below
    LPRT                 = Input_Opt%LPRT
	PPFER_PARTITIONING   = Input_Opt%PPFER_PARTITIONING

    ! Do we have to print debug output?
    prtDebug             = ( LPRT .and. am_I_Root )

    ! Point to chemical species array [kg]
    Spc                  => State_Chm%Species
	
	! Point to Temperature array
	T                    => State_Met%T
	
	! Point to AIRVOL array
	AIRVOL               => State_Met%AIRVOL

    !=================================================================
    ! Do chemistry for PAH on aerosol species 
    !=================================================================
	
	! ACE
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_ACEBCPO), Spc(:,:,:,id_ACEBCPI),   &
					Spc(:,:,:,id_ACEOCPO), Spc(:,:,:,id_ACEOCPI), RC )
					
	! ACY
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_ACYBCPO), Spc(:,:,:,id_ACYBCPI),   &
					Spc(:,:,:,id_ACYOCPO), Spc(:,:,:,id_ACYOCPI), RC )
					
	! ANT
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_ANTBCPO), Spc(:,:,:,id_ANTBCPI),   &
					Spc(:,:,:,id_ANTOCPO), Spc(:,:,:,id_ANTOCPI), RC )
	
    ! BAA
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_BAABCPO), Spc(:,:,:,id_BAABCPI),   &
					Spc(:,:,:,id_BAAOCPO), Spc(:,:,:,id_BAAOCPI), RC )
	
    ! BAP
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_BAPBCPO), Spc(:,:,:,id_BAPBCPI),   &
					Spc(:,:,:,id_BAPOCPO), Spc(:,:,:,id_BAPOCPI), RC )
	
    ! BAPB
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_BAPBBCPO), Spc(:,:,:,id_BAPBBCPI),   &
                                        Spc(:,:,:,id_BAPBOCPO), Spc(:,:,:,id_BAPBOCPI), RC )
    ! BAPC
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_BAPCBCPO), Spc(:,:,:,id_BAPCBCPI),   &
                                        Spc(:,:,:,id_BAPCOCPO), Spc(:,:,:,id_BAPCOCPI), RC )
    ! BBF
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_BBFBCPO), Spc(:,:,:,id_BBFBCPI),   &
					Spc(:,:,:,id_BBFOCPO), Spc(:,:,:,id_BBFOCPI), RC )
	
    ! BHP
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_BHPBCPO), Spc(:,:,:,id_BHPBCPI),   &
					Spc(:,:,:,id_BHPOCPO), Spc(:,:,:,id_BHPOCPI), RC )
	
    ! BKF
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_BKFBCPO), Spc(:,:,:,id_BKFBCPI),   &
					Spc(:,:,:,id_BKFOCPO), Spc(:,:,:,id_BKFOCPI), RC )
	
    ! CHR
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_CHRBCPO), Spc(:,:,:,id_CHRBCPI),   &
					Spc(:,:,:,id_CHROCPO), Spc(:,:,:,id_CHROCPI), RC )
					
    ! DAHA
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_DAHABCPO), Spc(:,:,:,id_DAHABCPI), &
					Spc(:,:,:,id_DAHAOCPO), Spc(:,:,:,id_DAHAOCPI), RC )
	
    ! FLA
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_FLABCPO), Spc(:,:,:,id_FLABCPI),   &
					Spc(:,:,:,id_FLAOCPO), Spc(:,:,:,id_FLAOCPI), RC )
	
    ! FLO
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_FLOBCPO), Spc(:,:,:,id_FLOBCPI),   &
					Spc(:,:,:,id_FLOOCPO), Spc(:,:,:,id_FLOOCPI), RC )
	
    ! ICDP
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_ICDPBCPO), Spc(:,:,:,id_ICDPBCPI), &
					Spc(:,:,:,id_ICDPOCPO), Spc(:,:,:,id_ICDPOCPI), RC )					
	
    ! NAP (Unstable, due to its low paritioning may be worth removing NAP on PM)
    ! CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_NAPBCPO), Spc(:,:,:,id_NAPBCPI), &
	!				Spc(:,:,:,id_NAPOCPO), Spc(:,:,:,id_NAPOCPI), RC )	
	
    ! PHE
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_PHEBCPO), Spc(:,:,:,id_PHEBCPI),   &
					Spc(:,:,:,id_PHEOCPO), Spc(:,:,:,id_PHEOCPI), RC )
	
    ! PYR
    CALL PAH_BC_OC( am_I_Root, Input_Opt, Spc(:,:,:,id_PYRBCPO), Spc(:,:,:,id_PYRBCPI),   &
					Spc(:,:,:,id_PYROCPO), Spc(:,:,:,id_PYROCPI), RC )					
					
	IF ( prtDebug ) THEN
	  CALL DEBUG_MSG( '### AFTER PAH DEP and PO-PI conversion ' )
	ENDIF
	!=================================================================
    ! Do partitioning
    !=================================================================
	
	IF (.NOT. PPFER_PARTITIONING) THEN
	  ! Original Dachs-Eisenreich scheme (SP - single parameter scheme)
	  ! ACE
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "ACE ",                &
							 Spc(:,:,:,id_ACE),      Spc(:,:,:,id_ACEBCPO),  Spc(:,:,:,id_ACEBCPI), &
							 Spc(:,:,:,id_ACEOCPO),  Spc(:,:,:,id_ACEOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! ACY
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "ACY ",                &
							 Spc(:,:,:,id_ACY),      Spc(:,:,:,id_ACYBCPO),  Spc(:,:,:,id_ACYBCPI), &
							 Spc(:,:,:,id_ACYOCPO),  Spc(:,:,:,id_ACYOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! ANT
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "ANT ",                &
							 Spc(:,:,:,id_ANT),      Spc(:,:,:,id_ANTBCPO),  Spc(:,:,:,id_ANTBCPI), &
							 Spc(:,:,:,id_ANTOCPO),  Spc(:,:,:,id_ANTOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! BAA
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BAA  ",               &
							 Spc(:,:,:,id_BAA),      Spc(:,:,:,id_BAABCPO),  Spc(:,:,:,id_BAABCPI), &
							 Spc(:,:,:,id_BAAOCPO),  Spc(:,:,:,id_BAAOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! BAP
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BAP  ",               &
							 Spc(:,:,:,id_BAP),      Spc(:,:,:,id_BAPBCPO),  Spc(:,:,:,id_BAPBCPI), &
							 Spc(:,:,:,id_BAPOCPO),  Spc(:,:,:,id_BAPOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! BAPB
          CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BAPB ",               &
                                                         Spc(:,:,:,id_BAPB),      Spc(:,:,:,id_BAPBBCPO),  Spc(:,:,:,id_BAPBBCPI), &
                                                         Spc(:,:,:,id_BAPBOCPO),  Spc(:,:,:,id_BAPBOCPI),  Spc(:,:,:,id_BCPO),    &
                                                         Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &                                               
                                                         AIRVOL(:,:,:),              T(:,:,:),               RC )
          ! BAPC
          CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BAPC ",               &
                                                         Spc(:,:,:,id_BAPC),      Spc(:,:,:,id_BAPCBCPO),  Spc(:,:,:,id_BAPCBCPI), &
                                                         Spc(:,:,:,id_BAPCOCPO),  Spc(:,:,:,id_BAPCOCPI),  Spc(:,:,:,id_BCPO),    &
                                                         Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &
                                                         AIRVOL(:,:,:),              T(:,:,:),               RC )
          ! BBF
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BBF  ",               &
							 Spc(:,:,:,id_BBF),      Spc(:,:,:,id_BBFBCPO),  Spc(:,:,:,id_BBFBCPI), &
							 Spc(:,:,:,id_BBFOCPO),  Spc(:,:,:,id_BBFOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! BHP
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BHP ",                &
							 Spc(:,:,:,id_BHP),      Spc(:,:,:,id_BHPBCPO),  Spc(:,:,:,id_BHPBCPI), &
							 Spc(:,:,:,id_BHPOCPO),  Spc(:,:,:,id_BHPOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! BKF
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BKF ",                &
							 Spc(:,:,:,id_BKF),      Spc(:,:,:,id_BKFBCPO),  Spc(:,:,:,id_BKFBCPI), &
							 Spc(:,:,:,id_BKFOCPO),  Spc(:,:,:,id_BKFOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! CHR
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "CHR ",                &
							 Spc(:,:,:,id_CHR),      Spc(:,:,:,id_CHRBCPO),  Spc(:,:,:,id_CHRBCPI), &
							 Spc(:,:,:,id_CHROCPO),  Spc(:,:,:,id_CHROCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! DAHA
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "DAHA",                &
							 Spc(:,:,:,id_DAHA),     Spc(:,:,:,id_DAHABCPO), Spc(:,:,:,id_DAHABCPI),&
							 Spc(:,:,:,id_DAHAOCPO), Spc(:,:,:,id_DAHAOCPI), Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! FLA
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "FLA ",                &
							 Spc(:,:,:,id_FLA),      Spc(:,:,:,id_FLABCPO),  Spc(:,:,:,id_FLABCPI), &
							 Spc(:,:,:,id_FLAOCPO),  Spc(:,:,:,id_FLAOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! FLO
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "FLO ",                &
							 Spc(:,:,:,id_FLO),      Spc(:,:,:,id_FLOBCPO),  Spc(:,:,:,id_FLOBCPI), &
							 Spc(:,:,:,id_FLOOCPO),  Spc(:,:,:,id_FLOOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	  ! ICDP
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "ICDP",                &
							 Spc(:,:,:,id_ICDP),     Spc(:,:,:,id_ICDPBCPO), Spc(:,:,:,id_ICDPBCPI),&
							 Spc(:,:,:,id_ICDPOCPO), Spc(:,:,:,id_ICDPOCPI), Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )
	
	  ! NAP (Unstable, due to its low paritioning may be worth removing NAP on PM)
	  ! CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "NAP",                 &
							 ! Spc(:,:,:,id_ICDP),     Spc(:,:,:,id_ICDPBCPO), Spc(:,:,:,id_ICDPBCPI),&
							 ! Spc(:,:,:,id_ICDPOCPO), Spc(:,:,:,id_ICDPOCPI), Spc(:,:,:,id_BCPO),    &
							 ! Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 ! AIRVOL(:,:,:),		     T(:,:,:),               RC )

	  ! PHE
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "PHE ",                &
							 Spc(:,:,:,id_PHE),      Spc(:,:,:,id_PHEBCPO),  Spc(:,:,:,id_PHEBCPI), &
							 Spc(:,:,:,id_PHEOCPO),  Spc(:,:,:,id_PHEOCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )

	  ! PYR
	  CALL PARTITIONING_PAH( am_I_Root,              Input_Opt,              "PYR ",                &
							 Spc(:,:,:,id_PYR),      Spc(:,:,:,id_PYRBCPO),  Spc(:,:,:,id_PYRBCPI), &
							 Spc(:,:,:,id_PYROCPO),  Spc(:,:,:,id_PYROCPI),  Spc(:,:,:,id_BCPO),    &
							 Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
							 AIRVOL(:,:,:),		     T(:,:,:),               RC )		
						   
      IF ( prtDebug ) THEN
        CALL DEBUG_MSG( '### SP PARTITIONING' )
	  ENDIF
	ENDIF
	 
	IF (PPFER_PARTITIONING) THEN
	  ! PPFER PARTITIONING SCHEME https://www.sciencedirect.com/science/article/pii/S0021967314012370
	  ! ACE
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "ACE ",                &
	                               Spc(:,:,:,id_ACE),      Spc(:,:,:,id_ACEBCPO),  Spc(:,:,:,id_ACEBCPI), &
                                   Spc(:,:,:,id_ACEOCPO),  Spc(:,:,:,id_ACEOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! ACY
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "ACY ",                &
	                               Spc(:,:,:,id_ACY),      Spc(:,:,:,id_ACYBCPO),  Spc(:,:,:,id_ACYBCPI), &
                                   Spc(:,:,:,id_ACYOCPO),  Spc(:,:,:,id_ACYOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! ANT
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "ANT ",                &
	                               Spc(:,:,:,id_ANT),      Spc(:,:,:,id_ANTBCPO),  Spc(:,:,:,id_ANTBCPI), &
                                   Spc(:,:,:,id_ANTOCPO),  Spc(:,:,:,id_ANTOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! BAA
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BAA  ",               &
	                               Spc(:,:,:,id_BAA),      Spc(:,:,:,id_BAABCPO),  Spc(:,:,:,id_BAABCPI), &
                                   Spc(:,:,:,id_BAAOCPO),  Spc(:,:,:,id_BAAOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! BAP
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BAP  ",               &
	                               Spc(:,:,:,id_BAP),      Spc(:,:,:,id_BAPBCPO),  Spc(:,:,:,id_BAPBCPI), &
                                   Spc(:,:,:,id_BAPOCPO),  Spc(:,:,:,id_BAPOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! BAPB
          CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BAPC ",               &
                                       Spc(:,:,:,id_BAPB),      Spc(:,:,:,id_BAPBBCPO),  Spc(:,:,:,id_BAPBBCPI), &
                                   Spc(:,:,:,id_BAPBOCPO),  Spc(:,:,:,id_BAPBOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &                                             
                                                           AIRVOL(:,:,:),                  T(:,:,:),               RC )
          ! BAPC
          CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BAPC ",               &
                                       Spc(:,:,:,id_BAPC),      Spc(:,:,:,id_BAPCBCPO),  Spc(:,:,:,id_BAPCBCPI), &
                                   Spc(:,:,:,id_BAPCOCPO),  Spc(:,:,:,id_BAPCOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &
                                                           AIRVOL(:,:,:),                  T(:,:,:),               RC )
          ! BBF
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BBF  ",               &
	                               Spc(:,:,:,id_BBF),      Spc(:,:,:,id_BBFBCPO),  Spc(:,:,:,id_BBFBCPI), &
                                   Spc(:,:,:,id_BBFOCPO),  Spc(:,:,:,id_BBFOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! BHP
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BHP ",                &
	                               Spc(:,:,:,id_BHP),      Spc(:,:,:,id_BHPBCPO),  Spc(:,:,:,id_BHPBCPI), &
                                   Spc(:,:,:,id_BHPOCPO),  Spc(:,:,:,id_BHPOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! BKF
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "BKF ",                &
	                               Spc(:,:,:,id_BKF),      Spc(:,:,:,id_BKFBCPO),  Spc(:,:,:,id_BKFBCPI), &
                                   Spc(:,:,:,id_BKFOCPO),  Spc(:,:,:,id_BKFOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! CHR
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "CHR ",                &
	                               Spc(:,:,:,id_CHR),      Spc(:,:,:,id_CHRBCPO),  Spc(:,:,:,id_CHRBCPI), &
                                   Spc(:,:,:,id_CHROCPO),  Spc(:,:,:,id_CHROCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! DAHA
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "DAHA",                &
	                               Spc(:,:,:,id_DAHA),     Spc(:,:,:,id_DAHABCPO), Spc(:,:,:,id_DAHABCPI),&
                                   Spc(:,:,:,id_DAHAOCPO), Spc(:,:,:,id_DAHAOCPI), Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! FLA
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "FLA ",                &
	                               Spc(:,:,:,id_FLA),      Spc(:,:,:,id_FLABCPO),  Spc(:,:,:,id_FLABCPI), &
                                   Spc(:,:,:,id_FLAOCPO),  Spc(:,:,:,id_FLAOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! FLO
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "FLO ",                &
	                               Spc(:,:,:,id_FLO),      Spc(:,:,:,id_FLOBCPO),  Spc(:,:,:,id_FLOBCPI), &
                                   Spc(:,:,:,id_FLOOCPO),  Spc(:,:,:,id_FLOOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	  ! ICDP
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "ICDP",                &
	                               Spc(:,:,:,id_ICDP),     Spc(:,:,:,id_ICDPBCPO), Spc(:,:,:,id_ICDPBCPI),&
                                   Spc(:,:,:,id_ICDPOCPO), Spc(:,:,:,id_ICDPOCPI), Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )
	 
	  ! NAP (Unstable, due to its low paritioning may be worth removing NAP on PM)
	  ! CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "NAP",                 &
	                               ! Spc(:,:,:,id_ICDP),     Spc(:,:,:,id_ICDPBCPO), Spc(:,:,:,id_ICDPBCPI),&
                                   ! Spc(:,:,:,id_ICDPOCPO), Spc(:,:,:,id_ICDPOCPI), Spc(:,:,:,id_BCPO),    &
                                   ! Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           ! AIRVOL(:,:,:),		   T(:,:,:),               RC )
								   
	  ! PHE
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "PHE ",                &
	                               Spc(:,:,:,id_PHE),      Spc(:,:,:,id_PHEBCPO),  Spc(:,:,:,id_PHEBCPI), &
                                   Spc(:,:,:,id_PHEOCPO),  Spc(:,:,:,id_PHEOCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )

	  ! PYR
	  CALL PPFER_PARTITIONING_PAH( am_I_Root,              Input_Opt,              "PYR ",                &
	                               Spc(:,:,:,id_PYR),      Spc(:,:,:,id_PYRBCPO),  Spc(:,:,:,id_PYRBCPI), &
                                   Spc(:,:,:,id_PYROCPO),  Spc(:,:,:,id_PYROCPI),  Spc(:,:,:,id_BCPO),    &
                                   Spc(:,:,:,id_BCPI),     Spc(:,:,:,id_OCPO),     Spc(:,:,:,id_OCPI),    &						   
						           AIRVOL(:,:,:),		   T(:,:,:),               RC )		
      IF ( prtDebug ) THEN
        CALL DEBUG_MSG( '### PPFER PARTITIONING' )
	  ENDIF
	ENDIF
	
    ! Free pointer
    Spc    => NULL()
	T      => NULL()
	AIRVOL => NULL()
	
	!Clean Up
	IF ( ALLOCATED( BCCONV) ) DEALLOCATE( BCCONV)
	IF ( ALLOCATED( OCCONV) ) DEALLOCATE( OCCONV) 
	
  END SUBROUTINE PAHCHEM 
!EOC	  
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PAH_BC_OC
!
! !DESCRIPTION: Subroutine PAH_BC_OC calls the PO to PI conversion and 
!               deposition for PAH on particulate
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PAH_BC_OC( am_I_Root, Input_Opt, TC1, TC2, TC3, TC4, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root             ! Root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
!
! !INPUT/OUTPUT PARAMETERS:
! 
    REAL(fp),       INTENT(INOUT) :: TC1(IIPAR,JJPAR,LLPAR) ! PAHBCPO [kg]
    REAL(fp),       INTENT(INOUT) :: TC2(IIPAR,JJPAR,LLPAR) ! PAHBCPI [kg]
    REAL(fp),       INTENT(INOUT) :: TC3(IIPAR,JJPAR,LLPAR) ! PAHOCPO [kg]
    REAL(fp),       INTENT(INOUT) :: TC4(IIPAR,JJPAR,LLPAR) ! PAHOCPI [kg]
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!


    !=================================================================
    ! CHEM_BCPO begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
	
	! Chemistry for hydrophobic PAHBCPO
    !BAP
	IF ( id_BAPBCPO > 0 ) THEN
      CALL CHEM_PAHBCPO( am_I_Root, Input_Opt, TC1, RC )
    ENDIF
	  
	! AS THINGS CURRENTLY STAND CHEM_PAHBCPO MUST BE FOLLOWED BY CHEM_PAHBCPI
	! FOR EACH SPECIES. AS THE SECOND USES THE BCCONV ARRAY CREATED BY THE
	! FIRST. THESE TWO FUCTIONS SHOULD BE COMBINED INTO ONE IN THE NEAR FUTURE.

    ! Chemistry for hydrophilic PAHBCPI
	!BAP
    IF ( id_BAPBCPI > 0 ) THEN
      CALL CHEM_PAHBCPI( am_I_Root, Input_Opt, TC2, RC )
    ENDIF

	! Chemistry for hydrophobic PAHOCPO
	!BAP
    IF ( id_BAPOCPO > 0 ) THEN
      CALL CHEM_PAHOCPO( am_I_Root, Input_Opt, TC3, RC )
    ENDIF
	  
    ! Chemistry for hydrophilic PAHOCPI
	!BAP
    IF ( id_BAPOCPI > 0 ) THEN
      CALL CHEM_PAHOCPI( am_I_Root, Input_Opt, TC4, RC )
    ENDIF


  END SUBROUTINE PAH_BC_OC
!EOC	  
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CHEM_PAHBCPO
!
! !DESCRIPTION: Subroutine CHEM\_BCPO converts hydrophobic PAHBC to hydrophilic 
!  PAHBC and calculates the dry deposition of hydrophobic PAHBC.
!  COPIED and EDITED chem_bcpo 
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_PAHBCPO( am_I_Root, Input_Opt, TC, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root             ! Root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
!
! !INPUT/OUTPUT PARAMETERS:
! 
    REAL(fp),       INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! H-phobic BC [kg]
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: I,      J,   L
    REAL(fp)              :: DTCHEM, KBC, FREQ, TC0, CNEW, RKT
	
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: BC_LIFE = 1.15e+0_fp

    !=================================================================
    ! CHEM_BCPO begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Initialize
    KBC       = 1.e+0_fp / ( 86400e+0_fp * BC_LIFE )
    DTCHEM    = GET_TS_CHEM() * 60e+0_fp
    BCCONV    = 0e+0_fp

    !=================================================================
    ! For species with dry deposition, the loss rate of dry dep is 
    ! combined in chem loss term.
    !
    ! Conversion from hydrophobic to hydrophilic:  
    ! e-folding time 1.15 days 
    ! ----------------------------------------
    ! Use an e-folding time of 1.15 days or a convertion rate 
    ! of 1.0e-5 /sec. 
    !
    ! Hydrophobic(2) --> Hydrophilic(1) ,  k  = 1.0e-5          
    ! Both aerosols are dry-deposited,     kd = Dvel/DELZ (sec-1)      
    !=================================================================
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initial BC mass [kg]
       TC0  = TC(I,J,L)

       ! Zero drydep freq
       ! ### NOTE: Remove this later, but need to make
       ! ### sure we don't incur numerical diffs (bmy, 6/12/15)
       FREQ = 0e+0_fp

       ! Amount of PAHBCPO left after chemistry and drydep [kg]
       RKT  = ( KBC + FREQ ) * DTCHEM
       CNEW = TC0 * EXP( -RKT )

       ! Prevent underflow condition
       IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

       ! Amount of PAHBCPO converted to PAHBCPI [kg/timestep]
       BCCONV(I,J,L) = ( TC0 - CNEW ) * KBC / ( KBC + FREQ )
	   
       ! Store new concentration back into species array
       TC(I,J,L) = CNEW
    ENDDO
    ENDDO
    ENDDO 
  END SUBROUTINE CHEM_PAHBCPO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CHEM_PAHBCPI
!
! !DESCRIPTION: Subroutine CHEM\_PAHBCPI calculates dry deposition of 
!  hydrophilic PAHBC.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_PAHBCPI( am_I_Root, Input_Opt, TC, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root             ! Root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
!
! !INPUT/OUTPUT PARAMETERS: 
!
    REAL(fp),       INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! H-philic BC [kg]
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: I,      J,      L
    REAL(fp) :: L_FRAC, TC0, CNEW, CCV

    !=================================================================
    ! CHEM_BCPI begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
	
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initial H-philic BC [kg]
       TC0 = TC(I,J,L)

       ! H-philic PAHBC that used to be H-phobic PAHBC [kg]
       CCV = BCCONV(I,J,L)
         
       ! Add the amount of converted PAHBCPO to PAHBCPI
       CNEW = TC0 + CCV
      
       ! Prevent underflow condition
       IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

       ! Save new concentration of H-philic IC in species array
       TC(I,J,L) = CNEW

    ENDDO
    ENDDO
    ENDDO

    !=================================================================
    ! Zero out the BCCONV array for the next iteration
    !=================================================================
    BCCONV = 0e+0_fp

  END SUBROUTINE CHEM_PAHBCPI
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CHEM_PAHOCPO
!
! !DESCRIPTION: Subroutine CHEM_BCPO converts hydrophobic PAHOC to hydrophilic 
!  PAHBC and calculates the dry deposition of hydrophobic PAHOC.
!  COPIED and EDITED chem_bcpo 
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_PAHOCPO( am_I_Root, Input_Opt, TC, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root             ! Root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
!
! !INPUT/OUTPUT PARAMETERS:
! 
    REAL(fp),       INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! H-phobic BC [kg]
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: I,      J,   L
    REAL(fp)              :: DTCHEM, KOC, FREQ, TC0, CNEW, RKT
	
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: OC_LIFE = 1.15e+0_fp

    !=================================================================
    ! CHEM_OCPO begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Initialize
    KOC       = 1.e+0_fp / ( 86400e+0_fp * OC_LIFE )
    DTCHEM    = GET_TS_CHEM() * 60e+0_fp
    OCCONV    = 0e+0_fp

    !=================================================================
    ! For species with dry deposition, the loss rate of dry dep is 
    ! combined in chem loss term.
    !
    ! Conversion from hydrophobic to hydrophilic:  
    ! e-folding time 1.15 days 
    ! ----------------------------------------
    ! Use an e-folding time of 1.15 days or a convertion rate 
    ! of 1.0e-5 /sec. 
    !
    ! Hydrophobic(2) --> Hydrophilic(1) ,  k  = 1.0e-5          
    ! Both aerosols are dry-deposited,     kd = Dvel/DELZ (sec-1)      
    !=================================================================
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initial BC mass [kg]
       TC0  = TC(I,J,L)

       ! Zero drydep freq
       ! ### NOTE: Remove this later, but need to make
       ! ### sure we don't incur numerical diffs (bmy, 6/12/15)
       FREQ = 0e+0_fp

       ! Amount of PAHBCPO left after chemistry and drydep [kg]
       RKT  = ( KOC + FREQ ) * DTCHEM
       CNEW = TC0 * EXP( -RKT )

       ! Prevent underflow condition
       IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

       ! Amount of PAHBCPO converted to PAHBCPI [kg/timestep]
       OCCONV(I,J,L) = ( TC0 - CNEW ) * KOC / ( KOC + FREQ )

       ! Store new concentration back into species array
       TC(I,J,L) = CNEW
    ENDDO
    ENDDO
    ENDDO 
  END SUBROUTINE CHEM_PAHOCPO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CHEM_PAHOCPI
!
! !DESCRIPTION: Subroutine CHEM\_PAHBCPI calculates dry deposition of 
!  hydrophilic PAHOC.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_PAHOCPI( am_I_Root, Input_Opt, TC, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root             ! Root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt             ! Input Options
!
! !INPUT/OUTPUT PARAMETERS: 
!
    REAL(fp),       INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR) ! H-philic BC [kg]
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC                    ! Success?
!
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: I,      J,      L
    REAL(fp) :: L_FRAC, TC0, CNEW, CCV

    !=================================================================
    ! CHEM_BCPI begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
	
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Initial H-philic BC [kg]
       TC0 = TC(I,J,L)

       ! H-philic PAHBC that used to be H-phobic PAHBC [kg]
       CCV = OCCONV(I,J,L)
         
       ! Add the amount of converted PAHBCPO to PAHBCPI
       CNEW = TC0 + CCV
      
       ! Prevent underflow condition
       IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

       ! Save new concentration of H-philic IC in species array
       TC(I,J,L) = CNEW

    ENDDO
    ENDDO
    ENDDO

    !=================================================================
    ! Zero out the BCCONV array for the next iteration
    !=================================================================
    OCCONV = 0e+0_fp

  END SUBROUTINE CHEM_PAHOCPI
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PARITITIONING_PAH
!
! !DESCRIPTION: Subroutine CHEM\_BCPO converts hydrophobic PAHBC to hydrophilic 
!  PAHBC and calculates the dry deposition of hydrophobic PAHBC.
!  COPIED and EDITED chem_bcpo 
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PARTITIONING_PAH( am_I_Root, Input_Opt, PAH_NAME, TC1, TC2, TC3,&
							   TC4, TC5, TC6, TC7, TC8, TC9, AIRVOL, T, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN) :: am_I_Root                 ! Root CPU?
    TYPE(OptInput),   INTENT(IN) :: Input_Opt                 ! Input Options
	CHARACTER(LEN=4), INTENT(IN) :: PAH_NAME                  ! Tracer name of PAH
	REAL(fp),         INTENT(IN) :: TC6(IIPAR,JJPAR,LLPAR)    ! BCPO
	REAL(fp),         INTENT(IN) :: TC7(IIPAR,JJPAR,LLPAR)    ! BCPI
	REAL(fp),         INTENT(IN) :: TC8(IIPAR,JJPAR,LLPAR)    ! OCPO
	REAL(fp),         INTENT(IN) :: TC9(IIPAR,JJPAR,LLPAR)    ! OCPI
	REAL(fp),         INTENT(IN) :: AIRVOL(IIPAR,JJPAR,LLPAR) ! Air volume [m3]
	REAL(fp),         INTENT(IN) :: T(IIPAR,JJPAR,LLPAR)      ! Temperature [K]
!!
! !INPUT/OUTPUT PARAMETERS: 
!
	REAL(fp),       INTENT(INOUT) :: TC1(IIPAR,JJPAR,LLPAR) ! PAH
	REAL(fp),       INTENT(INOUT) :: TC2(IIPAR,JJPAR,LLPAR) ! PAHBCPO
	REAL(fp),       INTENT(INOUT) :: TC3(IIPAR,JJPAR,LLPAR) ! PAHBCPI
	REAL(fp),       INTENT(INOUT) :: TC4(IIPAR,JJPAR,LLPAR) ! PAHOCPO
	REAL(fp),       INTENT(INOUT) :: TC5(IIPAR,JJPAR,LLPAR) ! PAHOCPI
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC       ! Success?
!
!-----------------------------------------------------------------------------
  !BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: I,        J,      L
	REAL(fp) :: TC6_C, TC7_C, TC8_C, TC9_C
	REAL(fp) :: TC_TOTAL_BC, TC_TOTAL_OC, RATIO_BC, RATIO_OC, TC_TOTAL_PAH
	REAL(fp) :: K_SA, K_OA, K_SA_T, K_OA_T, H_SOLV, R
	REAL(fp) :: SIGMA, KP, PAH_ON_PM, PAH_ON_BC, PAH_ON_OC, KP_RATIO
	
    !=================================================================
    ! PARTITIONING_PAH begins here!
    !=================================================================
    ! Assume success
    RC            = GC_SUCCESS
	
	! Initialize
	R             = 8.3144598e-3_fp       !kJ mol k-1

	! ACE
	IF (TRIM(PAH_NAME) .eq. 'ACE') THEN
		H_SOLV    = -2.5e+1_fp           ! Set to the same as FLO
		K_SA      = 7.41e+7_fp			 ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
		K_OA      = 2.19e+6_fp			 !https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! ACY
	IF (TRIM(PAH_NAME) .eq. 'ACY') THEN
		H_SOLV    = -2.5e+1_fp           ! Set to the same as FLO
		K_SA      = 7.41e+7_fp			 ! https://doi.org/10.1016/j.envint.2016.12.002 JMK 
		K_OA      = 3.31e+11_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! ANT
	IF (TRIM(PAH_NAME) .eq. 'ANT') THEN
		H_SOLV    = -2.9e+1_fp           ! http://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 3.16e+9_fp			 ! https://pubs.acs.org/doi/pdf/10.1021/es035337q
		K_OA      = 5.13e+7_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! BAA
	IF (TRIM(PAH_NAME) .eq. 'BAA') THEN
		H_SOLV    = -3.5e+1_fp           ! http://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 1.62e+11_fp			 ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
		K_OA      = 1.90e+10_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! BAP
	IF (TRIM(PAH_NAME) .eq. 'BAP') THEN
		H_SOLV    = -2.2e+1_fp           ! http://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 1.00e+13_fp			 ! https://pubs.acs.org/doi/pdf/10.1021/es035337q
		K_OA      = 3.63e+11_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
        ! BAPB - treat as PYR
        IF (TRIM(PAH_NAME) .eq. 'BAPB') THEN
                H_SOLV    = -3.7e+1_fp           ! https://pubs.acs.org/doi/pdf/10.1021/es051083s same as BKF
                K_SA      = 2.09e+12_fp                  ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
                K_OA      = 2.34e+11_fp                  ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
        ENDIF        ! BAPC - treat as PYR
        IF (TRIM(PAH_NAME) .eq. 'BAPC') THEN
                H_SOLV    = -3.7e+1_fp           ! https://pubs.acs.org/doi/pdf/10.1021/es051083s BKF
                K_SA      = 2.09e+12_fp                  ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
                K_OA      = 2.34e+11_fp                  ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
        ENDIF
	! BBF
	IF (TRIM(PAH_NAME) .eq. 'BBF') THEN
		H_SOLV    = -3.7e+1_fp           ! Set to the the same as BKF
		K_SA      = 2.00e+12_fp			 ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
		K_OA      = 2.19e+11_fp			 !https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! BHP
	IF (TRIM(PAH_NAME) .eq. 'BHP') THEN
		H_SOLV    = -3.3e+1_fp           ! https://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 2.82e+13_fp			 ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
		K_OA      = 3.55e+12_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! BKF
	IF (TRIM(PAH_NAME) .eq. 'BKF') THEN
		H_SOLV    = -3.7e+1_fp           ! https://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 2.09e+12_fp			 ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
		K_OA      = 2.34e+11_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! CHR
	IF (TRIM(PAH_NAME) .eq. 'CHR') THEN
		H_SOLV    = -3.5e+1_fp           ! Set to the same as BAA
		K_SA      = 1.26e+12_fp			 ! https://pubs.acs.org/doi/pdf/10.1021/es035337q
		K_OA      = 2.00e+10_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! DAHA
	IF (TRIM(PAH_NAME) .eq. 'DAHA') THEN
		H_SOLV    = -3.3e+1_fp           ! Set to the same as BHP
		K_SA      = 2.82e+13_fp			 ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
		K_OA      = 3.63e+12_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! FLA
	IF (TRIM(PAH_NAME) .eq. 'FLA') THEN
		H_SOLV    = -3.1e+1_fp           ! http://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 3.16e+10_fp			 ! https://pubs.acs.org/doi/pdf/10.1021/es035337q
		K_OA      = 5.75e+8_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! FLO
	IF (TRIM(PAH_NAME) .eq. 'FLO') THEN
		H_SOLV    = -2.5e+1_fp           ! http://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 3.98e+8_fp			 ! https://pubs.acs.org/doi/pdf/10.1021/es035337q
		K_OA      = 7.94e+6_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! ICDP
	IF (TRIM(PAH_NAME) .eq. 'ICDP') THEN
		H_SOLV    = -2.2e+1_fp           ! Set to the same as BAP
		K_SA      = 2.82e+13_fp			 ! https://doi.org/10.1016/j.envint.2016.12.002 JMK
		K_OA      = 2.69e+12_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! NAP
	IF (TRIM(PAH_NAME) .eq. 'NAP') THEN
		H_SOLV    = -2.2e+1_fp           ! http://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 1.00e+10_fp			 ! Set to an arbitrary v-low value until literature value is added 
		K_OA      = 1.00e+7_fp			 ! Set to an arbitrary v-low value until literature value is added
	ENDIF
	! PHE
	IF (TRIM(PAH_NAME) .eq. 'PHE') THEN
		H_SOLV    = -2.9e+1_fp           ! http://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 2.51e+9_fp			 ! https://pubs.acs.org/doi/pdf/10.1021/es035337q
		K_OA      = 4.79e+7_fp			 ! https://www.sciencedirect.com/science/article/pii/S1352231006005401
	ENDIF
	! PYR
	IF (TRIM(PAH_NAME) .eq. 'PYR') THEN
		H_SOLV    = -2.9e+1_fp           ! http://pubs.acs.org/doi/pdf/10.1021/es051083s
		K_SA      = 3.98e+10_fp			 ! https://pubs.acs.org/doi/pdf/10.1021/es035337q
		K_OA      = 2.00e+7_fp			 ! https://pubs.acs.org/doi/pdf/10.1021/es035337q (wet)
	ENDIF
	
	DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
	  
	  ! Convert kg to ug/m3
	  TC6_C  = (TC6(I,J,L) / AIRVOL(I,J,L))*1e+9_fp ! BCPO
	  TC7_C  = (TC7(I,J,L) / AIRVOL(I,J,L))*1e+9_fp ! BCPI
	  TC8_C  = (TC8(I,J,L) / AIRVOL(I,J,L))*1e+9_fp ! OCPO
	  TC9_C  = (TC9(I,J,L) / AIRVOL(I,J,L))*1e+9_fp ! OCPI
	
	  ! Total BC and OC
	  TC_TOTAL_BC    = TC6_C + TC7_C
	  TC_TOTAL_OC    = TC8_C + TC9_C
	  
	  ! Ratio of PO to PI
	  !RATIO_BC       = TC6_C / TC7_C
	  !RATIO_OC       = TC8_C / TC9_C
	  !JMK, 29/05/2019 - The above ratios were written by Peter Ivatt, but I think they should be inverted. 
          RATIO_BC       = TC7_C / TC6_C
          RATIO_OC       = TC9_C / TC8_C

	  ! Total mass of PAH [kg]
	  TC_TOTAL_PAH   = TC1(I,J,L) + TC2(I,J,L) + TC3(I,J,L) + TC4(I,J,L) + TC5(I,J,L)
	  
	  ! Set the ammount of PAH on BC that cannot parition[kg]
	  ! TC_TOTAL_UNREACT = (TC2(I,J,L) + TC3(I,J,L)) * 0.3_fp
	  
	  ! Remove the unreactable PAH from total PAH
	  ! TC_TOTAL_PAH = TC_TOTAL_PAH - TC_TOTAL_UNREACT
		
	  ! Temperature dependant Soot-Air and Octanal-Air partitioning coeffecients
	  K_SA_T = K_SA*EXP((-H_SOLV/R)*((1.0_fp/T(I,J,L))-(1.0_fp/298.15_fp)))
	  K_OA_T = K_OA*EXP((-H_SOLV/R)*((1.0_fp/T(I,J,L))-(1.0_fp/298.15_fp)))
	  
	  ! Partitioning coefficent (1e-12 converts the rates from L kg-1 to ug m3)
	  KP = 1e-12_fp*((K_SA_T * TC_TOTAL_BC) + (K_OA_T * TC_TOTAL_OC))
	  
	  ! PAH on BC vs Gas phase
	  SIGMA    = (1_fp+(1_fp/(KP)))**-1_fp

	  ! Paritioning onto the gas phase
	  TC1(I,J,L) = TC_TOTAL_PAH * (1_fp-SIGMA)
	  
	  ! Kp RATIO of BC to OC
	  KP_RATIO = (K_SA_T * TC_TOTAL_BC) / (K_OA_T * TC_TOTAL_OC)
	  
	  ! Paritioning onto the BC phase
	  PAH_ON_PM  = (TC_TOTAL_PAH * SIGMA)
	  
	  PAH_ON_BC = PAH_ON_PM * (1.0_fp-((1.0_fp+KP_RATIO)**-1.0_fp))
	  PAH_ON_OC = PAH_ON_PM * ((1.0_fp+KP_RATIO)**-1.0_fp)
	  
	  ! Add PAH that was unavaible for partitiong back on
	  ! PAH_ON_BC = PAH_ON_BC + TC_TOTAL_UNREACT
	  
	  ! Paritioning onto the BCPO phase
	  TC2(I,J,L) = PAH_ON_BC * ((1.0_fp+RATIO_BC)**-1.0_fp)
	  
	  ! Paritioning onto the BCPI phase
	  TC3(I,J,L) = PAH_ON_BC * (1.0_fp-((1.0_fp+RATIO_BC)**-1.0_fp))
	  
	  ! Paritioning onto the BCPO phase
	  TC4(I,J,L) = PAH_ON_OC * ((1.0_fp+RATIO_OC)**-1.0_fp)
	  
	  ! Paritioning onto the BCPI phase
	  TC5(I,J,L) = PAH_ON_OC * (1.0_fp-((1.0_fp+RATIO_OC)**-1.0_fp))
	  
    ENDDO
    ENDDO
    ENDDO
   
  END SUBROUTINE PARTITIONING_PAH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PPFER_PARTITIONING_PAH
!
! !DESCRIPTION: Subroutine CHEM\_BCPO converts hydrophobic PAHBC to hydrophilic 
!  PAHBC and calculates the dry deposition of hydrophobic PAHBC.
!  COPIED and EDITED chem_bcpo 
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PPFER_PARTITIONING_PAH( am_I_Root, Input_Opt, PAH_NAME, TC1, TC2, TC3, TC4, &
						             TC5, TC6, TC7, TC8, TC9, AIRVOL, T, RC )
!
! !USES:
!
    USE CMN_DIAG_MOD
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN) :: am_I_Root                 ! Root CPU?
    TYPE(OptInput),   INTENT(IN) :: Input_Opt                 ! Input Options
	CHARACTER(LEN=4), INTENT(IN) :: PAH_NAME                  ! Tracer name of PAH
	REAL(fp),         INTENT(IN) :: TC6(IIPAR,JJPAR,LLPAR)    ! BCPO
	REAL(fp),         INTENT(IN) :: TC7(IIPAR,JJPAR,LLPAR)    ! BCPI
	REAL(fp),         INTENT(IN) :: TC8(IIPAR,JJPAR,LLPAR)    ! OCPO
	REAL(fp),         INTENT(IN) :: TC9(IIPAR,JJPAR,LLPAR)    ! OCPI
	REAL(fp),         INTENT(IN) :: AIRVOL(IIPAR,JJPAR,LLPAR) ! Air volume [m3]
	REAL(fp),         INTENT(IN) :: T(IIPAR,JJPAR,LLPAR)      ! Temperature [K]
!!
! !INPUT/OUTPUT PARAMETERS: 
!
	REAL(fp),       INTENT(INOUT) :: TC1(IIPAR,JJPAR,LLPAR) ! PAH
	REAL(fp),       INTENT(INOUT) :: TC2(IIPAR,JJPAR,LLPAR) ! PAHBCPO
	REAL(fp),       INTENT(INOUT) :: TC3(IIPAR,JJPAR,LLPAR) ! PAHBCPI
	REAL(fp),       INTENT(INOUT) :: TC4(IIPAR,JJPAR,LLPAR) ! PAHOCPO
	REAL(fp),       INTENT(INOUT) :: TC5(IIPAR,JJPAR,LLPAR) ! PAHOCPI
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC       ! Success?
!
!-----------------------------------------------------------------------------
  !BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: I,        J,      L
	REAL(fp) :: TC6_C, TC7_C, TC8_C, TC9_C
	REAL(fp) :: TC_TOTAL_BC, TC_TOTAL_OC, RATIO_BC, RATIO_OC, TC_TOTAL_PAH
	REAL(fp) :: PARA_E, PARA_S, PARA_A, PARA_B, PARA_V, PARA_L
	REAL(fp) :: KP_OC, KP_BC, KP_WSOC, WSOC_RATIO, DEL_H, R
	REAL(fp) :: SIGMA, KP, PAH_ON_PM, PAH_ON_BC, PAH_ON_OC, KP_RATIO
	
    !=================================================================
    ! PARTITIONING_PAH begins here!
    !=================================================================
    ! Assume success
    RC            = GC_SUCCESS
	
	! Copy fields from INPUT_OPT to local variables for use below
    WSOC_RATIO    = Input_Opt%WSOC_RATIO
	
    ! Initialize
	R             = 8.3144598e-3_fp       !kJ mol k-1
	
	! PAH PARAMmeters https://www.sciencedirect.com/science/article/pii/S0021967314012370
	! ACE
	IF (TRIM(PAH_NAME) .eq. 'ACE') THEN
       PARA_E         = 1.453_fp  ! Excess molar refraction
	   PARA_S		  = 0.951_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.221_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.2586_fp ! McGowan molar volume
	   PARA_L		  = 6.709_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! ACY
	IF (TRIM(PAH_NAME) .eq. 'ACY') THEN
       PARA_E         = 1.553_fp  ! Excess molar refraction
	   PARA_S		  = 1.125_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.214_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.2156_fp ! McGowan molar volume
	   PARA_L		  = 6.395_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! ANT
	IF (TRIM(PAH_NAME) .eq. 'ANT') THEN
       PARA_E         = 1.981_fp  ! Excess molar refraction
	   PARA_S		  = 1.284_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.269_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.4544_fp ! McGowan molar volume
	   PARA_L		  = 7.735_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! BAA
	IF (TRIM(PAH_NAME) .eq. 'BAA') THEN
       PARA_E         = 2.735_fp  ! Excess molar refraction
	   PARA_S		  = 1.678_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.368_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.8234_fp ! McGowan molar volume
	   PARA_L		  = 10.124_fp ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! BAP
	IF (TRIM(PAH_NAME) .eq. 'BAP') THEN
       PARA_E         = 3.023_fp  ! Excess molar refraction
	   PARA_S		  = 1.846_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.418_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.9536_fp ! McGowan molar volume
	   PARA_L		  = 11.54_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
        ! BAPB - 1-nitropyrene (Ariyasena et al., 2014)
        IF (TRIM(PAH_NAME) .eq. 'BAPB') THEN
       PARA_E         = 2.809_fp  ! Excess molar refraction
           PARA_S                 = 2.074_fp  ! Dipolarity/polarizabitliy
           PARA_A                 = 0.0_fp    ! Solute H-bond acidity
           PARA_B                 = 0.325_fp  ! Solute H-bond basicity
           PARA_V                 = 1.7588_fp ! McGowan molar volume
           PARA_L                 = 10.456_fp  ! Log of solute hexadecane-air partitioning coefficent
        ENDIF
        ! BAPC - treat at 1-nitropyrene, despite this representing dinitropyrene
        IF (TRIM(PAH_NAME) .eq. 'BAPC') THEN
       PARA_E         = 2.809_fp  ! Excess molar refraction
           PARA_S                 = 2.074_fp  ! Dipolarity/polarizabitliy
           PARA_A                 = 0.0_fp    ! Solute H-bond acidity
           PARA_B                 = 0.325_fp  ! Solute H-bond basicity
           PARA_V                 = 1.7588_fp ! McGowan molar volume
           PARA_L                 = 10.456_fp  ! Log of solute hexadecane-air partitioning coefficent        
        ENDIF
	! BBF - JMK Abraham et al. (2007) - https://doi.org/10.1016/j.fluid.2006.11.007 
	IF (TRIM(PAH_NAME) .eq. 'BBF') THEN
       PARA_E         = 3.194_fp  ! Excess molar refraction
	   PARA_S		  = 1.820_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.400_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.9536_fp ! McGowan molar volume
	   PARA_L		  = 12.70_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! BHP
	IF (TRIM(PAH_NAME) .eq. 'BHP') THEN
       PARA_E         = 3.612_fp  ! Excess molar refraction
	   PARA_S		  = 2.110_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.436_fp  ! Solute H-bond basicity
	   PARA_V		  = 2.0838_fp ! McGowan molar volume
	   PARA_L		  = 12.707_fp ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! BKF - JMK Abraham et al. (2007) - https://doi.org/10.1016/j.fluid.2006.11.007
	IF (TRIM(PAH_NAME) .eq. 'BKF') THEN
       PARA_E         = 3.190_fp  ! Excess molar refraction
	   PARA_S		  = 1.91_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.33_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.9536_fp ! McGowan molar volume
	   PARA_L		  = 11.607_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! CHR
	IF (TRIM(PAH_NAME) .eq. 'CHR') THEN
       PARA_E         = 2.593_fp  ! Excess molar refraction
	   PARA_S		  = 1.660_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.294_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.8234_fp ! McGowan molar volume
	   PARA_L		  = 10.142_fp ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! DAHA
	IF (TRIM(PAH_NAME) .eq. 'DAHA') THEN
       PARA_E         = 3.827_fp  ! Excess molar refraction
	   PARA_S		  = 2.261_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.549_fp  ! Solute H-bond basicity
	   PARA_V		  = 2.1924_fp ! McGowan molar volume
	   PARA_L		  = 12.552_fp ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! FLA
	IF (TRIM(PAH_NAME) .eq. 'FLA') THEN
       PARA_E         = 2.348_fp  ! Excess molar refraction
	   PARA_S		  = 1.479_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.300_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.5846_fp ! McGowan molar volume
	   PARA_L		  = 8.733_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! FLO
	IF (TRIM(PAH_NAME) .eq. 'FLO') THEN
       PARA_E         = 1.660_fp  ! Excess molar refraction
	   PARA_S		  = 1.104_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.256_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.3565_fp ! McGowan molar volume
	   PARA_L		  = 6.948_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! ICDP - JMK Abraham et al. (2007) - https://doi.org/10.1016/j.fluid.2006.11.007
	IF (TRIM(PAH_NAME) .eq. 'ICDP') THEN
       PARA_E         = 3.610_fp  ! Excess molar refraction
	   PARA_S		  = 1.93_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.42_fp  ! Solute H-bond basicity
	   PARA_V		  = 2.0838_fp ! McGowan molar volume
	   PARA_L		  = 12.699_fp ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! NAP
	IF (TRIM(PAH_NAME) .eq. 'NAP') THEN
       PARA_E         = 1.230_fp  ! Excess molar refraction
	   PARA_S		  = 0.906_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.191_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.0854_fp ! McGowan molar volume
	   PARA_L		  = 5.157_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! PHE
	IF (TRIM(PAH_NAME) .eq. 'PHE') THEN
       PARA_E         = 1.917_fp  ! Excess molar refraction
	   PARA_S		  = 1.275_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.285_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.4544_fp ! McGowan molar volume
	   PARA_L		  = 7.712_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	! PYR
	IF (TRIM(PAH_NAME) .eq. 'PYR') THEN
       PARA_E         = 2.241_fp  ! Excess molar refraction
	   PARA_S		  = 1.475_fp  ! Dipolarity/polarizabitliy
	   PARA_A		  = 0.0_fp    ! Solute H-bond acidity
	   PARA_B		  = 0.283_fp  ! Solute H-bond basicity
	   PARA_V		  = 1.5846_fp ! McGowan molar volume
	   PARA_L		  = 8.974_fp  ! Log of solute hexadecane-air partitioning coefficent
	ENDIF
	
	DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
	
	  !convert kg to ug/m3
	  TC6_C  = (TC6(I,J,L) / AIRVOL(I,J,L))*1e+9_fp ! BCPO
	  TC7_C  = (TC7(I,J,L) / AIRVOL(I,J,L))*1e+9_fp ! BCPI
	  TC8_C  = (TC8(I,J,L) / AIRVOL(I,J,L))*1e+9_fp ! OCPO
	  TC9_C  = (TC9(I,J,L) / AIRVOL(I,J,L))*1e+9_fp ! OCPI
	
	  !Total BC and OC
	  TC_TOTAL_BC    = TC6_C + TC7_C
	  TC_TOTAL_OC    = TC8_C + TC9_C
	  
	  !Ratio of PO to PI
	  !RATIO_BC       = TC6_C / TC7_C
	  !RATIO_OC       = TC8_C / TC9_C
	  !JMK, 29/05/2019 - The above ratios were written by Peter Ivatt, but I think they should be inverted. 
          RATIO_BC       = TC7_C / TC6_C
          RATIO_OC       = TC9_C / TC8_C

	  !Total mass of PAH [kg]
	  TC_TOTAL_PAH   = TC1(I,J,L) + TC2(I,J,L) + TC3(I,J,L) + TC4(I,J,L) + TC5(I,J,L)
	  
	  !System PARAMmeters from suplimenatry info (https://pubs.acs.org/doi/suppl/10.1021/acs.est.6b02158)
	  !Kp BCPO & BCPI (m3(air) m-2(surface) units)
	  KP_BC = 10.0_fp**((PARA_A*2.7_fp)+(PARA_B*2.45_fp)+(PARA_L*1.09_fp)-8.47_fp)
	 ! Where 1.821E-5 is the surface area (m2 ug-1) of soot Jonker and Koelmans 2002 (Avg, traffic, wood,coal diesel) 
	  KP_BC = KP_BC  * 1.821e-5_fp 
	  
	  ! Kp (L(air)L-1(solvent) units)
	  ! OCTANOL surogate for the insoluable componenet of OC
	  KP_OC = 10.0_fp**((PARA_E*(-0.21_fp))+(PARA_S*0.56_fp)+(PARA_A*3.51_fp)+(PARA_B*0.75_fp)+(PARA_L*0.94_fp)-0.22_fp)
	  ! Convert units to (m3(air)L-1(solvent)
	  KP_OC = KP_OC / 1e+3_fp
	  ! Where 8.2E8 is the density (ug L-3) of octanol
	  KP_OC = (KP_OC / 8.2e+8_fp) 
	  
	  ! Dimethyl sulfoxide surogate for water soluable OC
	  KP_WSOC = 10.0_fp**((PARA_E*(-0.22_fp))+(PARA_S*2.90_fp)+(PARA_A*5.04_fp)+(PARA_L*0.72_fp)-0.56_fp)
	  ! Convert units to (m3(air)L-1(solvent)
	  KP_WSOC = KP_WSOC / 1e+3_fp
	  ! Where 8.2E8 is the density (ug L-3) of dimethyl sulfoxide
	  KP_WSOC = (KP_WSOC / 1.1e+9_fp) 
	  
      ! Temperature dependance from sorption enthaply derived from ppFELR PARAMmeters
	  ! https://pubs.acs.org/doi/pdf/10.1021/ac070265x
	  DEL_H   = (-12.8_fp*PARA_V)+(-4.3_fp*PARA_L)+(-46.6*PARA_A)+(-17.6*PARA_S)+2.7
	  KP_BC   = KP_BC*EXP((-DEL_H/R)*((1.0_fp/T(I,J,L))-(1.0_fp/298_fp)))
	  KP_OC   = KP_OC*EXP((-DEL_H/R)*((1.0_fp/T(I,J,L))-(1.0_fp/298_fp)))
      KP_WSOC = KP_WSOC*EXP((-DEL_H/R)*((1.0_fp/T(I,J,L))-(1.0_fp/298_fp)))
	  
	  ! Use the WSOC rate for the WSOC ratio of OC.
	  KP_OC = (KP_WSOC * WSOC_RATIO) + (KP_OC * (1-WSOC_RATIO))
	  
	  ! Total PAH in particulate phase is equal to the sum of all the KP multiplied by the tracer concentration. 
	  KP = (KP_BC * TC_TOTAL_BC) + (KP_OC * TC_TOTAL_OC)

	  ! PAH on particulate vs Gas phase
	  SIGMA    = (1.0_fp+(1.0_fp/KP))**-1.0_fp
	   
	  ! Paritioning onto the gas phase
	  TC1(I,J,L) = TC_TOTAL_PAH * (1_fp-SIGMA)
	  
	  ! Kp RATIO of BC to OC
	  KP_RATIO = (KP_BC * TC_TOTAL_BC) / (KP_OC * TC8_C)
	  
	  ! Paritioning onto the BC phase
	  PAH_ON_PM  = (TC_TOTAL_PAH * SIGMA)
	  
	  PAH_ON_BC = PAH_ON_PM * (1.0_fp-((1.0_fp+KP_RATIO)**-1.0_fp))
	  PAH_ON_OC = PAH_ON_PM * ((1.0_fp+KP_RATIO)**-1.0_fp)
	  
	  ! Paritioning onto the BCPO phase
	  TC2(I,J,L) = PAH_ON_BC * ((1.0_fp+RATIO_BC)**-1.0_fp)
	  
	  ! Paritioning onto the BCPI phase
	  TC3(I,J,L) = PAH_ON_BC * (1.0_fp-((1.0_fp+RATIO_BC)**-1.0_fp))
	   
	  ! Paritioning onto the BCPO phase
	  TC4(I,J,L) = PAH_ON_OC * ((1.0_fp+RATIO_OC)**-1.0_fp)
	  
	  ! Paritioning onto the BCPI phase
	  TC5(I,J,L) = PAH_ON_OC * (1.0_fp-((1.0_fp+RATIO_OC)**-1.0_fp))
	  
	  ! Prevent underflow condition
	  IF (TC1(I,J,L) < SMALLNUM ) TC1(I,J,L) = SMALLNUM
      IF (TC2(I,J,L) < SMALLNUM ) TC2(I,J,L) = SMALLNUM
	  IF (TC3(I,J,L) < SMALLNUM ) TC3(I,J,L) = SMALLNUM
	  IF (TC4(I,J,L) < SMALLNUM ) TC4(I,J,L) = SMALLNUM
	  IF (TC5(I,J,L) < SMALLNUM ) TC5(I,J,L) = SMALLNUM
	  
	  
    ENDDO
    ENDDO
    ENDDO
   
  END SUBROUTINE PPFER_PARTITIONING_PAH
END MODULE PAH_MOD
