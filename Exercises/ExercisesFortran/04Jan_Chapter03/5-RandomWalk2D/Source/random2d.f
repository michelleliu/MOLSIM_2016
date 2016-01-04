      Program Random2d
      Implicit None

      Include 'system.inc'


Cccccccccccccccccccccccccc
C     Random Walk 2d     C
Cccccccccccccccccccccccccc

      Integer I,J,K,L,X,Y,Sstmm,Dx,Dy,Xnew,Ynew,NumberOfJumps,
     &     LatticeSize,NumberOfInitializationCylcles,Nfixed,
     &     CycleMultiplication,EMPTY

      Parameter (EMPTY = 0)
      Parameter (CycleMultiplication = 1000)

      Double Precision M1,RandomNumber,Accepted,Attempts,Rm

Ccccccccccccccccccccccccc
C     Initialize Rng    C
Ccccccccccccccccccccccccc

      M1 = 0.001d0*Dble(Mod((10+10*Sstmm()),1000))

      If(M1.Lt.0.001d0) M1 = 0.001d0
      If(M1.Gt.0.999d0) M1 = 0.999d0

      Call InitializeRandomNumberGenerator(M1)
      Call Sample(1)
 
Cccccccccccccccccc
C     Read Input C
Cccccccccccccccccc

      Write(*,*) 'Number of Jumps (x ',CycleMultiplication,') ? '
      Read(*,*) NumberOfJumps

      Write(*,*) 'Lattice Size            ? '
      Read(*,*) LatticeSize

      If(LatticeSize.Lt.10.Or.LatticeSize.Ge.MaxLattice) then
        Write(*,*) 'Error in Lattice Size.....'
        Stop
      End if

      If(NumberOfParticles.gt.(LatticeSize*LatticeSize)) then
         Write(6,*) 'Too many particles for this lattice'
         stop
      End if

      Write(*,*) 'Number of Particles    ? '
      Read(*,*) NumberOfParticles

      If(NumberOfParticles.Lt.1.Or.NumberOfParticles.Ge.MaxParticles) then
        Write(*,*) 'Error in Number of particles....'
        Stop
      End if

      Do I=1,MaxLattice
         Do J=1,MaxLattice
            Lattice(I,J) = 0
         Enddo
      Enddo

      Do I=1,MaxParticles
         ParticlePosX(I) = 0
         ParticlePosY(I) = 0
         Mxx(I) = 0
         Myy(I) = 0
      Enddo

      Accepted   = 0.0d0
      Attempts   = 0.0d0
      NumberOfInitializationCylcles = NumberOfJumps/4

c     question 5: Fix a number of particles. Do not 
c     forget to declare Nfixed
c      Nfixed=int(NumberOfParticles/4)

Cccccccccccccccccccccccccccccccccccccccccccccc
C     Put particles on the lattice at random C
C     Look for an empty lattice site         C
Cccccccccccccccccccccccccccccccccccccccccccccc

      Do J=1,NumberOfParticles
 10      X = 1 + Idint(RandomNumber()*Dble(LatticeSize))
         Y = 1 + Idint(RandomNumber()*Dble(LatticeSize))

         If(Lattice(X,Y).Ne.0) Goto 10

         Lattice(X,Y) = J
         ParticlePosX(J) = X
         ParticlePosY(J) = Y
         Mxx(J) = X
         Myy(J) = Y 
      Enddo
         
Cccccccccccccccccccccccccccccccccc
C     Perform The Random Walk    C
Cccccccccccccccccccccccccccccccccc

      Do K=1,NumberOfJumps
         Do L=1,CycleMultiplication
            
Ccccccccccccccccccccccccccccccccccccccccccccc
C     Choose A Random Site (J)              C
C     Choose A Random Displacement (Rm):    C
C     - Up, Down, Left, Right               C
Ccccccccccccccccccccccccccccccccccccccccccccc

            J  = 1 + Idint(RandomNumber()*Dble(NumberOfParticles))
            Rm = 4.0d0*RandomNumber()

            If(Rm.Lt.1.0d0) Then
               Dx = 1
               Dy = 0
            Elseif(Rm.Lt.2.0d0) Then
               Dx = -1
               Dy = 0
            Elseif(Rm.Lt.3.0d0) Then
               Dx = 0
               Dy = 1
            Else
               Dx = 0
               Dy = -1
            Endif 

Cccccccccccccccccccccccccccccccccccccccccccc
C     New Position                         C
Cccccccccccccccccccccccccccccccccccccccccccc

            Xnew = ParticlePosX(J) + Dx
            Ynew = ParticlePosY(J) + Dy

CcccccccccccccccccccccccccccccccccccccccccccccccccccccC
c     If the new position is outside of the lattice   c
c     boundaries,  use periodic boundary conditions:  c
c     put particle back on the lattice.               c
CcccccccccccccccccccccccccccccccccccccccccccccccccccccC 

            If(Xnew.Lt.1) Then
               Xnew = Xnew + LatticeSize
            Elseif(Xnew.Gt.LatticeSize) Then
               Xnew = Xnew - LatticeSize
            Endif

            If(Ynew.Lt.1) Then
               Ynew = Ynew + LatticeSize
            Elseif(Ynew.Gt.LatticeSize) Then
               Ynew = Ynew - LatticeSize
            Endif

            Attempts = Attempts + 1.0d0

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Check if new lattice site is occupied     C
C     If it is empty, accept the new position   C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

            If(Lattice(Xnew,Ynew).Eq.EMPTY) Then
               Lattice(ParticlePosX(J),ParticlePosY(J)) = EMPTY
               ParticlePosX(J) = Xnew
               ParticlePosY(J) = Ynew
               Lattice(Xnew,Ynew) = J
               Mxx(J) = Mxx(J) + Dx
               Myy(J) = Myy(J) + Dy
               Accepted = Accepted + 1.0d0
            Endif

            If(K.Ge.NumberOfInitializationCylcles) Call Sample(2)
         Enddo
      Enddo

Ccccccccccccccccccccccc
C     Write Results   C
Ccccccccccccccccccccccc

      Call Sample(3)

      Write(6,*)
      Write(6,*) 'Fraction Accepted Jumps : ',Accepted/Attempts
      Write(6,*) 'Lattice Occupation      : ',
     &     Dble(NumberOfParticles)/Dble(LatticeSize*LatticeSize)

      Stop
      End
