      Subroutine Realspace(Ureal)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccc
C     Real Path; Also Direct Calculation   C
C     Loop Over All Particle Pairs         C
Cccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J
      Double Precision Dx,Dy,Dz,Ureal,R,R2,Ir,Dderfc

      Ureal = 0.0d0

C     Start Modification
c    1. For all particle pairs calculate the distance in x, y, and z.
c    2. Apply periodic boundary conditions where necessary.
c    3. Calculate the real-space contribution to the energy. 
c    (Use the "Dderfc(x)" function to calculate the Error-Function Complement.) 

C     End   Modification

      Return
      End
