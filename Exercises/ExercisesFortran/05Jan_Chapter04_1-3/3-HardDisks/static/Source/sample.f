      Subroutine Sample(Ichoise)
      Implicit None

      Include 'system.inc'

Cccccccccccccccccccccccccccccccccccccccccccccccccc
C     Samples The Radial Distribution Function   C
Cccccccccccccccccccccccccccccccccccccccccccccccccc

      Integer I,J,Maxx,Ichoise,A
      Double Precision BoxSize

      Parameter(Maxx = 500)
      Parameter(BoxSize = 10.0)

      Double Precision Ngr,Gg(Maxx),Delta,R2,DX,DY,R,ab,nid

      Save Ngr,Gg,Delta

      If(Ichoise.Eq.1) Then
         Do I=1,Maxx
            Gg(I) = 0.0d0
         Enddo

         Ngr   = 0.0d0
         Delta = 5.0d0/Dble(Maxx-1)
      
      Elseif(Ichoise.Eq.2) Then

Ccccccccccccccccccccccccccccccccccccccccccccccccc
C     Sample The Radial Distribution Function   C
C     Loop Over All Particle Pairs              C
C     See Frenkel/Smit P. 77                    C
C                                               C
C     Delta  = Binsize                          C
C     Ggt/Gg = Counters                         C
C     Maxx   = Maximum Number Of Bins           C
Ccccccccccccccccccccccccccccccccccccccccccccccccc

C     Start Modification

               
C     End   Modification

      Else

Ccccccccccccccccccccccccccccccc
C     Write Results To Disk   C
Ccccccccccccccccccccccccccccccc

         Do I=1,Maxx-1

            r=Delta*(I+0.5d0)
            ab=((i+1)**2-i**2)*Delta**2
            nid= 4.0d0*Datan(1.0d0)*Dble(Npart-1)/(BoxSize*BoxSize)*ab
            Gg(I)=Gg(I)/(Ngr*Dble(Npart)*nid)
         
            Write(21,*) r,Gg(I)
         Enddo
      Endif

      Return
      End
