Example G: a building with a set-back 

15                      % number of joints 
24			% number of members
1			% number of static load cases

% joint    x         y         z         r
 
 1      -100         0         0         0
 2       100         0         0         0
 3         0        70         0         0         
 4      -100         0        80         0
 5       100         0        80         0
 6         0        70        80         0         
 7      -100         0       180         0
 8       100         0       180         0
 9         0        70       180         0         
10      -100         0       310         0
11       100         0       310         0
12         0        70       310         0         
13      -100         0       510         0
14       100         0       510         0
15         0        70       510         0         

% m    j1 j2   Ax  Asy   Asz   Jxx    Iyy   Izz       E      G p

 1      1  4  100  60    60    500   1000  1000   10000   7000 0
 2      2  5  100  60    60    500   1000  1000   10000   7000 0
 3      3  6  100  60    60    500   1000  1000   10000   7000 0
 4      4  7  100  60    60    500   1000  1000   10000   7000 0
 5      5  8  100  60    60    500   1000  1000   10000   7000 0
 6      6  9  100  60    60    500   1000  1000   10000   7000 0
 7      7 10  100  60    60    500   1000  1000   10000   7000 0
 8      8 11  100  60    60    500   1000  1000   10000   7000 0
 9      9 12  100  60    60    500   1000  1000   10000   7000 0
10     10 13  100  60    60    500   1000  1000   10000   7000 0
11     11 14  100  60    60    500   1000  1000   10000   7000 0
12     12 15  100  60    60    500   1000  1000   10000   7000 0

13      4  5  100  60    60    500   1000  1000   10000   7000 0
14      5  6  100  60    60    500   1000  1000   10000   7000 0
15      6  4  100  60    60    500   1000  1000   10000   7000 0
16      7  8  100  60    60    500   1000  1000   10000   7000 0
17      8  9  100  60    60    500   1000  1000   10000   7000 0
18      9  7  100  60    60    500   1000  1000   10000   7000 0
19     10 11  100  60    60    500   1000  1000   10000   7000 0
20     11 12  100  60    60    500   1000  1000   10000   7000 0
21     12 10  100  60    60    500   1000  1000   10000   7000 0
22     13 14  100  60    60    500   1000  1000   10000   7000 0
23     14 15  100  60    60    500   1000  1000   10000   7000 0
24     15 13  100  60    60    500   1000  1000   10000   7000 0

3                               % number of joints with reactions
% J     x  y  z xx yy zz          1= fixed, 0=free
  1     1  1  1  1  1  1
  2     1  1  1  1  1  1
  3     1  1  1  1  1  1


1                               % 1: include shear deformation
1                               % 1: include geometric stiffness
/tmp/exG-msh                    % mesh data file name
exG.plt				% plot file name
2.0                             % exaggerate mesh deformations
1                               % 1: stiffness analysis, 0: data check only


				% Begin Static Load Case 1 of 1
1                               % number of loaded joints
% J       Fx         Fy     Fz     Mxx     Myy     Mzz
  15       0        200      0       0       0       0

12                               % number of distributed loads
%  j      wx         wy      wz
  13      0          0      -2.361
  14      0          0      -2.361
  15      0          0      -2.361
  16      0          0      -2.361
  17      0          0      -2.361
  18      0          0      -2.361
  19      0          0      -2.361
  20      0          0      -2.361
  21      0          0      -2.361
  22      0          0      -2.361
  23      0          0      -2.361
  24      0          0      -2.361

0                               % number of internal concentrated loads
0                               % number of members with temperature loads
0                               % number of joints with support settlements
				% End   Static Load Case 1 of 1


4                               % number of desired dynamic modes of vibration
1                               % 1: subspace Jacobi     2: Stodola
0                               % 0: consistent mass ... 1: lumped mass matrix
/tmp/exG-m                      % mode shape data file
1e-6                            % mode shape tolerance
0.0                             % shift value ... for unrestrained structures

 1     3e-7  0			% bar numbers, density, and extra mass
 2     3e-7  0
 3     3e-7  0
 4     3e-7  0
 5     3e-7  0
 6     3e-7  0
 7     3e-7  0
 8     3e-7  0
 9     3e-7  0
10     3e-7  0
11     3e-7  0
12     3e-7  0
13     3e-7  0
14     3e-7  0
15     3e-7  0
16     3e-7  0
17     3e-7  0
18     3e-7  0
19     3e-7  0
20     3e-7  0
21     3e-7  0
22     3e-7  0
23     3e-7  0
24     3e-7  0

0                             % number of joints with extra mass or inertia
% j    M        Ixx   Iyy    Izz -  joints and concentrated mass and inertia

4                               % number of modes to animate
 1  2  3  4                     % modes to animate
1                               % pan during animation


________________________________________________________________________________
FRAME3DD version: 20081230               http://frame3dd.sf.net/
GPL Copyright (C) 1992-2008, Henri P. Gavin 
FRAME3DD is distributed in the hope that it will be useful but with no warranty.
For details see the GNU Public Licence: http://www.fsf.org/copyleft/gpl.html
________________________________________________________________________________

Example G: a building with a set-back  
Tue Dec 30 15:25:43 2008
________________________________________________________________________________
15 JOINTS;    24 MEMBERS;    1 LOAD CASES;
3 FIXED JOINTS;   -1073826504 PRESCRIBED DISPLACEMENTS;
For 2D problems, the Y-axis is vertical. 
For 3D problems, the Z-axis is vertical. 
________________________________________________________________________________
J O I N T   D A T A                                         R E S T R A I N T S
  Joint      X              Y              Z         radius  Fx Fy Fz Mx My Mz
    1    -100.000000       0.000000       0.000000    0.000   1  1  1  1  1  1
    2     100.000000       0.000000       0.000000    0.000   1  1  1  1  1  1
    3       0.000000      70.000000       0.000000    0.000   1  1  1  1  1  1
    4    -100.000000       0.000000      80.000000    0.000   0  0  0  0  0  0
    5     100.000000       0.000000      80.000000    0.000   0  0  0  0  0  0
    6       0.000000      70.000000      80.000000    0.000   0  0  0  0  0  0
    7    -100.000000       0.000000     180.000000    0.000   0  0  0  0  0  0
    8     100.000000       0.000000     180.000000    0.000   0  0  0  0  0  0
    9       0.000000      70.000000     180.000000    0.000   0  0  0  0  0  0
   10    -100.000000       0.000000     310.000000    0.000   0  0  0  0  0  0
   11     100.000000       0.000000     310.000000    0.000   0  0  0  0  0  0
   12       0.000000      70.000000     310.000000    0.000   0  0  0  0  0  0
   13    -100.000000       0.000000     510.000000    0.000   0  0  0  0  0  0
   14     100.000000       0.000000     510.000000    0.000   0  0  0  0  0  0
   15       0.000000      70.000000     510.000000    0.000   0  0  0  0  0  0
M E M B E R   D A T A							(local)
  Member J1    J2     Ax   Asy   Asz    Jxx     Iyy     Izz       E       G roll
    1     1     4  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
    2     2     5  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
    3     3     6  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
    4     4     7  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
    5     5     8  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
    6     6     9  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
    7     7    10  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
    8     8    11  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
    9     9    12  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   10    10    13  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   11    11    14  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   12    12    15  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   13     4     5  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   14     5     6  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   15     6     4  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   16     7     8  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   17     8     9  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   18     9     7  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   19    10    11  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   20    11    12  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   21    12    10  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   22    13    14  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   23    14    15  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
   24    15    13  100.0  60.0  60.0  500.0  1000.0  1000.0  10000.0  7000.0   0
  Include shear deformations.
  Include geometric stiffness.

L O A D   C A S E   1   O F   1  ... 

   1 joints  with concentrated loads
  12 members with uniformly distributed loads
   0 members with concentrated point loads
   0 members with temperature loads
   0 joints  with prescribed displacements
 J O I N T   L O A D S  +  E Q U I V A L E N T   J O I N T   L O A D S  (global)
  Joint       Fx          Fy          Fz          Mxx         Myy         Mzz
     4       0.000       0.000    -380.198   -1681.148   10271.640       0.000
     5       0.000       0.000    -380.198   -1681.148  -10271.640       0.000
     6       0.000       0.000    -288.197    3362.296       0.000       0.000
     7       0.000       0.000    -380.198   -1681.148   10271.640       0.000
     8       0.000       0.000    -380.198   -1681.148  -10271.640       0.000
     9       0.000       0.000    -288.197    3362.296       0.000       0.000
    10       0.000       0.000    -380.198   -1681.148   10271.640       0.000
    11       0.000       0.000    -380.198   -1681.148  -10271.640       0.000
    12       0.000       0.000    -288.197    3362.296       0.000       0.000
    13       0.000       0.000    -380.198   -1681.148   10271.640       0.000
    14       0.000       0.000    -380.198   -1681.148  -10271.640       0.000
    15       0.000     200.000    -288.197    3362.296       0.000       0.000
 U N I F O R M   M E M B E R   L O A D S					(local)
  Member      Wx               Wy               Wz
    13       0.00000000       0.00000000      -2.36100000
    14       0.00000000       0.00000000      -2.36100000
    15       0.00000000       0.00000000      -2.36100000
    16       0.00000000       0.00000000      -2.36100000
    17       0.00000000       0.00000000      -2.36100000
    18       0.00000000       0.00000000      -2.36100000
    19       0.00000000       0.00000000      -2.36100000
    20       0.00000000       0.00000000      -2.36100000
    21       0.00000000       0.00000000      -2.36100000
    22       0.00000000       0.00000000      -2.36100000
    23       0.00000000       0.00000000      -2.36100000
    24       0.00000000       0.00000000      -2.36100000

E L A S T I C   S T I F F N E S S   A N A L Y S I S   via  L D L'  decomposition


L O A D   C A S E   1   O F   1  ... 

J O I N T   D I S P L A C E M E N T S					(global)
  Joint    X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     4    0.000518    2.556165   -0.023994   -0.054464    0.000726   -0.000046
     5   -0.000518    2.556165   -0.023994   -0.054464   -0.000726    0.000046
     6    0.0         2.548302   -0.288618   -0.044599    0.0         0.0     
     7    0.002870   10.955855   -0.036915   -0.094791   -0.006141    0.000033
     8   -0.002870   10.955855   -0.036915   -0.094791    0.006141   -0.000033
     9    0.0        10.958399   -0.578196   -0.079161    0.0         0.0     
    10    0.000974   27.045015   -0.048380   -0.120219   -0.012993    0.000108
    11   -0.000974   27.045015   -0.048380   -0.120219    0.012993   -0.000108
    12    0.0        27.057952   -0.828383   -0.097824    0.0         0.0     
    13    0.006767   54.937399   -0.070315   -0.100357   -0.008600    0.000303
    14   -0.006767   54.937399   -0.070315   -0.100357    0.008600   -0.000303
    15    0.0        54.973108   -0.994375   -0.057784    0.0         0.0     
M E M B E R   E N D   F O R C E S					(local)
  Member Joint      Nx          Vy         Vz         Txx        Myy        Mzz
     1      1    299.923c    -74.876     -6.383      2.017    165.916 -10077.469
     1      4   -299.923c     74.876      6.383     -2.017    344.548   3320.750
     2      2    299.923c    -74.876      6.383     -2.017   -165.916 -10077.469
     2      5   -299.923c     74.876     -6.383      2.017   -344.548   3320.750
     3      3   3607.723c    -50.248      0.0        0.0        0.0   -11108.943
     3      6  -3607.723c     50.248      0.0        0.0        0.0    -2104.496
     4      4    129.210c    -98.229     31.791     -2.772   -910.370  -9443.349
     4      7   -129.210c     98.229    -31.791      2.772  -2268.987  -1464.871
     5      5    129.210c    -98.229    -31.791      2.772    910.370  -9443.349
     5      8   -129.210c     98.229     31.791     -2.772   2268.987  -1464.871
     6      6   2895.781c     -3.542      0.0        0.0        0.0   -14976.137
     6      9  -2895.781c      3.542      0.0        0.0        0.0    -9731.857
     7      7     88.195c   -102.330     66.539     -2.004  -3804.406  -9292.657
     7     10    -88.195c    102.330    -66.539      2.004  -4845.540  -5429.243
     8      8     88.195c   -102.330    -66.539      2.004   3804.406  -9292.657
     8     11    -88.195c    102.330     66.539     -2.004   4845.540  -5429.243
     9      9   1924.517c      4.660      0.0        0.0        0.0   -16235.531
     9     12  -1924.517c     -4.660      0.0        0.0        0.0   -14142.489
    10     10    109.673c    -70.975     32.010     -3.412  -3412.905  -7670.240
    10     13   -109.673c     70.975    -32.010      3.412  -2989.645  -9583.763
    11     11    109.673c    -70.975    -32.010      3.412   3412.905  -7670.240
    11     14   -109.673c     70.975     32.010     -3.412   2989.645  -9583.763
    12     12    829.961c    -58.050      0.0        0.0        0.0   -15941.125
    12     15   -829.961c     58.050      0.0        0.0        0.0   -18837.427
    13      4      5.184c      0.0      236.100      0.0    -7797.511     -4.602
    13      5     -5.184c      0.0      236.100      0.0     7797.511      4.602
    14      5     40.420c     -0.213    -65.459    219.785  10362.588     -9.390
    14      6    -40.420c      0.213    353.656   -219.785  15206.510    -16.905
    15      6     40.420c      0.213    353.656   -219.785 -15206.510     16.905
    15      4    -40.420c     -0.213    -65.459    219.785 -10362.588      9.390
    16      7     28.698c      0.0      236.100      0.0    -8478.213      3.278
    16      8    -28.698c      0.0      236.100      0.0     8478.213     -3.278
    17      8      7.309c      0.110   -195.094    468.118  18090.177      4.046
    17      9     -7.309c     -0.110    483.291   -468.118  23309.616      9.463
    18      9      7.309c     -0.110    483.291   -468.118 -23309.616     -9.463
    18      7     -7.309c      0.110   -195.094    468.118 -18090.177     -4.046
    19     10      9.736c      0.0      236.100      0.0    -9165.120     10.717
    19     11     -9.736c      0.0      236.100      0.0     9165.120    -10.717
    20     11    -54.242t      0.303   -257.573    739.699  21786.001      9.309
    20     12     54.242t     -0.303    545.770   -739.699  27286.602     27.044
    21     12    -54.242t     -0.303    545.770   -739.699 -27286.602    -27.044
    21     10     54.242t      0.303   -257.573    739.699 -21786.001     -9.309
    22     13     67.669c      0.0      236.100      0.0    -8710.617     29.569
    22     14    -67.669c      0.0      236.100      0.0     8710.617    -29.569
    23     14   -122.342t      0.986   -126.439   1141.466  15081.155     32.980
    23     15    122.342t     -0.986    414.636  -1141.466  18055.244     83.298
    24     15   -122.342t     -0.986    414.636  -1141.466 -18055.244    -83.298
    24     13    122.342t      0.986   -126.439   1141.466 -15081.155    -32.980
R E A C T I O N S							(global)
  Joint       Fx          Fy          Fz         Mxx         Myy         Mzz
     1       6.383     -74.876     299.923   10077.469     165.916       2.017
     2      -6.383     -74.876     299.923   10077.469    -165.916      -2.017
     3       0.0       -50.248    3607.723   11108.943       0.0         0.0  
R M S   E Q U I L I B R I U M    E R R O R: 2.412e-04

B O D A L   A N A L Y S I S   R E S U L T S
  Total Mass:  9.919573e-02     Structural Mass:  9.919573e-02 
J O I N T   M A S S E S	(diagonal of the mass matrix)			(global)
  Joint X-mass      Y-mass      Z-mass      X-inrta     Y-inrta     Z-inrta
     1 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01
     2 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01
     3 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01 9.35022e+01
     4 5.28132e-03 5.56045e-03 5.39348e-03 6.25794e-01 3.08696e+00 2.82725e+00
     5 5.28132e-03 5.56045e-03 5.39348e-03 6.25794e-01 3.08696e+00 2.82725e+00
     6 4.54882e-03 4.64634e-03 4.52622e-03 7.92387e-01 1.14729e+00 1.05807e+00
     7 5.83673e-03 6.11586e-03 5.89348e-03 1.10922e+00 3.57039e+00 2.82975e+00
     8 5.83673e-03 6.11586e-03 5.89348e-03 1.10922e+00 3.57039e+00 2.82975e+00
     9 5.10423e-03 5.20175e-03 5.02622e-03 1.27582e+00 1.63071e+00 1.06057e+00
    10 6.94921e-03 7.22834e-03 6.89348e-03 3.11322e+00 5.57439e+00 2.83475e+00
    11 6.94921e-03 7.22834e-03 6.89348e-03 3.11322e+00 5.57439e+00 2.83475e+00
    12 6.21672e-03 6.31423e-03 6.02622e-03 3.27982e+00 3.63471e+00 1.06557e+00
    13 5.49787e-03 5.77700e-03 5.59348e-03 2.48031e+00 4.94147e+00 2.82825e+00
    14 5.49787e-03 5.77700e-03 5.59348e-03 2.48031e+00 4.94147e+00 2.82825e+00
    15 4.76538e-03 4.86289e-03 4.72622e-03 2.64690e+00 3.00180e+00 1.05907e+00
  Use consistent mass matrix.
N A T U R A L   F R E Q U E N C I E S   & 
M A S S   N O R M A L I Z E D   M O D E   S H A P E S 
 convergence tolerance: 1.000e-06 
  MODE     1:   f= 1.624893 Hz,  T= 0.615425 sec
		X- modal participation factor =   2.9600e-14 
		Y- modal participation factor =   2.5141e-01 
		Z- modal participation factor =   1.7968e-03 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     2   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     3   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     4   2.303e-04   3.286e-01   9.617e-03  -6.797e-03  -8.916e-04   7.893e-06
     5  -2.303e-04   3.286e-01   9.617e-03  -6.797e-03   8.916e-04  -7.893e-06
     6   4.407e-12   3.278e-01  -1.928e-02  -6.101e-03   3.624e-14  -5.812e-14
     7   6.423e-06   1.350e+00   1.919e-02  -1.094e-02  -1.648e-03   4.710e-05
     8  -6.423e-06   1.350e+00   1.919e-02  -1.094e-02   1.648e-03  -4.710e-05
     9   6.649e-12   1.350e+00  -3.846e-02  -9.614e-03  -2.465e-14  -8.284e-14
    10  -5.933e-05   3.087e+00   2.676e-02  -1.192e-02  -2.184e-03   1.086e-04
    11   5.933e-05   3.087e+00   2.676e-02  -1.192e-02   2.184e-03  -1.086e-04
    12  -5.938e-12   3.089e+00  -5.361e-02  -1.016e-02  -2.116e-14   6.616e-14
    13  -4.243e-04   5.610e+00   3.083e-02  -8.792e-03  -2.422e-03   1.753e-04
    14   4.243e-04   5.610e+00   3.083e-02  -8.792e-03   2.422e-03  -1.753e-04
    15   3.191e-12   5.610e+00  -6.175e-02  -6.303e-03   5.818e-15  -6.019e-14
  MODE     2:   f= 2.932488 Hz,  T= 0.341007 sec
		X- modal participation factor =  -2.2239e-01 
		Y- modal participation factor =   3.2940e-12 
		Z- modal participation factor =  -3.9588e-12 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     2   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     3   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     4  -4.063e-01   2.792e-01  -1.189e-02  -6.187e-03  -6.294e-03  -2.958e-03
     5  -4.063e-01  -2.792e-01   1.189e-02   6.187e-03  -6.294e-03  -2.958e-03
     6  -2.102e-01   2.024e-10  -3.264e-12  -1.518e-12  -3.066e-03  -2.934e-03
     7  -1.466e+00   1.039e+00  -2.338e-02  -8.629e-03  -8.646e-03  -1.038e-02
     8  -1.466e+00  -1.039e+00   2.338e-02   8.629e-03  -8.646e-03  -1.038e-02
     9  -7.387e-01  -2.569e-11  -6.701e-12   4.692e-12  -4.264e-03  -1.048e-02
    10  -3.159e+00   2.148e+00  -3.236e-02  -8.460e-03  -8.572e-03  -2.127e-02
    11  -3.159e+00  -2.148e+00   3.236e-02   8.460e-03  -8.572e-03  -2.127e-02
    12  -1.656e+00  -6.115e-11  -1.174e-11  -9.049e-12  -4.791e-03  -2.154e-02
    13  -5.594e+00   3.426e+00  -3.686e-02  -4.855e-03  -4.462e-03  -3.377e-02
    14  -5.594e+00  -3.426e+00   3.686e-02   4.855e-03  -4.462e-03  -3.377e-02
    15  -3.197e+00  -4.522e-11  -2.235e-11   1.440e-11  -2.859e-03  -3.429e-02
  MODE     3:   f= 3.666268 Hz,  T= 0.272757 sec
		X- modal participation factor =   1.2170e-01 
		Y- modal participation factor =   2.4238e-10 
		Z- modal participation factor =  -2.6039e-10 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     2   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     3   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     4   7.555e-02   5.519e-01   1.248e-02  -9.903e-03  -1.453e-04  -5.817e-03
     5   7.555e-02  -5.519e-01  -1.248e-02   9.903e-03  -1.453e-04  -5.817e-03
     6   4.616e-01   1.538e-08  -1.796e-10  -8.300e-11   5.806e-03  -5.737e-03
     7   2.388e-01   1.962e+00   2.417e-02  -1.294e-02  -1.306e-04  -1.963e-02
     8   2.388e-01  -1.962e+00  -2.417e-02   1.294e-02  -1.306e-04  -1.963e-02
     9   1.613e+00  -2.708e-09  -3.518e-10   2.378e-10   7.367e-03  -1.969e-02
    10   6.105e-01   3.953e+00   3.285e-02  -1.114e-02   5.413e-04  -3.920e-02
    11   6.105e-01  -3.953e+00  -3.285e-02   1.114e-02   5.413e-04  -3.920e-02
    12   3.378e+00  -5.067e-09  -6.630e-10  -2.701e-10   6.673e-03  -3.945e-02
    13   1.464e+00   6.103e+00   3.676e-02  -4.283e-03   7.178e-04  -6.036e-02
    14   1.464e+00  -6.103e+00  -3.676e-02   4.283e-03   7.178e-04  -6.036e-02
    15   5.736e+00  -5.239e-09  -1.163e-09   5.005e-10   2.699e-03  -6.077e-02
  MODE     4:   f= 5.869413 Hz,  T= 0.170375 sec
		X- modal participation factor =   6.0896e-08 
		Y- modal participation factor =   1.3164e-01 
		Z- modal participation factor =  -2.0427e-03 
  Joint   X-dsp       Y-dsp       Z-dsp       X-rot       Y-rot       Z-rot
     1   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     2   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     3   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00   0.000e+00
     4   4.640e-04   1.221e+00  -1.637e-04  -2.216e-02  -3.508e-03   5.441e-04
     5  -4.629e-04   1.221e+00  -1.637e-04  -2.216e-02   3.508e-03  -5.442e-04
     6   5.361e-06   1.221e+00   1.054e-04  -1.970e-02   4.420e-08  -7.016e-08
     7  -1.172e-03   3.674e+00  -8.748e-03  -1.575e-02  -3.265e-03   1.514e-03
     8   1.174e-03   3.674e+00  -8.748e-03  -1.575e-02   3.265e-03  -1.514e-03
     9   8.122e-06   3.679e+00   1.693e-02  -1.352e-02  -2.952e-08  -9.908e-08
    10  -1.004e-03   3.726e+00  -2.741e-02   1.721e-02   3.398e-03   1.459e-03
    11   1.002e-03   3.726e+00  -2.741e-02   1.721e-02  -3.398e-03  -1.459e-03
    12  -7.060e-06   3.727e+00   5.369e-02   1.532e-02  -2.514e-08   8.352e-08
    13   7.516e-04  -4.208e+00  -4.447e-02   3.445e-02   1.180e-02  -1.670e-03
    14  -7.524e-04  -4.208e+00  -4.447e-02   3.445e-02  -1.180e-02   1.670e-03
    15   4.185e-06  -4.209e+00   8.711e-02   2.293e-02   7.387e-09  -6.777e-08
M A T R I X    I T E R A T I O N S: 5
There are 4 modes below 5.869413 Hz. ... All 4 modes were found.
