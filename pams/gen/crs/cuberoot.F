*Used to get the real cuberoot of a neg number
      Real*4 Function CubeRoot(x)
      Implicit none
      Real*4 X
      CubeRoot=(X/Abs(x))*(Abs(X)**(1.0d0/3.0d0))
      End
