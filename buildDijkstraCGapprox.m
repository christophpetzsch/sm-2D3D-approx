% whenever it doesnt compile check if you inited and updated submodules :) 

mex dijkstraCGapprox/dijkstraCGapprox.cpp dijkstraCGapprox/MinHeapCGapprox.cpp -I"tools/eigen" ...
    CXXFLAGS='$CXXFLAGS -std=c++17'
% you might replace CXXFLAGS='$CXXFLAGS -std=c++17' with COMPFLAGS='$COMPFLAGS -std=c++17'

