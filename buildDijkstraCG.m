% whenever it doesnt compile check if you inited and updated submodules :) 

mex dijkstraCG/dijkstraCG.cpp dijkstraCG/MinHeapCG.cpp -I"dijkstraCG/eigen" ...
    CXXFLAGS='$CXXFLAGS -std=c++17'
% you might replace CXXFLAGS='$CXXFLAGS -std=c++17' with COMPFLAGS='$COMPFLAGS -std=c++17'

