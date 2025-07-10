g++ -shared -DCMAKE_BUILD_TYPE=Debug -fPIC -o libAnalysisModule.so *.cpp `root-config --cflags --libs` -lEve
