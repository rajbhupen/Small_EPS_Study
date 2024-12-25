LIB = -L${ROOTSYS}/lib -lCore  -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint  -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lMinuit -L/usr/lib64 -lc -lm -ldl -lcrypt -lpthread 

INC = -I${ROOTSYS}/include

sipm_mass_test:             
	g++ -g -pthread -m64 -Wno-deprecated -std=c++17 -o oscill_bin oscill_bin.C $(LIB) $(INC)
