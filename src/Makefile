OSYS = $(shell uname)

ifeq ($(OSYS),Darwin)
	FLAGS = -std=c++11 -I`root-config --incdir` -I$(HOME)/.local/include -O3 -g -mmacos-version-min=12.0
	CXX = clang++
endif
ifeq ($(OSYS),Linux)
	FLAGS = -std=c++11 -I`root-config --incdir` -I$(HOME)/.local/include -O3 -g
	CXX = g++
endif

LDFLAGS = -L$(HOME)/.local/lib -lpixie -lm

all : libFDSi

libFDSi : DEGAi.o LaBr.o DEGAiHit.o Event.o Gamma.o Implant.o
	$(CXX) $(FLAGS) -shared -o ../lib/libFDSi.so ../obj/DEGAi.o ../obj/LaBr.o ../obj/DEGAiHit.o ../obj/Event.o ../obj/Gamma.o ../obj/Implant.o $(LDFLAGS)

DEGAi.o : DEGAi.cc DEGAi.hh
	$(CXX) $(FLAGS) -fPIC -c -o ../obj/DEGAi.o DEGAi.cc

LaBr.o : LaBr.cc LaBr.hh
	$(CXX) $(FLAGS) -fPIC -c -o ../obj/LaBr.o LaBr.cc

DEGAiHit.o : DEGAiHit.cc DEGAiHit.hh
	$(CXX) $(FLAGS) -fPIC -c -o ../obj/DEGAiHit.o DEGAiHit.cc

Event.o : Event.cc Event.hh
	$(CXX) $(FLAGS) -fPIC -c -o ../obj/Event.o Event.cc

Gamma.o : Gamma.cc Gamma.hh
	$(CXX) $(FLAGS) -fPIC -c -o ../obj/Gamma.o Gamma.cc

Implant.o : Implant.cc Implant.hh
	$(CXX) $(FLAGS) -fPIC -c -o ../obj/Implant.o Implant.cc

clean :
	rm ../obj/*.o
	rm ../lib/libFDSi.*

install :
	cp *.hh ~/.local/include/libFDSi
	cp ../lib/libFDSi.so ~/.local/lib
