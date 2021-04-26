all: MissingMassNtupleAnalyzer

MissingMassNtupleAnalyzer: MissingMassNtupleAnalyzer.cc
	g++ -g --std=c++11 `root-config --libs` `root-config --cflags` \
        MissingMassNtupleAnalyzer.cc -o MissingMassNtupleAnalyzer

clean:
	rm MissingMassNtupleAnalyzer
