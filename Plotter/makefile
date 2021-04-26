all: Plotter EventList

Plotter: Plotter.cc 
	g++ -I`root-config --incdir` `root-config --libs` `root-config --cflags` \
	Plotter.cc -o Plotter

EventList: EventList.cc 
	g++ -I`root-config --incdir` `root-config --libs` `root-config --cflags` \
        EventList.cc -o EventList

clean:
	rm Plotter EventList
