CXX = g++
CXXFLAGS = -g -std=gnu++11

LDFLAGS = -L/usr/X11R6/lib -lglut -lGL -lGLU -lXi -lXmu -lcholmod

SRCS = main.cc ArcBall.cc uistate.cc Vector3.cc geodesic.cc
OBJS = $(SRCS:.cc=.o)
PROG = geodist

all: $(PROG) square_gen

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

geodist_tet: main_tet.cc geodesic_tet.cc ArcBall.cc uistate.cc Vector3.cc
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

square_gen: square_gen.c
	gcc square_gen.c -o square_gen

clean:
	rm -f *.o $(PROG) square_gen *~

.PHONY: clean all
