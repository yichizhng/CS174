CXX = g++
CXXFLAGS = -g -std=gnu++11

LDFLAGS = -L/usr/X11R6/lib -lglut -lGL -lGLU -lXi -lXmu -lcholmod

SRCS = main.cc ArcBall.cc uistate.cc Vector3.cc geodesic.cc
OBJS = $(SRCS:.cc=.o)
PROG = geodist

TETSRCS = main_tet.cc geodesic_tet.cc ArcBall.cc uistate.cc Vector3.cc
TETOBJS = $(TETSRCS:.cc=.o)
TETPROG = geodist_tet

all: $(PROG) $(TETPROG) square_gen

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(TETPROG): $(TETOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

square_gen: square_gen.c
	gcc square_gen.c -o square_gen

clean:
	rm -f *.o $(PROG) square_gen *~

.PHONY: clean all
