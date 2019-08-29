# Makefile for compiling a the idraw
CC=g++
LD=g++
RM=rm -f
#Production flags
CFLAGS=-Wno-unused-result -std=c++11 -fopenmp
LDFLAGS=
UBUNTUFLAGS=-lm $(LDFLAGS) -lstdc++ -std=c++11 -fopenmp
OPENSUSEFLAGS=-lm $(LDFLAGS)
DEBUG=-O2
OBJS:=moea.o RealLife-MOPs.o cec09.o wfg/WFG1/WFG1.o wfg/WFG2/WFG2.o wfg/WFG3/WFG3.o wfg/WFG4/WFG4.o wfg/WFG5/WFG5.o wfg/WFG6/WFG6.o wfg/WFG7/WFG7.o wfg/WFG8/WFG8.o wfg/WFG9/WFG9.o 
MAIN=Ejecutable
all:$(MAIN)
ubuntu:
	$(MAKE) "OPENSUSEFLAGS="
opensuse:
	$(MAKE) "UBUNTUFLAGS="
$(MAIN):$(OBJS)
	$(LD) $(DEBUG) $(OPENSUSEFLAGS) $(OBJS) -o $(MAIN) $(UBUNTUFLAGS)
%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(DEBUG) -c $<
%.o: %.c %.h
	$(CC) $(CFLAGS) $(DEBUG) -c $<	
%.o: %.c
	$(CC) $(CFLAGS) $(DEBUG) -c $<	
%.o: %.cpp 
	$(CC) $(CFLAGS) $(DEBUG) -c $<	
wfg/WFG1/WFG1.o:  wfg/WFG1/WFG1.cpp
	g++ -c $< -o wfg/WFG1/WFG1.o 
wfg/WFG2/WFG2.o:  wfg/WFG2/WFG2.cpp
	g++ -c $< -o wfg/WFG2/WFG2.o 
wfg/WFG3/WFG3.o:  wfg/WFG3/WFG3.cpp
	g++ -c $< -o wfg/WFG3/WFG3.o 
wfg/WFG4/WFG4.o:  wfg/WFG4/WFG4.cpp
	g++ -c $< -o wfg/WFG4/WFG4.o 
wfg/WFG5/WFG5.o:  wfg/WFG5/WFG5.cpp
	g++ -c $< -o wfg/WFG5/WFG5.o 
wfg/WFG6/WFG6.o:  wfg/WFG6/WFG6.cpp
	g++ -c $< -o wfg/WFG6/WFG6.o 
wfg/WFG7/WFG7.o:  wfg/WFG7/WFG7.cpp
	g++ -c $< -o wfg/WFG7/WFG7.o 
wfg/WFG8/WFG8.o:  wfg/WFG8/WFG8.cpp
	g++ -c $< -o wfg/WFG8/WFG8.o 
wfg/WFG9/WFG9.o:  wfg/WFG9/WFG9.cpp
	g++ -c $< -o wfg/WFG9/WFG9.o 



#$(CC) wfg/WFG1/WFG1.h $(CFLAGS) $(DEBUG) -c -o wfg/WFG1/WFG1.o
debug:
	$(MAKE) "DEBUG=-g"
clean:
	$(RM) $(OBJS) *.*~
cleanall:
	$(RM) $(OBJS) $(MAIN) *.*~
	
