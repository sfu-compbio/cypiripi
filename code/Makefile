CPLEXDIR?=~/Applications/opt/cplex12_5_1
CPLEXARCH?=x86-64_sles10_4.1

CPLEXINC=$(CPLEXDIR)/cplex/include
CPLEXLIB=$(CPLEXDIR)/cplex/lib/$(CPLEXARCH)/static_pic
CONCERTINC=$(CPLEXDIR)/concert/include
CONCERTLIB=$(CPLEXDIR)/concert/lib/$(CPLEXARCH)/static_pic
CPLEXFLAGS=-L $(CPLEXLIB) -L $(CONCERTLIB)  -Wl,--start-group -lconcert -lilocplex -lcplex

CC:=g++
DF=-O3 -DNDEBUG

CFLAGS=-c $(DF) -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -std=c++11 
CFLAGS+=-DIL_STD -I $(CPLEXINC) -I $(CONCERTINC) 
LDFLAGS= -lm  $(CPLEXLIB)/libilocplex.a $(CPLEXLIB)/libcplex.a $(CONCERTLIB)/libconcert.a -lpthread
#LDFLAGS= -Wl,-Bstatic -libc++ -Wl,-Bdynamic $(LF) -lm -lpthread $(CPLEXFLAGS) 

OBJECTS=$(SOURCES:.cc=.o)
SOURCES=$(wildcard *.cc)
EXECUTABLE=allele

all: $(SOURCES) $(EXECUTABLE) 

debug: DF=-g
debug: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(EXECUTABLE)

# allele2:
# 	$(CC) allele2.cc $(CFLAGS) -DIL_STD -I $(CPLEXINC) -I $(CONCERTINC) $(LDFLAGS) $(CPLEXFLAGS) -o allele2

# allele2a:
# 	$(CC) allele2a.cc $(CFLAGS)  $(LDFLAGS) -o allele2a

# alleleNew:
# 	$(CC) alleleNew.cc $(CFLAGS) -DIL_STD -I $(CPLEXINC) -I $(CONCERTINC) $(LDFLAGS) $(CPLEXFLAGS) -o alleleNew
# alleleNewG:
# 	$(CC) alleleNew.cc -g -std=c++11 -DIL_STD -I $(CPLEXINC) -I $(CONCERTINC) $(LDFLAGS) $(CPLEXFLAGS) -o alleleNew


# plot:
# 	$(CC) plot.cc $(CFLAGS) $(LDFLAGS) -o plot

# clean:
# 	rm -f allele plot
