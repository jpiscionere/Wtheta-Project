NRi=/home/piscioja/NRC/include
NRd=/home/piscioja/NRC/recipes
LIB_DIR=/data0/sys


CLINK=gcc
FLINK=f77
CFLAGS= -g03 -Wall -lm  
GSL_FLAGS= -L$(LIB_DIR)/lib -lgsl -lgslcblas -lm

NRCd=/home/piscioja/NRC/recipes
libd=/home/piscioja/lib
maskd=/home/piscioja/SDSS-Mask

all:  test covariance_fitting covariance_fitting2

OBJS15= covariance_fitting.o  $(NRCd)/ran2.o  $(NRCd)/sort2ii.o $(NRCd)/nrutil.o $(NRCd)/gaussj.o
covariance_fitting: $(OBJS15)
	$(CLINK) $(CFLAGS) $(GSL_FLAGS) -o covariance_fitting $(OBJS15)

OBJS16= covariance_fitting2.o  $(NRCd)/ran2.o  $(NRCd)/sort2ii.o $(NRCd)/nrutil.o $(NRCd)/gaussj.o
covariance_fitting2: $(OBJS16)
	$(CLINK) $(CFLAGS) $(GSL_FLAGS) -o covariance_fitting2 $(OBJS16)

OBJS17= test.o  $(NRCd)/ran2.o  $(NRCd)/sort2ii.o $(NRCd)/nrutil.o $(NRCd)/gaussj2.o
test: $(OBJS17)
	$(CLINK) $(CFLAGS) -o test $(OBJS17)



clean: 
	rm -f $(OBJS15) $(OBJS16) $(OBJS17)
