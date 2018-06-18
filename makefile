# simple make file
SOURCES=2Ar_PES_GP_Symm.f90 PES_GP_Symm.f90
PRODUCT=2Ar_PES.out


all: $(PRODUCT)

$(PRODUCT) : $(SOURCES)
	gfortran -o $(PRODUCT) $(SOURCES)
