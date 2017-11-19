CC = mpicc
CFLAGS = 
LDFLAGS = 
TARGET = xrandpfile

all: bin/$(TARGET)

obj:
	mkdir obj
bin:
	mkdir bin
output:
	mkdir output

obj/%.o: src/%.o $(filter-out $(wildcard obj), obj)
	mv src/$*.o obj/

obj/$(TARGET): obj/$(TARGET).o


bin/$(TARGET): obj/$(TARGET) $(filter-out $(wildcard bin), bin)
	cp obj/$(TARGET) bin/$(TARGET)

clean:
	rm -rf obj

purge: clean
	rm -rf bin
	rm -rf output

../01_writing_mpi_files/output/file.bin:
	$(error Error: the file ../01_writing_mpi_files/output/file.bin does not exists - run the associated make submit)

submit: all $(filter-out $(wildcard output), output) $(filter-out $(wildcard ../01_writing_mpi_files/output/file.bin), ../01_writing_mpi_files/output/file.bin ) 
	qsub job.sh
