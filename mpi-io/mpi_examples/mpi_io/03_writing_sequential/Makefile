CC = mpicc
CFLAGS = 
LDFLAGS = 
TARGET = mkdrandfiles

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

submit: all $(filter-out $(wildcard output), output)
	qsub job.sh
