KC= kotlinc-jvm
KFLAGS= -include-runtime
SRC= Sortlib.kt Main.kt
PROG= TestSort.jar

all:
	$(KC) $(SRC) $(KFLAGS) -d $(PROG)

.PHONY: clean

clean:
	rm -rf $(PROG)