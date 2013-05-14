CC=g++
LDFLAGS=
CFLAGS=-Wall -std=c++11
SRC=main.cpp
OBJ=$(SRC:.cpp=.o)
	BIN=program
all: $(SRC) $(BIN)
$(BIN): $(OBJ)
	$(CC) $(LDFLAGS)  $(OBJ) -o $@
.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@
launch:
	./$(BIN)

pdf:
	pdflatex article.tex
	rm *.log *.aux 2>/dev/null
