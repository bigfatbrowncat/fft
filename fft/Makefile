CXXFLAGS =	-O3 -Wall -fmessage-length=0

OBJS =		fft_core.o fft_test.o fft_main.o

LIBS =

TARGET =	fft.exe

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
 