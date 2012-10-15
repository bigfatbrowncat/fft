CXXFLAGS =	-O3 -Wall -fmessage-length=0 -ffast-math

OBJS =		fft_core.o fft_test.o fft_main.o

LIBS =

TARGET =	fft.exe

$(TARGET):	$(OBJS)
	$(CXX) -static-libgcc -static-libstdc++ -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
 