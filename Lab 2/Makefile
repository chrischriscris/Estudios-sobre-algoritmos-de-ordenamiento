KC=	kotlinc
KFLAG=	-cp
LIB=	lib

all:	\
	SwingJfxDemoFrontendContext.class\
	PlotRuntimeKt.class\
	MainKt.class

SwingJfxDemoFrontendContext.class: SwingJfxDemoFrontendContext.kt
	$(KC) $(KFLAG) $(LIB)/lets-plot-jfx-2.0.2.jar:$(LIB)/lets-plot-kotlin-api-2.0.1.jar SwingJfxDemoFrontendContext.kt

PlotRuntimeKt.class: PlotRuntime.kt
	$(KC) $(KFLAG) $(LIB)/lets-plot-kotlin-api-2.0.1.jar:. PlotRuntime.kt

MainKt.class: Main.kt
	$(KC) $(KFLAG) . Sortlib.kt Main.kt

clean:
	rm -rf *.png *.class META-INF