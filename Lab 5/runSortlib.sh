#!/bin/bash

export JAVA_OPTS="-Xss10m"
kotlin -cp lib/lets-plot-jfx-2.0.2.jar:lib/lets-plot-kotlin-api-2.0.1.jar:lib/lets-plot-image-export-2.0.2.jar:lib/batik-all-1.12.jar:lib/w3c.jar:lib/kotlin-logging-1.12.4.jar:lib/slf4j-api-1.7.30.jar:lib/slf4j-simple-1.7.30.jar:lib/jaxp-1.4.jar:lib/sac.jar:lib/xmlgraphics-commons-1.5.jar:lib/javafx-swing-11.jar:lib/javafx-graphics-11.jar:lib/javafx-base-11.jar:. MainKt $*
