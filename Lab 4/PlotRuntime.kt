/**
 * Copyright (c) 2021. Universidad Simón Bolívar
 *
 * Este código es publicado bajo la MIT license    
 *
 * Descripción: Librería para la generación de gráficas de
 * rendimiento de algoritmos. El objetivo es el
 * poder comparar el tiempo de ejecución de varios 
 * algoritmo, con respecto a varios tamaños de entrada
 * 
 * @author Guillermo Palma <gvpalma@usb.ve>
 * @since 0.1 
 */

import jetbrains.letsPlot.LetsPlot
import jetbrains.letsPlot.geom.geomPoint
import jetbrains.letsPlot.geom.geomLine
import jetbrains.letsPlot.geom.geomErrorBar
import jetbrains.letsPlot.letsPlot
import jetbrains.letsPlot.export.ggsave
import jetbrains.letsPlot.label.ggtitle
import jetbrains.letsPlot.label.xlab
import jetbrains.letsPlot.label.ylab
import jetbrains.letsPlot.scale.scaleColorBrewer
import SwingJfxDemoFrontendContext

/**
 * Genera un gráfico que muestra los trazados
 * correspondientes a los tiempos de ejecución
 * de varios algoritmos con respecto a varias
 * entradas de diverso tamaño. Esta función
 * muestra un panel con las gráficas del  
 * comportamiento de los algoritmos. También 
 * crea una imagen que salva en archivo .png de
 * la gráfica mostrada en el panel. Se supone que cada
 * entrada es ejecutada más de una vez por un algoritmo.
 * Por ello se solicita el tiempo promedio, el menor tiempo 
 * y el mayor tiempo de ejecución. 
 *
 * Por ejemplo, suponga que se corren dos algoritmos
 * en tres tamaños de entrada. El formato del arreglo 
 * con las etiquetas de los algoritmos es: 
 * ["Algo1", "Algo1", "Algo1", "Algo2", "Algo2", "Algo2"]
 * El formato de los tamaños de la entradas es:
 * [n1, n2, n3, n1, n2, n3]
 * El formato del arreglo con  los tiempos de ejecución es:
 * [t_n1_Alg1, t_n2_Alg1, t_n3_Alg1, t_n1_Alg2, t_n2_Alg2, t_n3_Alg2]
 *
 * @param[windowTitle] Título del panel
 * @param[imgPath] Camino a donde va a ser guardada la imagen
 * @param[imgName] Nombre de la imagen con la gráfica
 * @param[title] Título de la gráfica
 * @param[xLabel] Etitqueta del eje X
 * @param[yLabel] Etiquetea del eje Y
 * @param[algorithmsLabels] Etiqueta de los algoritmos
 * @param[numElements] Número de elementos de entrada
 * @param[minTimes] Menores tiempos obtenidos por los algoritmos
 * @param[averageTimes] Tiempos promedios obtenidos por los algoritmos
 * @param[maxTimes] Tiempos máximos obtenidos por los algoritmos
 */
fun plotRuntime(windowTitle: String, 
		imgPath: String,
		imgName: String,
		title: String,
		xLabel: String,
		yLabel: String,
		algorithmsLabels: Array<String>,
		numElements: Array<Int>,
		minTimes: Array<Double>,
		averageTimes: Array<Double>,
		maxTimes: Array<Double>) {
    val ctx = SwingJfxDemoFrontendContext(title=windowTitle)
    LetsPlot.frontendContext = ctx
    val data = mapOf(
	"Algoritmo" to algorithmsLabels, 
	"n" to numElements,
	"tiempo_avg" to averageTimes,
	"tiempo_min" to minTimes,
	"tiempo_max" to maxTimes
    )
    val p = letsPlot(data) {x="n"; color="Algoritmo"} +
    geomErrorBar(width=.1) {ymin="tiempo_min"; ymax="tiempo_max"} +
    geomLine {y="tiempo_avg"} + geomPoint {y="tiempo_avg"} +
    ggtitle(title) +
    xlab(xLabel) +
    ylab(yLabel) +
    scaleColorBrewer(palette = "Dark2")
    p.show()
    ggsave(p, imgName, path=imgPath)
    ctx.showAll()
}
